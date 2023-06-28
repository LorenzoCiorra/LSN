/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <string>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "MD_MC.h"

using namespace std;

int main() {

  cout << "Have you cleaned the .dat files? [y/n]" << endl;
  char answer;
  cin >> answer;
  if (answer == 'y' ) { 
  }else if (answer == 'n') { 
    return 0; 
  } 
  cout << "Which phase of the argon would you like to simulate? Write it down in Italian choosing from solido, liquido or gas: ";
  cin >> phase;

  Input(); //Inizialization
  int nconf = 1;
  for(int j=0; j<2000; j++) Move();
  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
      if(istep%10 == 0){
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  Print_gfunction();
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

//Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

//Read input informations
  ReadInput.open(phase +"/input.in");

  ReadInput >> iNVET;
  ReadInput >> restart;

  if(restart) Seed.open(phase + "/equilibrio/seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Calcoliamo le Tail Corrections visto che sono costanti e andranno solo aggiunte
  vtail = 8.0*M_PI*rho *(1.0/(9.0*pow(rcut,9)) - 1.0/(3.0*pow(rcut,3)));
  ptail = 32.0*M_PI*rho*rho*(1.0/(9.0*pow(rcut,9)) - 1.0/(6.0*pow(rcut,3)));

//Prepare arrays for measurements
  iv = 0; //Potential energy
  it = 1; //Temperature
  ik = 2; //Kinetic energy
  ie = 3; //Total energy
  ipress = 4; //Pression of the gas 
  n_props = 5; //Number of observables

  bin_size = (box/2.0)/(double)nbins; //perche' oltre L/2 non mi interessa andare, tanto copro tutto la scatola cosi'


//Read initial configuration
  cout << "Read initial configuration" << endl << endl;
  if(restart)
  {
    ReadConf.open(phase + "/equilibrio/config.out");
    ReadVelocity.open(phase + "/equilibrio/velocity.out");
    for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
  }
  else 
  {
    ReadConf.open("config.in");
    cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i)
    {
      vx[i] = rnd.Gauss(0.,sqrt(temp));
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i)
    {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    cout << "velocity scale factor: " << fs << endl << endl;
    for (int i=0; i<npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();

  for (int i=0; i<npart; ++i)
  {
    if(iNVET)
    {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else
    {
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  
//Evaluate properties of the initial configuration
  Measure();

//Print initial values for measured properties
  cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  cout << "Initial temperature      = " << walker[it] << endl;
  cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << endl;
  cout << "Initial total energy     = " << walker[ie]/(double)npart << endl;

  return;
}


void Move()
{
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) // Monte Carlo (NVT) move
  {
    for(int i=0; i<npart; ++i)
    {
    //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
      o = (int)(rnd.Rannyu()*npart);

    //Old
      energy_old = Boltzmann(x[o],y[o],z[o],o);

    //New
      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

      energy_new = Boltzmann(x[o],y[o],z[o],o);

    //Metropolis test
      p = exp(beta*(energy_old-energy_new));
      if(p >= rnd.Rannyu())  
      {
      //Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
        accepted = accepted + 1.0;
      } else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted = attempted + 1.0;
    }
  } else // Molecular Dynamics (NVE) move
  {
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i){ //Force acting on particle i lungo la coordinata (x, y, z) = (0, 1, 2)
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    for(int i=0; i<npart; ++i){ //Verlet integration scheme

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0;
      attempted = attempted + 1.0;
    }
  }
  return;
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
 
  return f;
}

void Measure() //Properties measurement
{
  double v = 0.0, kin=0.0, p = 0.0;
  double vij, pij;
  double dx, dy, dz, dr, r; //r is the left limit of the bins

//need to be sure that all the bins are empty before starting to fill them up
for (int j =0; j<nbins; j++) g_function[j] =0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {
// distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr); //this's the distance between two particles

      if(dr < rcut)
      {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        v += vij;
		pij = 1.0/pow(dr,12) - 0.5/pow(dr,6); //sigma e epsilon si semplificano passando alle unita` naturali(di lennard-jones)
		p += pij;
      }
      for (int k =0; k<nbins; k++) {//cycle over every bin to update the histogram
        r = bin_size * k;
        if(r<dr && dr< r+bin_size) g_function[k] +=2;
      }
    }          
  }

  for (int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  walker[iv] = 4.0 * v; // Potential energy //pot_tail qui come *(v+v_tail)?
  walker[ik] = kin; // Kinetic energy
  walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature
  walker[ie] = 4.0 * v + kin;  // Total energy;
  walker[ipress] = rho*walker[it] + (16./vol) * p; //Total pressure 
  //those are all 
  return;
}


void Reset(int iblk) { //Reset block averages
   
  if(iblk == 1) {
    for(int i=0; i<n_props; ++i){
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
    for(int j=0; j<nbins; j++){
      glob_av_g[j] = 0;
      glob_av2_g[j] = 0;
    }    
  }

  for(int i=0; i<n_props; ++i) {
    blk_av[i] = 0;
  }
  for(int j=0; j<nbins; j++) {
    blk_av_g[j] = 0;
  }  
  blk_norm = 0;
  blk_norm_g=0;
  attempted = 0;
  accepted = 0;
}


void Accumulate(void) { //Update block averages, also for the g function

  for(int i=0; i<n_props; ++i){
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
  for(int j=0; j<nbins; j++) {
    blk_av_g[j] = blk_av_g[j] + g_function[j];
  }
  blk_norm_g = blk_norm_g +1.0;
}


void Averages(int iblk) { //Print results for current block
    
  ofstream Epot, Ekin, Etot, Temp, Press;
  const int wd=12;
  double r=0.0, delta_v=0.0;
    
  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted/attempted *100 << endl << endl;
    
  Epot.open(phase + "/output_epot.dat",ios::app);
  Ekin.open(phase + "/output_ekin.dat",ios::app);
  Temp.open(phase + "/output_temp.dat",ios::app);
  Etot.open(phase + "/output_etot.dat",ios::app);
  Press.open (phase + "/output_pressure.dat",ios::app);
    
  stima_pot = blk_av[iv]/blk_norm/(double)npart +vtail; //Potential energy //oppure pot_tail qui sommato al calcolo?
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot*stima_pot;
  err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
  stima_kin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy
  glob_av[ik] += stima_kin;
  glob_av2[ik] += stima_kin*stima_kin;
  err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

  stima_etot = blk_av[ie]/blk_norm/(double)npart; //Total energy
  glob_av[ie] += stima_etot;
  glob_av2[ie] += stima_etot*stima_etot;
  err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

  stima_temp = blk_av[it]/blk_norm; //Temperature
  glob_av[it] += stima_temp;
  glob_av2[it] += stima_temp*stima_temp;
  err_temp=Error(glob_av[it],glob_av2[it],iblk);

	//Pressure
  stima_press = blk_av[ipress]/blk_norm +ptail; //Stima nel blocco
  glob_av[ipress] += stima_press;  //Stima progressiva all`avanzare dei blocchi
  glob_av2[ipress] += stima_press*stima_press;
  err_press=Error(glob_av[ipress],glob_av2[ipress],iblk);
	
  for(int j = 0; j <nbins; j++) {
	r = bin_size * j;
	delta_v= (4./3.)*M_PI*(pow(r+bin_size,3)-pow(r,3));
	stima_g[j] = blk_av_g[j] /(blk_norm_g*rho*npart*delta_v);
	glob_av_g[j] += stima_g[j];
    glob_av2_g[j] += stima_g[j]*stima_g[j];
    err_g[j]=Error(glob_av_g[j],glob_av2_g[j],iblk);
    //stampo average value di g di ogni bin nel blocco, media progressiva di ogni bin e errore associato
	//G_out << stima_g[j] << " " << glob_av[j]/(double)iblk << " " << err_g[j] << endl; //se uso 30 blocchi, sulla prima riga avro 300 stime per il primo blocco 
    //perche' 100 stime della media di ogni bin, 100 medie progressive per ogni bin e 100 errori associati
  }


//Potential energy per particle
  //Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;        
  Epot << " " << iblk <<  " " << stima_pot << " " << glob_av[iv]/(double)iblk << " " << err_pot << endl; //c'e' uno spazio iniziale quindi queste sono 5 colonne in totale
//Kinetic energy
  Ekin << " " << iblk <<  " " << stima_kin << " " << glob_av[ik]/(double)iblk << " " << err_kin << endl;
//Total energy
  Etot << " " << iblk <<  " " << stima_etot << " " << setprecision(8) << glob_av[ie]/(double)iblk << " " << setprecision(8) << err_etot << endl;
//Temperature
  Temp << " " << iblk <<  " " << stima_temp << " " << glob_av[it]/(double)iblk << " " << err_temp << endl;
//Pressure
  Press << " " << iblk <<  " " << stima_press << " " << glob_av[ipress]/(double)iblk << " " << err_press << endl;

  cout << "----------------------------" << endl << endl;

  Epot.close();
  Ekin.close();
  Etot.close();
  Temp.close();
  Press.close();
}

void Print_gfunction() {
  ofstream G_out;
  double r=0.0;
  G_out.open (phase + "/output_gfunction.dat",ios::app);
  cout << "Printing final average of the g function to file output_gfunction.dat" << endl << endl;

  //scrivo la media progressiva del data blocking per l'ultimo blocco
  for(int j = 0; j <nbins; j++) {
	r = bin_size * j;
	G_out << r << " " << glob_av_g[j]/(double)nblk << " " << err_g[j] << endl; 
  }
  G_out.close();
}


void ConfFinal(void)
{
  ofstream WriteConf, WriteVelocity, WriteSeed;

  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open(phase + "/config.out");
  WriteVelocity.open(phase +"/velocity.out");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk); //fabs valore assoluto di un double 
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
