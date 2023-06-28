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
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main() {  

  cout << "Have you cleaned the .out files? [y/n]" << endl;
  char answer;
  cin >> answer;
  if (answer == 'y' ) { 
  }else if (answer == 'n') { 
    return 0; 
  } 
/* 
  cout << "Which temperature of the Ising Model would you like to simulate? Write it down as for example 1.3: ";
  cin >> temperatura_str;
  temp = stod(temperatura_str);
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;
*/

  Input(); //Inizialization
  for(double T = tmin; T<tmax; T=T+delta_T) { //ciclo su tutte le temperature da voler studiare
    temperatura_str= to_string(T);
    temperatura_str.resize(temperatura_str.size() - 5); //cosi' da avere i numeri in formato 2.0
    beta = 1.0/T;
    cout << "Temperature = " << temperatura_str << endl;
    if(restart!=1){ //se non ripartiamo da una config di equilibrio, termalizziamo il sistema facendo 2000 passi
      //cout << "Termalizzo..." <<endl;
      for (int i = 0; i< 2000; i++) Move(metro);
    }


    for(int iblk=1; iblk <= nblk; ++iblk) { //Simulation 
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep) {
        Move(metro);
        Measure();
        Accumulate(); //Update block averages
      }
      Averages(iblk);   //Print results for current block
    }
    Averages_final(T);
    ConfFinal(); //Write final configuration
  }
  return 0;
}


void Input(void) {
  ifstream ReadInput, ReadConf, ReadVelocity,Seed;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  //Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  
  //Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> restart;
  if(restart) Seed.open("seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


  //Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

  //initial configuration
  if(restart) {
    ReadConf.open("config.final");
    for (int i=0; i<nspin; ++i) {
      ReadConf >> s[i];
    }
  }
  else {
    for (int i=0; i<nspin; ++i) {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = -1;
    }
  }
  ReadConf.close();
  
  //Evaluate energy etc. of the initial configuration
  Measure();

  //Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro) {
  int o;
  double probability;
  //double energy_up, energy_down, energy_old, energy_new, sm;
  double A, delta_E;
  int s_gibbs;

  for(int i=0; i<nspin; ++i) {
    //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) { //Metropolis
      attempted ++;
      delta_E = 2.0*Boltzmann(-1*s[o], o);
      A = min(1.0, exp(-beta*delta_E));
      if (A == 1) {
        s[o] = -1*s[o];
        accepted++;
      } 
      else {
        if(rnd.Rannyu() < A){
          s[o] = -1*s[o];
          accepted++;
        }
      }
    }
    else { //Gibbs sampling
      accepted =1;
      attempted=1; //cause we accept every single move therefore there is no acceptance rate
      if(rnd.Rannyu() >= 0.5) s_gibbs=1;
      else s_gibbs =-1;
      //s_gibbs e' 1 oppure -1

      delta_E = -2.0*Boltzmann(s_gibbs, o); //- per neutralizzare il meno di boltzmann
      probability = 1./(1.+exp(-beta*delta_E));
      if(rnd.Rannyu() < probability) {
        s[o] = s_gibbs;
      } 
      else{
        s[o] = -1*s_gibbs;
      } //accettiamo qualunque mossa proposta
    }
  }
}

double Boltzmann(int sm, int ip) { //calcola il DeltaE di energia
//sm e' lo spin che considero, s[o], mentre ip e' proprio o quindi la sua posizione nell'array. 
//Ricordati che calcoli il deltaE ipotizzando di aver flippato lo spin quindi sm e' -1*spin considerato

  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure() {
  int bin;
  double u = 0.0, m = 0.0;
  //double u_sqr =0., m_sqr =0.;

  //cycle over spins
  for (int i=0; i<nspin; ++i) {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
	  //u_sqr += u*u;

	  m += s[i];
	  //m_sqr += m*m;
  }
  walker[iu] = u; //questa non e' ancora la media di H, quello si fa in avareges quindi la riga sotto e' sbagliata
  //walker [ic] = beta * beta * (u_sqr - u*u/nspin) /nspin;
  walker [ic] = u*u;
  walker [im] = m;
  walker [ix] =m*m ;
}


void Reset(int iblk) { //Reset block averages
  if(iblk == 1){
    for(int i=0; i<n_props; ++i) {
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for(int i=0; i<n_props; ++i) {
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}


void Accumulate(void) { //Update block averages
  for(int i=0; i<n_props; ++i) {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) { //Print results for current block
//in output ho i risultati per tutte le temperature della media a blocchi delle quantita' interessanti ad ogni blocco. 
//Stampo blocco, media nel blocco, media progressiva e suo errore
    
  ofstream ene, heat, mag, chi;
  const int wd=12;
    
  //cout << "Block number " << iblk << endl;
  //cout << "Acceptance rate " << accepted/attempted << endl << endl;
  if (h==0.0) {
    ene.open("output/" + temperatura_str + "/output.ene.out",ios::app); //temperatura_str e' una stringa ma sara' poi un numero come nome del file
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy: this is <H>/nspin
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    ene << " " << iblk <<  " " << stima_u << " " << glob_av[iu]/(double)iblk << " " << err_u << endl;
    ene.close();

    heat.open("output/" + temperatura_str + "/output.heat.out",ios::app); //Heat capacity
    stima_c = beta*beta* (blk_av[ic]/blk_norm/(double)nspin - stima_u*stima_u*(double)nspin); //o analogamente beta*beta* (blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm,2))/(double)nspin
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    heat << " " << iblk <<  " " << stima_c << " " << glob_av[ic]/(double)iblk << " " << err_c << endl;
    heat.close();

    chi.open("output/" + temperatura_str + "/output.chi.out",ios::app);
    stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //Magnetic susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    chi << " " << iblk <<  " " << stima_x << " " << glob_av[ix]/(double)iblk << " " << err_x << endl;
    chi.close();
  }
  if (h!=0) { 
    mag.open("output/" + temperatura_str + "/output.mag.out",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    mag << " " << iblk <<  " " << stima_m << " " << glob_av[im]/(double)iblk << " " << err_m << endl;
    mag.close();
  }
  //cout << "----------------------------" << endl << endl;
}

void Averages_final(double temp) { //Print results for the final block
//stampo i risulati migliori, quelli dell'ultimo blocco, ad ogni temperatura. 
//Stampo temperatura, media progressiva e errore associato

  ofstream ene, heat, mag, chi;
  if (metro==1) {
    if(h==0) {
      ene.open("Risultati_finali/Metropolis/ene_metro.out",ios::app);
      ene << temp << " " << glob_av[iu]/(double)nblk << " " << err_u << endl;
      ene.close();

      heat.open("Risultati_finali/Metropolis/heat_metro.out",ios::app);
      heat << temp << " " << glob_av[ic]/(double)nblk << " " << err_c << endl;
      heat.close();

      chi.open("Risultati_finali/Metropolis/chi_metro.out",ios::app);
      chi << temp << " " << glob_av[ix]/(double)nblk << " " << err_x << endl;
      chi.close();

    }else if(h!=0) {
      mag.open("Risultati_finali/Metropolis/mag_metro.out",ios::app);
      mag << temp << " " << glob_av[im]/(double)nblk << " " << err_m << endl;
      mag.close();
    }
  } else if (metro==0) { //gibbs
    if(h==0) {
      ene.open("Risultati_finali/Gibbs/ene_gibbs.out",ios::app);
      ene << temp << " " << glob_av[iu]/(double)nblk << " " << err_u << endl;
      ene.close();

      heat.open("Risultati_finali/Gibbs/heat_gibbs.out",ios::app);
      heat << temp << " " << glob_av[ic]/(double)nblk << " " << err_c << endl;
      heat.close();

      chi.open("Risultati_finali/Gibbs/chi_gibbs.out",ios::app);
      chi << temp << " " << glob_av[ix]/(double)nblk << " " << err_x << endl;
      chi.close();
    }else if(h!=0) {
      mag.open("Risultati_finali/Gibbs/mag_gibbs.out",ios::app);
      mag << temp << " " << glob_av[im]/(double)nblk << " " << err_m << endl;
      mag.close();
    }

  }

}


void ConfFinal(void) {
  ofstream WriteConf;

  //cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("output/" + temperatura_str + "/config.final");
  for (int i=0; i<nspin; ++i){
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i) {  //Algorithm for periodic boundary conditions

    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk) {
  if(iblk==1) return 0.0;
  else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
