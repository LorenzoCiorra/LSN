#include <iostream>
#include <string>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "main.h"

using namespace std;

int main() {
	input();
	cout << "Do you want to proceed with a simulated annealing simulation? [y/n]" << endl;
  	cin >> answer;
  	if (answer == 'y' ) { 
  	}else if (answer == 'n') { 
    	cout << "Doing a simple simulation with mu = "<< mu << " and sigma = " << sigma << endl;
		H_average(mu, sigma);
		return 0;
  	} 

//Simulated Annealing:
//voglio ottimizzare la mu e la sigma quindi piano piano alzo sempre di piu la beta finche mu e sigma non si stabilizzano
//ad ogni T faccio 100 steps: ad ognuno calcolo H_average con nuovi mu e sigma estratti e col metro valuto se i due nuovi parametri sono buoni o no. 

	ofstream SA_out, PerStep_out;
	SA_out.open("SimulatedAnnealing_results.out");
	PerStep_out.open("PerStep_results.out");
	while(beta<beta_max) { 
		for(int i=0; i<beta_steps; i++) { //fisso beta e faccio tanti beta_steps, poi su file salvo il risultato all'ultimo beta_step per ogni beta
			H=(double)sum_energy/blocks; //stima di H all'ultimo blocco precedente
			mu_new = mu + rnd.Rannyu(-L_mu, L_mu)/sqrt(beta); //voglio spostarmi sempre meno man mano che alzo la T perche' mi aspetto che piu' la alzo, piu' mi avvicino ai parametri ottimali
			sigma_new = sigma + rnd.Rannyu(-L_sigma, L_sigma)/sqrt(beta);
			H_average(mu_new, sigma_new); //qui stampo gli energy_results ma non mi interessano per il SA
			attempted_beta ++;

			H_new =(double)sum_energy/blocks; //la nuova ultima stima di H all'ultimo blocco
			H_error = error(sum_energy,sum2_energy, blocks);
			double delta_H = H_new - H;

			double H_acceptance = min(1., exp(-beta*delta_H));
			if(H_acceptance==1) {
				H = H_new;
				mu = mu_new;
				sigma = sigma_new;
				accepted_beta ++;
			} else{
				if(rnd.Rannyu () < H_acceptance) {
					H = H_new;
					mu = mu_new;
					sigma = sigma_new;
					accepted_beta ++;
				}
			}
			PerStep_out << progressive_steps << " " << mu << " " << sigma << " " << H << endl; //se voglio plottare H, mu, sigma man mano che evolvo ad ogni beta_step
			progressive_steps++;
		}
		cout << "Acceptance ratio for the Annealing Metropolis at beta = " << beta <<" is: " << accepted_beta/attempted_beta *100 << endl;
		SA_out << beta << " " << mu << " " << sigma << " " << H << " " << H_error << endl; //stampo mu, sigma, H e err alla fine dei beta_steps per ogni beta
		beta = beta+delta_beta;
	}
	SA_out.close();
	PerStep_out.close();
	cout << "Minimo di H = " << H << " +/- " << H_error << endl;
    cout << "realizzato con mu = " << mu <<   " e sigma = " << sigma << endl;
	return 0;
}



void input() {
	int seed[4];
   	int p1, p2;
   	ifstream Primes("randomgen/Primes");
   	if (Primes.is_open()){
       Primes >> p1 >> p2 ;
   	} else cerr << "PROBLEM: Unable to open Primes" << endl;
   	Primes.close();

   	ifstream input("randomgen/seed.in");
   	string property;
   	if (input.is_open()){
      	while ( !input.eof() ){
         input >> property;
        if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
        	}
    	}
    input.close();
   	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	rnd.SaveSeed();


	ifstream ReadInput;
    ReadInput.open("input.in");
	ReadInput >> mu;
	cout << "Mu = " << mu << endl;
	ReadInput >> sigma;
	cout << "Sigma = " << sigma << endl;
	ReadInput >> x;
	cout << "Starting point for the simulation = " << x << endl;

	ReadInput >> L;
	ReadInput >> blocks;
	ReadInput >> M;
	steps=int(M/blocks);
	cout << "The program perform Metropolis moves with uniform translations" << endl;
	cout << "Moves parameter = " << L << endl;
	cout << "Number of blocks = " << blocks << endl;
	cout << "Number of steps in one block = " << steps << endl << endl;

	ReadInput >> beta;
	ReadInput >> delta_beta;
	ReadInput >> beta_steps;
	ReadInput >> beta_max;
	ReadInput >> L_mu;
	ReadInput >> L_sigma;

	cout << "With starting Beta = " << beta << endl;
	cout << "Moves parameter for Beta steps = " << beta_steps << endl;
	cout << "And maximum beta = " << beta_max << endl;

	ReadInput.close();

}

double psi_test (double x, double mu, double s) { //ansatz variazionale
	return exp(-((x-mu)*(x-mu))/(2.*s*s)) + exp(-((x+mu)*(x+mu))/(2.*s*s));  
}

double psi_test_2derivative (double x, double mu, double s) { //derivata seconda dell'ansatz necessaria per il calcolo del valore di aspettazione dell'hamiltoniana
	double psi_1 = exp(-(pow(mu-x,2))/(2*s*s))/pow(s,4);
	double psi_2 = exp(-(pow(mu+x,2))/(2*s*s))/pow(s,4);
	return psi_1*(mu*mu-2.*mu*x-s*s+x*x) + psi_2*(mu*mu+2.*mu*x-s*s+x*x);
}

double epot (double x) {
	return pow(x,4) - 2.5*x*x;
}

double energy (double x, double mu, double sigma) {
	return epot(x) - 0.5*psi_test_2derivative(x, mu, sigma)/psi_test(x, mu, sigma); 
}

double acceptance(double mu, double sigma, double x_try, double x_old) {
	return min(1., pow(psi_test(x_try, mu, sigma),2)/pow(psi_test(x_old, mu, sigma),2) );
}

void H_average(double mu, double sigma) {
	ofstream energy_results;
   	energy_results.open("energy_results.out");
	sum_energy=0.;
	sum2_energy=0.;
	for (int i =0; i<blocks; i++) { //cycle on the blocks	
		block_ave =0.;
		attempted =0;
		accepted =0;
		for(int j=0; j<steps; j++) { //cycle on the throws in each block
			attempted ++;
			x_try = x + rnd.Rannyu(-L/2., L/2.);
			double A = acceptance (mu, sigma, x_try, x);
			if(A==1) {
				x = x_try;
				accepted ++;
			} else{
				if(rnd.Rannyu () < A) {
					x = x_try;
					accepted ++;
				}
			}
			block_ave += energy(x, mu, sigma);
		}
		//cout << attempted << endl;
		//cout << accepted << endl;
		//cout << "Block: " << i+1 << " Acceptance ratio: " << accepted/attempted *100 << endl;

		block_ave = (double)(block_ave/steps);
		sum_energy += block_ave;
		sum2_energy += block_ave*block_ave;
		
		// saving on an external file result of the block average, progressive average and its error for the expected value of H
        energy_results << i+1 << " " << block_ave << " " << sum_energy/(double)(i+1) << " " << error(sum_energy,sum2_energy, i+1) <<endl;
	}
	energy_results.close();
}

double error(double sum, double sum2, int iblk) {
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk); //fabs valore assoluto di un double 
}

