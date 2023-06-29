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
	H_average();

	return 0;
}


void input() {
    //checks for the functioning of the random numbers generator
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

    //reads all the input form the input.in file
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

double energy (double x, double mu, double sigma) { //sarebbe V+K ma K ha -0.5 come fattore moltiplicativo
	return epot(x) - 0.5*psi_test_2derivative(x, mu, sigma)/psi_test(x, mu, sigma); 
}

//Col metropolis vogliamo campionare il mod quadro della psi_test con prob di trasnzione uniforme: 
//quindi l'acceptance diventa il rapporto tra il mod quadro di psi_test calcolato in x_try fratto calcolata in x_old
double acceptance(double mu, double sigma, double x_try, double x_old) {
	return min(1., pow(psi_test(x_try, mu, sigma),2)/pow(psi_test(x_old, mu, sigma),2) );
}

void H_average() {
	ofstream energy_results, pos_results;
   	energy_results.open("energy_results.out");
	pos_results.open("position_results.out");
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
			pos_results << x << endl; //stampo l'evoluzione delle configurazioni campionate, ossia delle posizioni 
		}
		//cout << attempted << endl;
		//cout << accepted << endl;
		cout << "Block: " << i+1 << " Acceptance ratio: " << accepted/attempted *100 << endl;

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

