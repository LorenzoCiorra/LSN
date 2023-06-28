#include "randomgen/random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

double error(double av, double av_squares, int n);

int main(int argc, char *argv[]) {

  // checks for the functioning of the random numbers generator
	Random rnd;
  	int seed[4];
  	int p1, p2;
	ifstream Primes("randomgen/Primes");
	if (Primes.is_open()) {
		Primes >> p1 >> p2;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	
	ifstream input("randomgen/seed.in");
	string property;
	if (input.is_open()) {
		while (!input.eof()) {
	    	input >> property;
	      	if (property == "RANDOMSEED") {
	        	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	        	rnd.SetRandom(seed, p1, p2);
	   		}
		}
	    input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	rnd.SaveSeed();
	
	ofstream discrete_results;
	discrete_results.open("discrete_results.out");
	ofstream continous_results;
	continous_results.open("continous_results.out");

	double a =1.; //lattice spacing
	int nr_walks = 10000;
  	int steps = 100;
  	int blocks = 100;
  	int walks_in_each_block = (int)(nr_walks / blocks); // we have 100 walsks in each block
  	double mean_distance_sqr_discrete[steps][blocks]; // matrix 100x100: in each coloumn we have all 100 steps for each block
	double mean_distance_sqr_continous[steps][blocks];
	
	double sum_discrete=0;
	double error_discrete=0;
	double sum_continous=0;
	double error_continous=0;

	for(int i=0; i<steps; i++) {
		for(int j=0; j<blocks; j++){
			mean_distance_sqr_discrete[i][j] =0;
			mean_distance_sqr_continous[i][j] =0;
		}
	}

  	for (int i = 0; i < blocks; i++) { // cycle over the blocks

    	for (int j = 0; j < walks_in_each_block;j++) { // cycle over each walks in each block
      		double discrete_position[3] = {0., 0., 0.};
			double continous_position[3] = {0., 0., 0.};
      		for (int k = 0; k < steps; k++) {
                //Discrete rw
        		int direction = (int)(rnd.Rannyu(0., 3.)); // 0=x, 1=y, 2=z
        		int back_forth = (int)(rnd.Rannyu(1., 3.)); // 1=back, 2=forth; 0 is not good because then we want to sum all the outcomes and 0 would mean no movement
       			if (back_forth == 1) {
          			back_forth = -1;
        		} else if (back_forth == 2) {
         		back_forth = 1;
        		}
        		discrete_position[direction] += back_forth;

                //Continous rw
				double theta = rnd.Rannyu(0., M_PI);
				double phi = rnd.Rannyu(0., 2*M_PI); //now we calculate x,y,z in spherical coordinates
				//ricordati di sommare! La poszione si evolve nel tempo partendo dal punto precedente
				continous_position[0] += a * sin(theta) * cos(phi); //x
				continous_position[1] += a * sin(theta) * sin(phi); //y
				continous_position[2] += a * cos(theta);

		        // now we have generated a direction and a movement in that direction, we just need to calculate the distance from the origin.
				//Be careful: we need the mean distance in each block therefore we need to divide by the number of walks in each block
        		mean_distance_sqr_discrete[k][i] += (double)((pow(discrete_position[0], 2) + pow(discrete_position[1], 2) + pow(discrete_position[2], 2)) /walks_in_each_block); //attenzione: anche qui += peche' mi "allontano" sempre di piu
				mean_distance_sqr_continous[k][i] += (double)((pow(continous_position[0], 2) + pow(continous_position[1], 2) + pow(continous_position[2], 2)) /walks_in_each_block);
      		}
    	}
  	} // now we have the matrix for the mean distances of each block and each step
	//we want to calculate the total mean distance and write them down as a function of the steps. Therefore, being fixed i index of the step, we make an average over the blocks
	//praticamente vogliamo la media di tutti i blocchi dopo uno step, poi la media dopo il secondo step e cosi' via infatti poi plottiamo in funzione degli step, non del nr di blocchi

	for(int i =0; i<steps; i++) {
		sum_discrete =0; //reinitilize to 0 because otherwise we would sum values for different steps
		error_discrete=0;
		sum_continous=0;
		error_continous=0;
		for(int j=0; j<blocks; j++) {
			mean_distance_sqr_discrete[i][j] = (double)(mean_distance_sqr_discrete[i][j]/blocks);
			sum_discrete += mean_distance_sqr_discrete[i][j];
			error_discrete += blocks*pow(mean_distance_sqr_discrete[i][j],2); //se non moltiplico per blocks, ho error=nan

			mean_distance_sqr_continous[i][j] = (double)(mean_distance_sqr_continous[i][j]/blocks);
			sum_continous += mean_distance_sqr_continous[i][j];
			error_continous += blocks*pow(mean_distance_sqr_continous[i][j],2);
		}
		if(i==0) {
			discrete_results << sqrt(sum_discrete) << " " << 0 << endl; 
			continous_results << sqrt(sum_continous) << " " << 0 << endl;
		} else if(i!=0) {
			discrete_results << sqrt(sum_discrete) << " " << error(sum_discrete, error_discrete, blocks)/(2.*sqrt(sum_discrete)) << endl;  //the error needs to be normalized (it's an error propagation)
			continous_results << sqrt(sum_continous) << " " << error(sum_continous, error_continous, blocks)/(2.*sqrt(sum_continous)) << endl; 
		}
	}

	discrete_results.close();
  	continous_results.close();
  	return 0;
}

double error(double av, double av_squares, int n) { //attenzione che in questo esercizio raccogliamo tutte le distanze medie in una quindi n=blocks ed e' fix, 
													//non gira perche' ogni dato della posizione medie contiene tutti i blocchi 
    if (n == 0) {
        return 0;
    } else {
        return sqrt((av_squares - av * av) / n);
    }
}
