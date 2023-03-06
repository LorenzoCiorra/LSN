#include <iostream>
#include <fstream>
#include <string>
#include "randomgen/random.h"
#include <vector>
#include <cmath>

using namespace std;

double error(double av, double av_squares, int n);
 
int main (int argc, char *argv[]){

	//checks for the functioning of the random nmbers generator
   	Random rnd;
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

	//means of <r> and sigma^2 and their errors
   	int M=100000; //number of throws
   	int N=100; //number of blocks
   	int L=int(M/N); //number of throws in each block
  	ofstream meanresults;
   	meanresults.open("meanresults.out"); 
	ofstream sigmaresults;
   	sigmaresults.open("sigmaresults.out"); 
	

	double mean_tot=0; //total mean therefore the sum of the mean in each block
	double mean_error=0; //to calculate the sigma of <r>
	double sigma_tot=0;
	double sigma_error=0;

	for (int i =0; i<N; i++) { //cycle on the blocks
		double mean =0; //mean on each one of the blocks
		double sigma =0;
		for(int j=0; j<L; j++) { //cycle on the throws in each block
			double r = rnd.Rannyu();
			mean +=r;
			sigma += pow(r-0.5,2);
		}
		mean = mean/L;
		mean_tot += mean;
		mean_error += mean*mean;

		sigma = sigma/L;
		sigma_tot += sigma;
		sigma_error += sigma*sigma;
		
		// saving on an external file mean and its error
        meanresults << mean_tot/(i+1) << " " << error(mean_tot/(i+1),mean_error/(i+1),i) <<endl;
		// saving on an external file sigma squared and its error
        sigmaresults << sigma_tot/(i+1) << " " << error(sigma_tot/(i+1),sigma_error/(i+1),i) <<endl;
	}
   	meanresults.close();
	sigmaresults.close();

/*	
	//chi_squared calculation MANCA CAPIRE CHI E` ni
	M = 100;
	double n = 10000;
	double ni[M];
	ofstream chiresults;
	chiresults.open ("chiresults.out");
	
	for(int i=0; i<M; i++) { //cycles on the intervals
		double chi_2 =0;
		for(int j=0; j<n; j++) { //cycle on the throws in each interval
			double r = rnd.Rannyu();
			for(int k=0; k<M; k++) { //cycle to calcultate the ni, proably an array of doubles and with that we calculate the chi^2 in each block
				//??
				chi_2 += pow(ni[i] - n/M,2)/(n/M); 
			}
		}
		//now i have the total chi^2 for this interval and i can write it down
		chiresults << chi_2 << endl;
	}
	chiresults.close();
*/
   	return 0;
}

double error(double av, double av_squares, int n) {
   if(n==0){
      return 0;
   }
   else{
      return sqrt((av_squares-av*av)/n);
   }
}



