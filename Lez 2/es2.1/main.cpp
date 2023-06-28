#include <iostream>
#include <fstream>
#include <string>
#include "randomgen/random.h"
#include <vector>
#include <cmath>

using namespace std;

double error(double av, double av_squares, int n);
 
int main (int argc, char *argv[]){

	//checks for the functioning of the random numbers generator
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


   	int M=100000; //number of throws
   	int N=100; //number of blocks
   	int L=int(M/N); //number of throws in each block
  	ofstream uniform_results;
   	uniform_results.open("uniform_results.out"); 
	ofstream impsampl_results;
   	impsampl_results.open("impsampl_results.out"); 
	

	double sum_unif=0;
	double error_unif=0; 
	
	double sum_impsampl=0;
	double error_impsampl=0;

	for (int i =0; i<N; i++) { //cycle on the blocks
		double integral_unif =0; //mean on each one of the blocks
		double integral_impsampl =0;
		for(int j=0; j<L; j++) { //cycle on the throws in each block
			double x = rnd.Rannyu();
			double f_x = (M_PI/2.)*cos((M_PI/2.)*x);
			integral_unif += f_x;
			double x_impsampl = rnd.cosine_density();
			integral_impsampl += ((M_PI/2.)*cos((M_PI/2.)*x_impsampl))/(2.-2.*x_impsampl); //to get the starting function, we need to divide by 2-2x	
		}
		integral_unif = (double)(integral_unif/L); //here I should *(b-a) the lenght of the interval in which we integrate, but in our case it's equals to 1
		sum_unif += integral_unif;
		error_unif += integral_unif*integral_unif;

		integral_impsampl = (double)(integral_impsampl/L);
		sum_impsampl += integral_impsampl;
		error_impsampl += integral_impsampl*integral_impsampl;
		
		// saving on an external file result of the integral and its error sampling a uniform distribution 
        uniform_results << sum_unif/(i+1) << " " << error(sum_unif/(i+1),error_unif/(i+1),i) <<endl;
		// saving on an external file result of the integral and its error sampling a non-uniform distribution 
        impsampl_results << sum_impsampl/(i+1) << " " << error(sum_impsampl/(i+1),error_impsampl/(i+1),i) <<endl;
	}
   	uniform_results.close();
	impsampl_results.close();
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
