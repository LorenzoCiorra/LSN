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
   	ifstream Primes("Primes");
   	if (Primes.is_open()){
       Primes >> p1 >> p2 ;
   	} else cerr << "PROBLEM: Unable to open Primes" << endl;
   	Primes.close();

   	ifstream input("seed.in");
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

	double M = 10000; //total number of datas for each N 
	double N[4] = {1, 2, 10, 100}; // number of throws of each one of the four experiments

	//wanna create 3 .out, one for each distribution, with 40000 values: the first 10000 are for N=1, the second 10000 for N=2 and so on
	ofstream std_dice;
   	std_dice.open("std_dice.out"); 
	ofstream exp_dice;
   	exp_dice.open("exp_dice.out"); 
	ofstream lor_dice;
	lor_dice.open ("lor_dice.out");
	
	for (int i=0; i<4; i++){
 		for(int j=0; j<M; j++) {
			double std_acc = 0;
			double exp_acc =0;
			double lor_acc =0;
			
			for (int k=0; k<N[i]; k++) {
				double x = rnd.Rannyu();
				std_acc += x;
				double y = rnd.Exponential(1);
				exp_acc +=y;
				double z = rnd.Lorentzian(1);
				lor_acc +=z;
			}
			//now we have the sum for N numbers and we want the mean
			std_acc = std_acc/N[i]; // N[i] dovrebbe gia` essere un double
			exp_acc = exp_acc/N[i]; 
			lor_acc = lor_acc/N[i]; 

			std_dice << std_acc << endl;
			exp_dice << exp_acc << endl;
			lor_dice << lor_acc << endl;
		}
	}

	std_dice.close();
	exp_dice.close();
	lor_dice.close();
	
   	return 0;
}




