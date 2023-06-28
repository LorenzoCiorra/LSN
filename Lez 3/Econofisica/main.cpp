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
   	int blocks=100; //number of blocks
   	int steps=int(M/blocks); //number of throws in each block

	double S_0 = 100.;
	double delivery_time =1.;
	double strike_price =100.; //strike price
	double interest_rate =0.1; //this is r and also mu the mean of tha gaussian
	double volatility = 0.25;

	double S_T=0;

	//sampling directly the final asset price
  	ofstream directly_call_results;
   	directly_call_results.open("directly_call_results.out"); 
	ofstream directly_put_results;
   	directly_put_results.open("directly_put_results.out"); 
	
	double sum_directly_call=0;
	double error_directly_call=0; 
	double sum_directly_put=0;
	double error_directly_put=0;

	//sampling the discretized gbm path of the asset price
  	ofstream discrete_call_results;
   	discrete_call_results.open("discrete_call_results.out"); 
	ofstream discrete_put_results;
   	discrete_put_results.open("discrete_put_results.out"); 
	
	double sum_discrete_call=0;
	double error_discrete_call=0; 
	double sum_discrete_put=0;
	double error_discrete_put=0;
	double total_time =100;

	
	for (int i =0; i<blocks; i++) { //cycle on the blocks
		
		double directly_call=0;
		double directly_put=0;

		double discrete_call=0;
		double discrete_put=0;
		
		for(int j=0; j<steps; j++) { //cycle on the throws in each block
			
			double z = rnd.Gauss(0,1);
			S_T = S_0 * exp( (interest_rate - pow(volatility,2)/2.)*delivery_time + volatility*z*sqrt(delivery_time));
			directly_call += exp(-delivery_time*interest_rate)*max(0., S_T - strike_price);
			directly_put += exp(-delivery_time*interest_rate)*max(0., strike_price - S_T);

			//cycle for discretized path of the asset price
			for(int k=0; k<total_time; k++ ){ //when k=total_time-1, it's calculating s(T)=s(100) using s(99) 
				double z_2 = rnd.Gauss(0,1);
				double time_difference = delivery_time/100.; //t_(i+1)-t_(i) is always equal to 1/100
				if(k==0) {
					S_T = S_0 * exp( (interest_rate - pow(volatility,2)/2.)*time_difference + volatility*z_2*sqrt(time_difference));
				} else if(k!=0) {
					S_T = S_T * exp( (interest_rate - pow(volatility,2)/2.)*time_difference + volatility*z_2*sqrt(time_difference));
				}
			}
			discrete_call += exp(-delivery_time*interest_rate)*max(0., S_T - strike_price);
			discrete_put += exp(-delivery_time*interest_rate)*max(0., strike_price - S_T);
		}
		//sampling directly the final asset price
		directly_call = (double)(directly_call/steps);
		sum_directly_call += directly_call;
		error_directly_call += directly_call*directly_call;

		directly_put = (double)(directly_put/steps);
		sum_directly_put += directly_put;
		error_directly_put += directly_put*directly_put;
		
		// saving on an external file result of the directly sampled call adn put option price  
        directly_call_results << sum_directly_call/(i+1) << " " << error(sum_directly_call/(i+1),error_directly_call/(i+1),i) <<endl;
        directly_put_results << sum_directly_put/(i+1) << " " << error(sum_directly_put/(i+1),error_directly_put/(i+1),i) <<endl;

		
		//sampling the discretized gbm path of the asset price
		discrete_call = (double)(discrete_call/steps);
		sum_discrete_call += discrete_call;
		error_discrete_call += discrete_call*discrete_call;

		discrete_put = (double)(discrete_put/steps);
		sum_discrete_put += discrete_put;
		error_discrete_put += discrete_put*discrete_put;
		
		// saving on an external file result of the directly sampled call adn put option price  
        discrete_call_results << sum_discrete_call/(i+1) << " " << error(sum_discrete_call/(i+1),error_discrete_call/(i+1),i) <<endl;
        discrete_put_results << sum_discrete_put/(i+1) << " " << error(sum_discrete_put/(i+1),error_discrete_put/(i+1),i) <<endl;
	}
	
   	directly_call_results.close();
	directly_put_results.close();
	discrete_call_results.close();
	discrete_put_results.close();
	
	
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
