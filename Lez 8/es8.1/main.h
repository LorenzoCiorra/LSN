#ifndef __Haverage__
#define __Haverage__

#include "randomgen/random.h"
Random rnd;

const int hbar=1;
const int m=1;
	
int M; //number of throws
int blocks; //number of blocks
int steps; //number of throws in each block
double L;
double accepted =0;
double attempted =0;
double block_ave =0.;
double sum_energy =0.;
double sum2_energy =0.;
double mu;
double sigma; //per ora parto con valori dei parametri presi ad occhio dal plot del GS sul jupyter 
double x, x_try =0.;

void input();
double psi_test (double x, double mu, double s);
double psi_test_2derivative (double x, double mu, double sigma);
double epot (double x);
double energy (double x, double mu, double sigma);
double acceptance(double mu, double sigma, double x_try, double x_old);
void H_average();
double error(double sum, double sum2, int iblk);


#endif