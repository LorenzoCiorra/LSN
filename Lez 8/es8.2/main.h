#ifndef __Haverage__
#define __Haverage__

#include "randomgen/random.h"
Random rnd;

const int hbar=1;
const int m=1;

char answer;	
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

double beta, delta_beta, beta_max;
int beta_steps, progressive_steps=1;
double mu_new, sigma_new;
double L_mu, L_sigma;
double H, H_new, H_error;
double accepted_beta, attempted_beta;

void input();
double psi_test (double x, double mu, double s);
double psi_test_2derivative (double x, double mu, double sigma);
double epot (double x);
double energy (double x, double mu, double sigma);
double acceptance(double mu, double sigma, double x_try, double x_old);
void H_average(double mu, double sigma);
double error(double sum, double sum2, int iblk);


#endif