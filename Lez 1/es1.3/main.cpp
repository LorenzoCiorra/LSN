#include "randomgen/random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

double error(double av, double av_squares, int n);

int main(int argc, char *argv[]) {

    // checks for the functioning of the random nmbers generator
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("randomgen/Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else
        cerr << "PROBLEM: Unable to open Primes" << endl;
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

    double d = 5.0; // distance between two vertical lines: in this experiment we suppose to use only 2 lines
    double l = 3.;
    double num = 0.;
    double pi = 0.;
    double x2 =1; //1 e non 0 perche cosi non triggera l'if del n_hit anche se non calcolo x2 perche la norma e' >1

    double pi_tot = 0.;
    double pi_sqr = 0.;
    int M = 1000000;    // total number of throws
    int N = 100;        // number of blocks
    int T = int(M / N); // number of throws in each block

    ofstream pi_results;
    pi_results.open("pi_results.out");

    for (int i = 0; i < N; i++) {
        int n_hit = 0;
        num=0;
        for (int j = 0; j < T; j++) {
            // now we can define what is the needle and which are its coordinates
            double x1 = rnd.Rannyu(0., d); // the starting point of the needle
            double a = rnd.Rannyu(-1., 1.); //to generate the end of the needle: it must be on a circle of radius l with centre x1
            double b = rnd.Rannyu(-1., 1.);
            double norm = a*a + b*b;
            if(norm <= 1){
                x2 = x1+ l* a/sqrt(norm); //x coordinate of the end of the needle; we don't need its y
                // now that we have the needle on our table, we should check if it touches, or not, one of the two vertical planes
                if ( x2 >= d || x2 <= 0 ) {
                    n_hit++;
                }
            } else{ //the end of the needle is not on the circle: impossibile in reality! We try throwing the needle again wothout considering this unrealistic event
                j--;
            }
        }
        // now we can calcultate pi following Laplace anylis and utilizing the block method to make a more precise estimation
        num = (2.0 * l * (double)T) / (double)(n_hit);
        pi =  num/d; // quello che per la <r> era la media, ora e` il calcolo di pi. Questa riga cambia in base a che calcolo sto facendo; il seguente rimane uguale
        pi_tot += pi;
        pi_sqr += pi * pi;
        pi_results << pi_tot / (i + 1) << " "<< error(pi_tot / (i + 1), pi_sqr / (i + 1), i) << endl;
    }

    pi_results.close();

    return 0;
}

double error(double av, double av_squares, int n) {
    if (n == 0) {
        return 0;
    } else {
        return sqrt((av_squares - av * av) / n);
    }
}
