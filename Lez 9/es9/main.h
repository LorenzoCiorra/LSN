#ifndef __GeneticAlgorithm__
#define __GeneticAlgorithm__

#include "randomgen/random.h"
#include <iostream>
#include <string>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;


Random rnd;

std::string geometry, chromosome_number_string; 
double x, y, theta;
int nr_cities=34;
int nr_genes=34;
int nr_chromosomes=2000; //500 e 2000
double crossover_prob = 0.7, mutation_prob = 0.05, selection_prob = 0.1; //p basso cosi' da selezionare i cromosomi migliori
unsigned int nr_max_generations = 300; 
int generation_index = 0;
double average_distance=0.;
//double average_elegante;
ofstream out_cartesian_path, out_distance, out_average_distance;
//ofstream out_average_elegante;


class city { //calsse citta' per dare una posizione nello spazio ad ognuna di esse
	private:
	double m_y, m_x;

	public:
	city () {m_x=0; m_y=0; };
	city(double x, double y) {m_x=x; m_y=y;};
	~city() {;};

	double distance (city other) {
		return sqrt(pow(other.m_x-m_x,2) + pow(other.m_y-m_y,2)); //norma il L(2)
	}

	double get_x () {return m_x;}; 
	double get_y () {return m_y;};
};


class fitness { //prende un vettore di citta' cosi' da calcolare per un cromosoma, che passo come argomento della funzione path_fitness, la sua funzione fitness
	private:
	vector <city> m_city;

	public:
  	fitness(vector <city> city) {m_city = city;}; //i vector hanno = per copiarne una copia
  	~fitness(){;};

  	double path_fitness(vector <unsigned int> chromosome) {
    	double distance = 0;
		for(auto i =0; i<chromosome.size(); i++) {
			//if(i!=chromosome.size()-2) { //perche' i+1 con l'ultima citta' non ha senso
            if(chromosome[i]!= chromosome.back()) { //analogo a sopra, ma piu' elegante
				distance += m_city[chromosome[i+1]].distance(m_city[chromosome[i]]);
			}
		}
		//distanza tra l'ultima e la prima
		distance += m_city[chromosome.back()].distance(m_city[chromosome[0]]);
	
		return 1./distance; //vogliamo minimizzare la funzione costo quindi massimizzare il fitness 
	    //percio' il cromosoma con distanza tot piu' piccola e' quello con migliore (piu' alto) fitness
	}
};


class genetic_functions { //la creo su un cromosoma e ci posso applicare le mutazioni. Per il crossover gli do' come argomento il secondo cromosoma
	private:
	vector <unsigned int> m_chromosome;
	Random m_rnd;

	public:
	genetic_functions (vector <unsigned int> chromosome) {m_chromosome=chromosome; m_rnd = rnd;}; 
	~genetic_functions () {;};

	vector <unsigned int> get_chromosome() { return m_chromosome;};

	void pair_permutation (double probability) { //scambio due elementi tra loro. Tutti gli altri invariati
		if(m_rnd.Rannyu() < probability) {
			unsigned int l = (unsigned int)m_rnd.Rannyu(1, 34); //parte da 1 visto che la prima citta' deve rimanere sempre invariata
			unsigned int k = (unsigned int)m_rnd.Rannyu(1, 34);
			while(l==k) {
				k = (unsigned int)m_rnd.Rannyu(1, 34);
			}
			swap(m_chromosome[l], m_chromosome[k]);
		}
	}

	void shift_of_n (double probability) {
		if(m_rnd.Rannyu() < probability) {
			unsigned int n = (unsigned int)m_rnd.Rannyu(1, m_chromosome.size());
			rotate(m_chromosome.begin()+1,m_chromosome.begin()+1+n, m_chromosome.end()); //parto dal secondo gene e lo ruoto di n posizioni. Ripeto fino ad aver ruotato n geni
			//il +1 serve per non muovere mai il primo gene (funziona! prove_mutazioni/shift_n)
		}
	}

	void m_permutation (double probability) { //permuto a blocchi, praticamente e' uno swap tra blocchi con una lunghezza diversa da 1 (i.e. il singolo elemento)
		if(m_rnd.Rannyu() < probability) {
			unsigned int block_start = (unsigned int) m_rnd.Rannyu(1, m_chromosome.size()-1); //parto da 1 cosi' escludo direttamente il prendere il primo gene
			unsigned int block_finish = (unsigned int) m_rnd.Rannyu(block_start +1, m_chromosome.size());
			unsigned int n = (unsigned int) m_rnd.Rannyu(1, block_finish - block_start); //se n fosse block_finish ritroverei il vettore di partenza
			rotate (m_chromosome.begin()+block_start, m_chromosome.begin()+n+block_start, m_chromosome.begin()+block_finish); //(funziona! prove_mutazioni/m_permutation)
			//ruoto a sx il blocco che inizio in start e finisce in finish di n posizioni rispetto pero' a dove comincia il blocco che stiamo studiando, non rispetto al gene[1]
		}
	}

	void inversion (double probability) { //inverto i geni di un blocco di lunghezza estratta
		unsigned int block_start = (unsigned int) m_rnd.Rannyu(1, m_chromosome.size()-1);
		unsigned int block_finish = (unsigned int) m_rnd.Rannyu(block_start +1, m_chromosome.size());
		reverse (m_chromosome.begin()+block_start, m_chromosome.begin()+block_finish); //(funziona! prove_mutazioni/inversion)
	}


	vector <vector<unsigned int>> crossover(vector <unsigned int> chromosome_2) { //la probabilita' la metto nel main.cpp
		unsigned int path_cut = (unsigned int) m_rnd.Rannyu(1, 34); //tagliare a 0 non ha senso perche' vorrebbe dire prima del primo elemento

        vector <unsigned int> goodpart_1st, goodpart_2nd;
		vector <unsigned int> cutpart_1st(m_chromosome.size()-path_cut), cutpart_2nd(m_chromosome.size()-path_cut);
		vector <unsigned int> crossparts_1st, crossparts_2nd;
		for(auto i =0; i<path_cut; i++) { //salvo la prima parte dei genitori
			goodpart_1st.push_back(m_chromosome[i]);
			goodpart_2nd.push_back(chromosome_2[i]);
		}

		for(auto i = path_cut; i<m_chromosome.size(); i++) { //salvo la parte tagliata dei genitori
			cutpart_1st[i-path_cut]=m_chromosome[i]; 
			cutpart_2nd[i-path_cut]=chromosome_2[i];
		}

		for(auto i =0; i<m_chromosome.size(); i++) { //questo for deve andare su tutto il vettore iniziale nel caso in cui lo stesso numero si trovi nelle parti tagliate di entrambi i vettori
			for(auto j =0; j<m_chromosome.size()-path_cut; j++) { //m_chromosome.size()-path_cut e' la lunghezza delle parti tagliate
				if(m_chromosome[i] == cutpart_2nd[j]) {
					crossparts_2nd.push_back(cutpart_2nd[j]);
				}				
                if(chromosome_2[i] == cutpart_1st[j]) {
					crossparts_1st.push_back(cutpart_1st[j]);
				}

			}
		}

		goodpart_1st.insert(goodpart_1st.end(), crossparts_1st.begin(), crossparts_1st.end());
		goodpart_2nd.insert(goodpart_2nd.end(), crossparts_2nd.begin(), crossparts_2nd.end());

        vector < vector<unsigned int>> progeny = {goodpart_1st, goodpart_2nd};
        return progeny;
	}
};
#endif