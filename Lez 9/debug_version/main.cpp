#include "main.h"

using namespace std;

void input();
bool check(vector <unsigned int> v);
unsigned int rigged_selection ();
void print_c (vector<unsigned int> chromosome);



int main() {
    input();
    vector <city> county;
    cout << "Do you want to put the city on a circumference or inside a square? [c/s]" << endl;
    char answer;
    cin >> answer;
    if (answer == 'c' ) { 
        geometry = "circle";
        cout << geometry;
        for(auto i =0; i<nr_cities; i++) {
            theta = rnd.Rannyu (0, 2*M_PI);
            county.push_back(city(cos(theta), sin(theta))); 
        }
    }else if (answer == 's') { 
        geometry = "square";
        for(auto i =0; i<nr_cities; i++) {
            x = rnd.Rannyu(-1, 1);
            y = rnd.Rannyu(-1, 1);
            county.push_back(city(x, y)); 
        }
    } else {
        cout << "Input, not valid. Try again!" << endl;
        return 1;
    }

    fitness fitness_function(county);
    vector <unsigned int> chromosome(nr_genes);
    vector< vector<unsigned int> > population (nr_chromosomes);//la popolazione e' un vettore di 500 cromosomi

	vector<unsigned int> the_best_path; //qui salvero' in ogni generazione il miglior percorso
  	vector<unsigned int> lad_1, lad_2, mother, father, child_1, child_2; //definiamo i cromosomi su cui lavoreremo con le funzioni definite in genetica
  	vector<vector<unsigned int>> next_population; //evoluzione della popolazione 
    vector <vector<unsigned int>> progeny(2);


	for(int i=0; i<(int)chromosome.size(); i++) {
		chromosome[i]=i;
	} 
    check(chromosome);
    print_c(chromosome);

	//per creare la popolazione parto da chromosome generato, lo miscio e creo population[1], lo rimischio e ho population[2] ...
	random_device rd;
	mt19937 g(rd());
	for(int i=0; i<nr_chromosomes; i++) {
		shuffle(chromosome.begin()+1, chromosome.end(), g); //shuffle dal secondo elemento cosi' il primo e' sempre 0 in tutti
		population[i] = chromosome;
        check(population[i]);
        print_c(chromosome);
	}

        cout << "problemi col sort" << endl;
	//ordino i cromosomi in popolazione in base al loro fitness: l'ultimo e' il cromosoma migliore
	sort(population.begin(), population.end(), [&fitness_function](vector <unsigned int> &a, vector <unsigned int> &b){return fitness_function.path_fitness(a)<fitness_function.path_fitness(b);}); 
    //cosi' il sort funziona, se no creare fitness_function come var globale a costruttore nullo e fare un fitness.set_cities nel main con county
    
    for(int i=0; i<nr_chromosomes; i++) {
        cout << fitness_function.path_fitness(population[i]) << " ";
        print_c(population[i]);
    }
    cout << endl;
    cout << "qui gia va bene" << endl;


  	out_cartesian_path.open(geometry + "_cartesian_path_results.out", ios::app);
  	out_distance.open(geometry + "_distance_results.out", ios::app);
  	out_average_distance.open(geometry + "_average_distance_results.out", ios::app);
//      out_average_elegante.open(geometry + "_average_elegante_results.out", ios::app);

	//for(auto i=0; i<nr_max_generations; i++ ) {
    for(auto i=0; i<100; ++i ) {
		next_population.clear(); //ho bisogno di cancellare tutto cio' che c'e' dentro e non solo settarlo a cromosomi tutti di zero
        progeny.clear();
		average_distance=0.;
//      average_elegante=0.;
		while(next_population.size()<nr_chromosomes) {

			//scegliamo due individui che poi diventeranno i genitori
            int selection_1 = rigged_selection();
            int selection_2 = rigged_selection();
			mother = population[selection_1];
			father = population[selection_2];
            while(selection_1==selection_2) {
                selection_2 = rigged_selection();
                father = population[selection_2];
            }
			check(mother);
//            cout << "father"<< endl;
//            print_c(father);
//            cout << "pre crossover selection" << endl;
            check(father);
//            cout << "father"<< endl;
//            print_c(father);
//            cout << "pre crossover selection" << endl;

            if(rnd.Rannyu() < crossover_prob) {
                genetic_functions incubator (father);
                progeny = incubator.crossover(mother);
//                cout << "progeny[0] pre mutazione" << endl;
//                print_c(progeny[0]);
                child_1 = progeny[0];
                child_2 = progeny[1];
                check(child_1);
                check(child_2);
            } else {
                child_1 = father;
                child_2 = mother;
            }
            cout << "figlio 1 pre mutazione" << endl;
            print_c(child_1);
//            cout << "figlio 2 pre mutazione" << endl;
//            print_c(child_2);


			//definiamo l'oggetto evoluzione e mutiamo i figli
			genetic_functions evolution_1 (child_1);
			evolution_1.pair_permutation(mutation_prob);
			evolution_1.shift_of_n(mutation_prob);      
			evolution_1.m_permutation(mutation_prob);
			evolution_1.inversion(mutation_prob);
            child_1 = evolution_1.get_chromosome();
                cout << "primo figlio post mutazione" << endl;
                print_c(child_1);
			check(child_1);
//            cout << "primo figlio" << endl;
			
            //facciamo lo stesso per il figlio minore
            genetic_functions evolution_2(child_2);
			evolution_2.pair_permutation(mutation_prob);
			evolution_2.shift_of_n(mutation_prob);      
			evolution_2.m_permutation(mutation_prob);
			evolution_2.inversion(mutation_prob);            
            child_2 = evolution_2.get_chromosome();
//              cout << "secondo figlio" << endl;
//              print_c(child_2);
			check(child_2);
            
//              cout << "secondo figlio" << endl;
            next_population.push_back(child_1);            
			next_population.push_back(child_2);
            

            mother.clear();
            father.clear();
            child_1.clear();
            child_2.clear();
            progeny.clear();

		}
        population = next_population; //prossima generazione: i piccoli diventano adulti e sono pronti ad avere figli!
        //e il ciclo della vita ricomincia con la nuova generazione...

//          cout << "population" << endl;
//          for(int b=0; b<6; b++) {
//              print_c(next_population[b]);
//          }

        //prima di andare avanti con la prossima generazione abbiamo bisogno di salvare i risultati ottenuti
        //prima sortiamo la popolazione in base al loro fitness_score
        sort(population.begin(), population.end(), [&fitness_function](vector <unsigned int> &a, vector <unsigned int> &b){return fitness_function.path_fitness(a)<fitness_function.path_fitness(b);}); 

        if(generation_index == 0) {
            the_best_path = population.back();
        }
        if(fitness_function.path_fitness(the_best_path)<fitness_function.path_fitness(population.back())) {
            the_best_path = population.back();
        }
        out_distance << generation_index << " " << 1./fitness_function.path_fitness(the_best_path) << endl;
        for(auto k= population.size()/2; k<population.size(); k++) {
            average_distance += 1./fitness_function.path_fitness(population[k]);
        }
        out_average_distance << generation_index << " " << average_distance/((double)nr_chromosomes/2.) << endl;
/*        
        for(auto k =population.begin() + nr_chromosomes/2.; k!=population.end(); k++ ) { //stessa cosa di sopra,ma piu' elegante
            average_elegante += 1./fitness_function.path_fitness(*k); // *k e' il valore della cella puntata da k quindi e' un population[k]
        }
        out_average_elegante << generation_index << " " << average_elegante/((double)nr_chromosomes/2.) << endl;
*/


        generation_index ++;
        cout << "GENERAZIONE" <<generation_index << endl;
	}
    //arrivati al numero massimo di generazioni, abbiamo in the_best_path il cammino migliore possibile e vogliamo salvarlo come coordinate cartesiane delle citta' corrispondenti
    for(int i=0; i< (int) the_best_path.size(); i++) {
        out_cartesian_path << county[the_best_path[i]].get_x() << " " << county[the_best_path[i]].get_y() <<endl;
    } //stampiamo anche la prima (ultima) citta' cosi' il plot in jupyter unisce bene tutte le citta'
    out_cartesian_path << county[the_best_path.front()].get_x() << " " << county[the_best_path.front()].get_y() <<endl; //front rida' il primo elemento, begin un puntatore al primo

	out_cartesian_path.close();
	out_distance.close();
	out_average_distance.close();

	return 0;
}


unsigned int rigged_selection () {
    unsigned int j = (unsigned int)(nr_chromosomes*pow(rnd.Rannyu(),selection_prob));
    if(j>=nr_chromosomes) {
        j=nr_chromosomes-1;
    }
    if(j<0) {
        j=0;
    }
    return j;
}
/*
bool check(vector <int> v) { //voglio vedere se quello che passo e' una permutazione di un vettore che ha tutti i numeri tra 0 e 33
    if(v[0]!=0) {
        cout << "Something went wrong, the first gene is not 0!" << endl;
        exit(1);
    }
	vector <int> comparison;
	for(int i=0; i<(int)v.size(); i++) {
		comparison[i]=i;
	}
	return is_permutation(v.begin(), v.end(), comparison.begin());  //se e' una permutazione va bene perche' contiene tutti elementi diversi. True=accettabile
}
*/
bool check(vector<unsigned int>v) {
    if(v.front()!=0) {
    cout << "Something went wrong, the first gene is not 0!" << endl;
        exit(1);
    }
    int sum =0.;
    for(int i=0; i<(int)v.size(); i++) {
        sum+=v[i];
    }
    if(sum!=561) {
        cout << "Problemone!" << endl;
        print_c(v);
        exit(1);
    }
    return 1;
}

void input() {
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
}

void print_c (vector<unsigned int> chromosome) {
    for(int a=0; a<(int)chromosome.size(); a++) {
        cout << chromosome[a] << " " ;
    }
    cout << endl;
}

