#include "main.h"
#include "mpi.h"

using namespace std;

void input();
void declaration(int size, int rank);
bool check(vector <unsigned int> v);
unsigned int rigged_selection ();
void print_c (vector<unsigned int> chromosome);
void genetic_algorithm_core();
void print_distances(std::string rank_string, fitness fitness_function);
void print_cartesian_path (std::string rank_string, vector<city> county, vector<unsigned int> the_best_path);


int main(int argc, char* argv[]) {
    input();
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat1, stat2;
    MPI_Request request;
    cout << "size: " << size << endl;
    rnd.set_prime(rank); //ogni nodo viene settato su una coppia di numeri primi diversi
    rank_string = to_string(rank);
    cout <<"rank" <<rank_string << endl;
    declaration(size, rank);
    vector <int> rnd_list_cores(size);


    vector <city> county;
    int borders_int=0;
    if (rank==0) { //il primo core mi chiede se permettere le migrazioni o no, poi lui lo fa sapere agli altri
        cout << "Do you want to allow migrations or not? [y/n]" << endl;
        char answer;
        cin >> answer;
        if (answer == 'y' ) { 
            borders_int = 1;
        }else if (answer == 'n') { 
            borders_int = 2;
        } else {
            cout << "Input, not valid. Try again!" << endl;
            return 1;
        }
    }
    MPI_Bcast(&borders_int, 1, MPI_INT, 0, MPI_COMM_WORLD); //0 fa sapere a tutti come sono i confini //bcast di string è un incubo
    if(borders_int==1) {
        borders = "open";
    } else if(borders_int==2) {
        borders = "closed";
    }

    american_capitals.open("American_capitals.dat");
    list_of_capitals.open("list_of_capitals.dat");
    vector <string> list_capitals;
    std::string temporary;
    if(american_capitals.is_open() && list_of_capitals.is_open()){
        for(int i = 0; i<nr_cities; i++){
            american_capitals >> x >> y;
            county.push_back(city(x, y));
            list_of_capitals >> temporary;
            list_capitals.push_back(temporary); //cosi' poi potrei stampare il migliore percorso come lista di citta'
        }
    }else{
        cout << "Unable to read American_capitals.dat" << endl; 
        exit(1);
    }
    chromosome_number_string = to_string(nr_chromosomes);

    fitness fitness_function(county);

	for(int i=0; i<(int)chromosome.size(); i++) {
		chromosome[i]=i;
	} 
    check(chromosome);

	//per creare la popolazione parto da chromosome generato, lo mischio e creo population[1], lo rimischio e ho population[2] ...
	random_device rd;
	mt19937 g(rd());
	for(int i=0; i<nr_chromosomes; i++) {
		shuffle(chromosome.begin()+1, chromosome.end(), g); //shuffle dal secondo elemento cosi' il primo e' sempre 0 in tutti
		population[i] = chromosome;
        check(population[i]);
	}
    for(int b=0; b<size; b++) { //creo una lista dei cores da usare per le migrazioni
        rnd_list_cores[b]=b;
    }

	//ordino i cromosomi in popolazione in base al loro fitness: l'ultimo e' il cromosoma migliore
	sort(population.begin(), population.end(), [&fitness_function](vector <unsigned int> &a, vector <unsigned int> &b){return fitness_function.path_fitness(a)<fitness_function.path_fitness(b);}); 

	for(auto i=0; i<nr_max_generations; i++ ) {
		next_population.clear(); //ho bisogno di cancellare tutto cio' che c'e' dentro e non solo settarlo a cromosomi tutti di zero
        progeny.clear();
		average_distance=0.;
        if(borders == "closed" || i == 0 || migration_index%migration_generation!=0 ) { //continenti non si scambiano informazioni
            genetic_algorithm_core();
        } else if(migration_index%migration_generation==0) { //i continenti si parlano tra loro
            if(rank==0) {
                cout << "Migration process in action!" << endl;
            }
            //int addressee_core=0, sender_core=1; //devo usare MPI_UNSIGNED
            int flag_1=1, flag_2=2;

            //ora mischio la lista dei core e le migrazioni avvengono tra il primo e il secondo della lista, tra il terzo e il quarto e cosi' via
            //in questo modo a coppie tutti i continenti si scambiano un cromosoma, ma con chi lo scambiano e' casuale e puo' cambiare ad ogni migrazione
            if(rank==0) {
                shuffle(rnd_list_cores.begin(), rnd_list_cores.end(), g);  
            }
            //adesso il core 0 fa sapere a tutti gli altri la lista dei core  da seguire per fare la migrazione
            MPI_Bcast(&rnd_list_cores[0], size, MPI_INTEGER, 0, MPI_COMM_WORLD);
            for(int a=0; a<size-1; a=a+2) {//se i core sono in numero totale pari tutti subiscono uno scambio, se in numero dispari l'ultimo non scambia. 
            //l'aggiornamento e' a+2 perche' cosi' a diventa il primo elemento della nuova coppia
                if(rank == a) { //lui e' il sender
                    ambassador_sender = population.back();
                    MPI_Isend(&ambassador_sender[0], ambassador_sender.size(), MPI_UNSIGNED, a+1, flag_1, MPI_COMM_WORLD, &request);
                    MPI_Recv(&ambassador_addressee[0], ambassador_addressee.size(), MPI_UNSIGNED, a+1, flag_2, MPI_COMM_WORLD, &stat2);
                    population.back() = ambassador_addressee;
                    check(population.back());
                } else if(rank == a+1) { //lui il destinatario
                    ambassador_addressee = population.back();
                    MPI_Send(&ambassador_addressee[0], ambassador_addressee.size(), MPI_UNSIGNED, a, flag_2, MPI_COMM_WORLD);
                    MPI_Recv(&ambassador_sender[0], ambassador_sender.size(), MPI_UNSIGNED, a, flag_1, MPI_COMM_WORLD, &stat1);
                    population.back() = ambassador_sender;
                    check(population.back());                
                }
 
                //sorto cosi' da rimettere in ordine dopo gli scambi
                sort(population.begin(), population.end(), [&fitness_function](vector <unsigned int> &a, vector <unsigned int> &b){return fitness_function.path_fitness(a)<fitness_function.path_fitness(b);}); 
                genetic_algorithm_core();
            }
        }

        population = next_population; //prossima generazione: i piccoli diventano adulti e sono pronti ad avere figli!
        //e il ciclo della vita ricomincia con la nuova generazione...

        //prima di andare avanti con la prossima generazione abbiamo bisogno di salvare i risultati ottenuti
        //prima sortiamo la popolazione in base al loro fitness_score
        sort(population.begin(), population.end(), [&fitness_function](vector <unsigned int> &a, vector <unsigned int> &b){return fitness_function.path_fitness(a)<fitness_function.path_fitness(b);}); 
        the_best_path = population.back();
        print_distances(rank_string, fitness_function);

        generation_index ++;
        migration_index ++;
        if(rank==0) {
            cout << "\nGeneration: " <<generation_index << endl << endl;
            cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
        }
	}
    //arrivati al numero massimo di generazioni, abbiamo in the_best_path il cammino migliore possibile e vogliamo salvarlo come coordinate cartesiane delle citta' corrispondenti
    print_cartesian_path(rank_string, county, the_best_path);

    MPI_Finalize();

    cout << "Simulation finished!" << endl;
	return 0;
}


//~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~~~//

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

bool check(vector<unsigned int>v) {
    if(v.front()!=0) {
    cout << "Something went wrong, the first gene is not 0!" << endl;
        exit(1);
    }
    int sum =0.;
    for(int i=0; i<(int)v.size(); i++) {
        sum+=v[i];
    }
    if(sum!=1225) {
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

void declaration(int size, int rank) {
    if(rank==0) {
        cout << endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
        cout << "Parellalized Genetic Algorithm Simulation of the Traveling Salesman problem        " << endl<<endl;
        cout << "Performed with       " << size << " cores and " <<endl;
        cout << nr_chromosomes <<" chromosomes" << endl;
        cout << "each one with " << nr_genes << " genes" << endl;
        cout << "representing the " << nr_cities << " American Capitals which the Salesman must visit" << endl << endl;
        cout << "The simulation will resolve the problem using " << nr_max_generations << " generations." << endl << endl;
        cout << "It will use the following probabilities: " << endl;
        cout<< "_Crossover probability: " << crossover_prob << endl;
        cout<< "_Selection probability: " << selection_prob << endl;  
        cout<< "_Mutation probability: " << mutation_prob << endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" <<endl;
    }
}

void genetic_algorithm_core() {
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
        check(father);

        if(rnd.Rannyu() < crossover_prob) {
            genetic_functions incubator (father);
            progeny = incubator.crossover(mother);
            child_1 = progeny[0];
            child_2 = progeny[1];
            check(child_1);
            check(child_2);
        } else {
            child_1 = father;
            child_2 = mother;
        };

		//definiamo l'oggetto evoluzione e mutiamo i figli
		genetic_functions evolution_1 (child_1);
		evolution_1.pair_permutation(mutation_prob);
		evolution_1.shift_of_n(mutation_prob);      
		evolution_1.m_permutation(mutation_prob);
		evolution_1.inversion(mutation_prob);
        child_1 = evolution_1.get_chromosome();
		check(child_1);
			
        //facciamo lo stesso per il figlio minore
        genetic_functions evolution_2(child_2);
		evolution_2.pair_permutation(mutation_prob);
		evolution_2.shift_of_n(mutation_prob);      
		evolution_2.m_permutation(mutation_prob);
		evolution_2.inversion(mutation_prob);            
        child_2 = evolution_2.get_chromosome();;
		check(child_2);
            
        next_population.push_back(child_1);            
		next_population.push_back(child_2);
            

        mother.clear();
        father.clear();
        child_1.clear();
        child_2.clear();
        progeny.clear();
    }
}

void print_distances(std::string rank_string, fitness fitness_function) {

  	out_distance.open(borders + "/" + rank_string +  "_distance_results.out", ios::app);
  	out_average_distance.open(borders + "/" + rank_string +  "_average_distance_results.out", ios::app);

    out_distance << generation_index << " " << 1./fitness_function.path_fitness(the_best_path) << endl;
    for(auto k= population.size()/2; k<population.size(); k++) {
        average_distance += 1./fitness_function.path_fitness(population[k]);
    }
    out_average_distance << generation_index << " " << average_distance/((double)nr_chromosomes/2.) << endl;

    out_distance.close();
    out_average_distance.close();
}

void print_cartesian_path (std::string rank_string, vector<city> county, vector<unsigned int> the_best_path) {
    out_cartesian_path.open(borders + "/" + rank_string + "_cartesian_path_results.out", ios::app);
    for(int i=0; i< (int) the_best_path.size(); i++) {
        out_cartesian_path << county[the_best_path[i]].get_x() << " " << county[the_best_path[i]].get_y() <<endl;
    } //stampiamo anche la prima (ultima) citta' cosi' il plot in jupyter unisce bene tutte le citta'
    out_cartesian_path << county[the_best_path.front()].get_x() << " " << county[the_best_path.front()].get_y() <<endl;
    out_cartesian_path.close();
}

