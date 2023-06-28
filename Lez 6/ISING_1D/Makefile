CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_ISING_1D.exe : Monte_Carlo_ISING_1D.o random.o
	$(CC) random.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D.exe
Monte_Carlo_ISING_1D.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clear :
	rm *.o Monte_Carlo_ISING_1D.exe
esegui :
	./Monte_Carlo_ISING_1D.exe
clear_results :
	rm -rf *.out
clear_gibbs :
	rm -rf Risultati_finali/Gibbs/*.out
clear_metro :
	rm -rf Risultati_finali/Metropolis/*.out
clear_temp :
	rm -rf output/*/*.out
clear_all : clear_temp clear_gibbs clear_metro clear