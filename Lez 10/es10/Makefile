CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o
	mpicxx random.o main.o -o main.exe
main.o : main.cpp main.h
	mpicxx -c main.cpp -o main.o $(CFLAGS)
random.o : randomgen/random.cpp randomgen/random.h
	$(CC) -c randomgen/random.cpp -o random.o $(CFLAGS)
clear :
	rm *.o main.exe seed.out
esegui :
	mpiexec -np 10 ./main.exe
clear_results :
	rm -rf *.out
clear_results_closed :
	rm -rf closed/*.out
clear_results_open :
	rm -rf open/*.out
clear_all : clear clear_results

esegui_closed : clear_results_closed esegui
esegui_open : clear_results_open esegui
