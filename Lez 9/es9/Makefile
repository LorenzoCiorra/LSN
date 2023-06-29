CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp main.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : randomgen/random.cpp randomgen/random.h
	$(CC) -c randomgen/random.cpp -o random.o $(CFLAGS)
clear :
	rm *.o main.exe seed.out
esegui :
	./main.exe
clear_results :
	rm -rf *.out
clear_results_cerchio_500 :
	rm -rf 500/circle_*.out
clear_results_cerchio_2000 :
	rm -rf 2000/circle_*.out
clear_results_quadrato_500 :
	rm -rf 500/square_*.out
clear_results_quadrato_2000 :
	rm -rf 2000/square_*.out
clear_all : clear clear_results

esegui_c_500 : clear_results_cerchio_500 esegui
esegui_s_500 : clear_results_quadrato_500 esegui
esegui_c_2000 : clear_results_cerchio_2000 esegui
esegui_s_2000 : clear_results_quadrato_2000 esegui