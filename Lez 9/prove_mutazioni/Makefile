CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o 
	$(CC) main.o -o main.exe
clear :
	rm main.exe main.o
esegui : main.exe
	./main.exe
make shift_n :
	$(CC) -c shift_n.cpp -o main.o $(CFLAGS)
make inversion :
	$(CC) -c inversion.cpp -o main.o $(CFLAGS)
make m_permutation :
	$(CC) -c m_permutation.cpp -o main.o $(CFLAGS)
