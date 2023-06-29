CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : randomgen/random.cpp randomgen/random.h
	$(CC) -c randomgen/random.cpp -o random.o $(CFLAGS)
clear :
	rm *.o main.exe seed.out
esegui :
	./main.exe
clear_results :
	rm *.out
clear_all : clear clear_results
clean:
	rm -f main main-debug