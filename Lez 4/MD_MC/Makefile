CC = g++
CFLAGS = -Wall -O3 --std=c++11

NVE_NVT.exe : MD_MC.o random.o
	$(CC) random.o MD_MC.o -o NVE_NVT.exe
MD_MC.o : MD_MC.cpp
	$(CC) -c MD_MC.cpp -o MD_MC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clear :
	rm *.o NVE_NVT.exe
esegui :
	./NVE_NVT.exe
clear_results_solido :
	rm -rf solido/*.dat
clear_results_liquido :
	rm -rf liquido/*.dat
clear_results_gas :
	rm -rf gas/*.dat
clear_frames :
	rm -rf frames/*.xyz
esegui_solido : clear_results_solido esegui
esegui_liquido : clear_results_liquido esegui
esegui_gas : clear_results_gas esegui