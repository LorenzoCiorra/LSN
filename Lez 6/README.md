# Istruzioni di Funzionamento

Per andare ad eseguire i codici e' sufficiente compilare con  
```
make
```
ed eseguire il codice con 
```
make esegui 
```

\
Inoltre, il comando 
```
make clear
```
permette di cancellare tutti i file oggetto e l'eseguibile main.exe.

\
Dato che tutti i file di output verranno aperti in append, e' necessario pulirli prima di far partire una nuova simulazione. \
Con
```
make clear_temp
```
possiamo pulire la cartella output contentente i valori delle grandezze in esame per ogni temperatura e per ogni blocco. \
Con 
```
make clear_&algoritmo&
```
sostituendo ad &algoritmo& gibbs o metro, possiamo pulire la cartella risultati_finali/&algoritmo& dai valori finali delle grandezze in funzione della temperatura.\
Il comando 
```
make clear_all
```
raggruppa tutti i comandi di clear precedenti sia per i risultati ottenuti con l'algoritmo di Gibbs che con quello di Metropolis.
