# Istruzioni di Funzionamento 

Per andare ad eseguire i codici e' necessario compilare con  
```
make
```
ed eseguire il codice esplicitando la fase di interesse con 
```
make esegui_&borders_situation&
```
rimpiazzando &borders_situation& con open o closed a seconda che si voglia fare una simulazione in cui le migrazioni tra continenti sono permesse oppure no.
In questo modo, il makefile prima di eseguire il codice automaticamente pulira' i risultati della simulazione precedente riguardanti quella situazione. Pulire gli output prima di far ripartire una simulazione e' necessario dato che i file verranno aperti in append e, quindi, se non si pulissero i file di output conterrebbero misure di simulazioni differenti. \
Inoltre, nel makefile sotto la voce `esegui` e' possibile modificare il numero di core da usare nella simulazione.

Una volta eseguito il codice in questo modo, il codice richiedera', attraverso input da tastiera,  se si vuole eseguire una simulazione con migrazioni o senza: in questo modo sapra' se permettere le migrazioni o no e dove andare a depositare i risultati della simulazione.

\
Inoltre, il comando 
```
make clear
```
permette di cancellare tutti i file oggetto e l'eseguibile main.exe, mentre
```
make clear_results_&borders_situation&
```
da' la possibilita' di cancellare i risultati della simulazione desiderata nella cartella corrispondente sostituendo &borders_situation& con open o closed.
