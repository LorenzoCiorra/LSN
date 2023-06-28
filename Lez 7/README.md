# Istruzioni di Funzionamento

Per andare ad eseguire i codici e' necessario compilare con  
```
make
```
ed eseguire il codice esplicitando la fase di interesse con 
```
make esegui_&fase&
```
rimpiazzando &fase& con solido, liquido o gas.
In questo modo, il makefile prima di eseguire il codice automaticamente pulira' i risultati della simulazione precedente riguardanti quella fase. Pulire gli output prima di far ripartire una simulazione e' necessario visto che i file verranno aperti in append e, quindi, se non si pulissero i file di output conterrebbero misure di simulazioni differenti.

Una volta eseguito il codice in questo modo, il codice chiedera' due input da tastiera:
- il primo per sollecitare l'utente al controllo della pulizia degli output precedenti;
- il secondo per chiedere di inserire nuovamente la fase da voler simulare: in questo modo il codice leggera' il file di input corretto, quello presente nella cartella associata alla fase scelta, e sapra' in che cartella mettere i risultati.

\
Inoltre, il comando 
```
make clear
```
permette di cancellare tutti i file oggetto e l'eseguibile main.exe, mentre
```
make clear_results_&fase&
```
da' la possibilita' di cancellare i risultati della simulazione della fase desiderata nella cartella corrispondente sostituendo &fase& con solido, liquido o gas.
