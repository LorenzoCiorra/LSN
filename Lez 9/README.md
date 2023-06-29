# Istruzioni di Funzionamento 

Per andare ad eseguire i codici e' necessario compilare con  
```
make
```
ed eseguire il codice esplicitando la simulazione di interesse con 
```
make esegui_&cities-position&_&number-of-chromosomes&
```
rimpiazzando &cities-position& con s o c a seconda che si voglia fare una simulazione in cui le citta' sono disposte all'interno del quadrato o sul cerchio e sostituendo &number-of-chromosomes& con 500 o 2000 in base a quanto numerosa si vuole la popolazione di cromosomi.
In questo modo, il makefile prima di eseguire il codice automaticamente pulira' i risultati della simulazione precedente riguardanti quella situazione. Pulire gli output prima di far ripartire una simulazione e' necessario dato che i file verranno aperti in append e, quindi, se non si pulissero i file di output conterrebbero misure di simulazioni differenti. \
Per ottenere, quindi, tutti i file di risultati necessari per i plot nel jupyter notebook e' necessario eseguire 4 volte il codice.

Una volta eseguito il codice in questo modo, esso nuovamente richiedera', attraverso input da tastiera, se si vuole eseguire una simulazione con le citta' poste sulla circonferenza o nel quadrato: in questo modo sapra' dove posizionare le citta'.

\
Inoltre, il comando 
```
make clear
```
permette di cancellare tutti i file oggetto e l'eseguibile main.exe, mentre
```
make clear_results_&cities-position&_&number-of-chromosomes&
```
da' la possibilita' di cancellare i risultati della simulazione con citta' sul cerchio o nel quadrato, sostituendo a &cities-position& c o s, presenti nella cartella corrispondente al numero di cromosomi esplicitato in &number-of-chromosomes&, 500 o 2000.
