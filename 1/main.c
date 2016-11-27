#include "read_lp.h"
#include "fouriermotzkin.h"

int main(){

	//Eingabe für unterschiedliche Instanzen
	//
	//TODO besser, vielleicht über Konsolen-Eingabe des Dateipfades 
	char* in = malloc(256);
	printf("Programm zur Bestimmung der Lösbarkeit eines linearen Programmes\nDatei-Pfad der zu analysierenden Instanzen eingeben:\n");
	scanf(" %s",in);

	//Deklaration aller notweniden Felder
	double** A;
	double* b;
	double* c;
	int m;
	int n;
	if(read_LP(in,&m,&n,&A,&b,&c)){printf("Probleme beim Datei auslesen\n");return EXIT_FAILURE;}

	//Da sowohl m sich ändert, als auch die aller erste Zeilenanzahl gebrauch wird, wird es hier nochmal gespeichert
	int initialm = m;

	//nocheinmal eine kurze Ausgabe, um zu bestätigen was engelesen wurde. Copy-Pasted aus test_it(...);
	printf("INPUT:\n----------------------\n");

	//Vektorausgabe
	printf("n: %d\nm: %d\ntranspose(b) = [",n,m);
	for (int i = 0; i < m; i++) {
		printf("%g%s", b[i], (i == m-1 ? "]\n\n" : ", "));
	}

	//Matrixausgabe
	for (int i = 0; i < m; i++) {
		printf("A[%d] = [", i);
		for (int j = 0; j < n; j++) {
			printf("%g%s", A[i][j], (j == n-1 ? "]\n" : ", "));
		}
	}

	//Hier gehts um die Wurst (also Erfüllbarkeit)
	printf("----------------------\n\nBERECHNUNG:\n----------------------\n");
	
	//Eine Matrix, inder die Linearkombinationen, durch die die Zeilen abhängig von den aller ersten zeilen erzeugt werden
	//Wird gebraucht, um das Zertifikat am Ende auszugeben
	double** delta = malloc(m*sizeof(double*));

	//Am Anfang sind die Linearkombinationen der Zeilen einfach die Einheitsmatrix 1_m
	for(int i=0;i<m;++i){
		delta[i] = calloc(m,sizeof(double));
		delta[i][i]=1;
	}

	//Fourier-Motzkin-Elimination wird solange ausgeführt, bis die Matrix leer ist
	int zaehler = 0;
	while(n>0){
		++zaehler;
		
		//Ein Fourier-Motzkin-Eliminations-Aufruf, der selbst Ausgaben macht, über was passiert
		FM(&A,&b,&m,&n,&delta,initialm);
	}

	//Ausgabe darüber, ob eine Lösung existiert, oder falls nicht der Gegenbeweis nach Farakas
	printf("%d loops\n----------------------\n\nAUSGABE:\n----------------------\n",zaehler);

	//Überprüfung und Ausgabe, ob einer der Einträge aus dem resultierenden Vektor b negativ ist
	int flag = 0;
	for(int i=0;i<m;++i){
		//Falls ein negativer Eintrag gefunden wird, wird dies entsprechend mit Farakas-Zertifikat ausgegeben
		if(b[i]<0){
			printf("keine Lösung existiert\n");
		
			printf("transpose(zertifikat) = [");
			for (int j = 0; j < initialm; j++) {
				printf("%g%s", delta[i][j], (j == initialm-1 ? "]\n" : ", "));
			}
		
			flag = 1;	
			break;
		}
	}

	//Ansonsten wird ausgegeben, dass Lösungen existieren
	if(flag==0)printf("Lösung existiert\n");
	printf("----------------------\n");

	//Entfernen, von allen übrigen Feldern
	release_memory(m,A,b,c);
	release_memory(m,delta,(void*)0,(void*)0);
}

