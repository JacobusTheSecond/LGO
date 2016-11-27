// |  a  b  c |   |x|   |i|					|	  | 1 0 0 | 	
// | -d  e  f | * |y| ≤ |j|  mit a,d > 0 (nur um zu verdeu-	| delta : | 0 1 0 | repräsentiert die Lineardarstellung der Zeilen der 
// |  0  g  h |   |z|   |k|	tlichen, was passiert, wenn	|	  | 0 0 1 | aktuellen Matrix durch Zeilen aus der Anfangs-
//				wenn eins <0 und eins >0 ist)	|	            matrix AFM. Die erste Zeile ist | 1 0 0 |. Das bedeuted,
// =>   ax + by + cz ≤ i   <=>   x ≤ -(b/a)y - (c/a)z + i/a	|		    dass die erste Zeile in der aktuellen Matrix erzeugt wurde
// =>  -dx + ey + fz ≤ j   <=>   x ≥  (e/d)y + (f/d)z - j/d	|		    durch 1*[erste Zeile der AFM] + 0*[zweite...] + ...
//								|
// =>  (e/d)y + (f/d)z - j/d ≤ -(b/a)y - (c/a)z + i/a		| Um die erste Spalte zu eliminieren, wird die erste Zeile I*1/a genommen,
//<=>  (e/d + b/a)y + (f/d + c/a)z ≤ j/d + i/a			| und die zweite II*(1/-d) -> die neue Zeile ist 1/a*I - 1/d*II + 0*III
//								|  => 
//<=> | (e/d + b/a) (f/d + c/a) |   |y|   | j/d + i/a |		| delta : | (1/a) -(1/d) 0 |
//    |      g           h      | * |z| ≤ |     k     |		|	  |   0     0    1 |

#include "fouriermotzkin.h"
#include <stdlib.h>
#include <stdio.h>
typedef struct intarray intarray;

struct intarray { 
	int* data;
	int length;
};

/* Fourier-Motzkin-Elimination, die in A, b, m und n die neuen werte schreibt.
 * delta wird benutzt, um sich die linearkombinationen in der resultierenden
 * Matrix zu merken, um ein Zertifikat nach Farakas Lemma auszugeben */
int FM(double*** A, double** b, int *m, int *n,double*** delta, int initialm){

	// Im folgenden ist A immer eine Matrix mit m Zeilen und n Spalten

	// Z, P, N sind Indexmengen die alle Zeilen mit dem Ersten Eintrag x =0/>0/<0 aufzählen
	// upperbound für die Anzahl der Zeilen mit <condition> ist <= m
	intarray* _Z = malloc(sizeof(intarray));
	intarray* _P = malloc(sizeof(intarray));
	intarray* _N = malloc(sizeof(intarray));

	_Z->data = malloc(*m*sizeof(int));
	_P->data = malloc(*m*sizeof(int));
	_N->data = malloc(*m*sizeof(int));

	_Z->length = 0;
	_P->length = 0;
	_N->length = 0;

	//Z,P,N mit information füllen, und Z_length,P_length,N_length entsprechend setzen
	for(int i=0;i<*m;++i){
		if((*A)[i][0]>0){
			_P->data[_P->length]=i;
			++_P->length;
		}else if((*A)[i][0]<0){
			_N->data[_N->length]=i;
			++_N->length;
		}else{
			_Z->data[_Z->length]=i;
			++_Z->length;
		}
	}

	

	//neue Matrix, delta und Vektor
	double** _RA = malloc((_Z->length+_N->length*_P->length)*sizeof(double*));
	double* _Rb = malloc((_Z->length+_N->length*_P->length)*sizeof(double));
	double** _Rdelta = malloc((_Z->length+_N->length*_P->length)*sizeof(double*));

	//die ersten po*ne Zeilen füllen

	//i iteriert durch die P-Indexmenge
	for(int i=0;i<_P->length;++i){
		printf("i:%d,_P[i]:%d\n" , i , _P->data[i]);
	
		//j iteriert durch die N-Indexmenge
		for(int j=0;j<_N->length;++j){
			printf("    j:%d,_N[j]:%d\n",j,_N->data[j]);

			//in dem double* array müssen double arrays mallocd werden
			_RA[i*_N->length+j] = malloc((*n-1)*sizeof(double));
			_Rdelta[i*_N->length+j] = malloc(initialm*sizeof(double));

			//k iteriert durch n-1 Spalten
			for(int k=1;k<*n;++k){

				//Rechnung, um den k-ten Koeffizienten in der i*ne+j -ten Spalte zu berechnen, durch Kombination der
				//iten Zeile aus der Indexmenge _P und der Indexmenge _N wie in der Formel am Anfang
				_RA[i*_N->length+j][k-1]  = ((*A)[_P->data[i]][k] / (*A)[_P->data[i]][0]); 
				_RA[i*_N->length+j][k-1] -= ((*A)[_N->data[j]][k] / (*A)[_N->data[j]][0]);

				//Ausgabe
				printf("        k:%d\n",k);
				printf("            R[%d][%d] = %g = ",i*_N->length+j,k-1,_RA[i*_N->length+j][k-1]);
				printf("(%g / %g) - (%g / %g)\n",(*A)[_P->data[i]][k],(*A)[_P->data[i]][0],(*A)[_N->data[j]][k],(*A)[_N->data[j]][0]);

			}

			//Linearkombinationen der gerade berechneten Zeile updaten, mit den linearkombinationen der verwendeten
			//Zeilen aus der Alten Matrix multipliziert mit dem jeweiligen Sakalar
			for(int l=0;l<initialm;++l){
				_Rdelta[i*_N->length+j][l]  = (*delta)[_P->data[i]][l] / (*A)[_P->data[i]][0];
				_Rdelta[i*_N->length+j][l] -= (*delta)[_N->data[j]][l] / (*A)[_N->data[j]][0];
			}
			
			//Rechnung um den i*ne+j -ten Eintrag des Vektors zu berechnen
			_Rb[i*_N->length+j]  = ((*b)[_P->data[i]] / (*A)[_P->data[i]][0]);
			_Rb[i*_N->length+j] -= ((*b)[_N->data[j]] / (*A)[_N->data[j]][0]);
			
			//Ausgabe
			printf("        b[%d] = %g = ",i*_N->length+j,_Rb[i*_N->length+j]);
			printf("(%g / %g) - (%g / %g)\n",(*b)[_P->data[i]],(*A)[_P->data[i]][0],(*b)[_N->data[j]],(*A)[_N->data[j]][0]);
		}
	}

	//die restlichen ze Zeilen füllen
	int _Moffset = _P->length*_N->length;

	//i iteriert durch die Z-Indexmenge
	for(int i = 0;i<_Z->length;++i){
		printf("i:%d,Z[i]:%d\n",i,_Z->data[i]);

		_RA[_Moffset+i] = malloc((*n-1)*sizeof(double));
		_Rdelta[_Moffset+i] = malloc(initialm*sizeof(double));
		
		//k iteriert durch n-1 Spalten
		for(int k=1;k<*n;++k){

			//Ausgabe
			printf("    k:%d\n",k);
			printf("        RA[%d][%d] = %g\n",_Moffset+i,k-1,(*A)[_Z->data[i]][k]);
			
			//Überträgt eine Zeile die im ersten Eintrag eine 0 hat in die neue matrix, aber ohne die 0
			_RA[_Moffset+i][k-1] = (*A)[_Z->data[i]][k];
		}

		//Überträgt die Linearkombinationen aus der alten in die neue Matrix
		for(int l=0;l<initialm;++l){
			_Rdelta[i+_Moffset][l] = (*delta)[_Z->data[i]][l];
		}

		//Ausgabe
		printf("    Rb[%d] = %g\n",_Moffset+i,(*b)[_Z->data[i]]);
		
		//Überträgt den Eintrag aus b in das neue b
		_Rb[_Moffset+i] = (*b)[_Z->data[i]];
	}

	//free'd alle nichtmehr benötigten Variablen/Felder, und ändert A,b,m und n auf die neuen Werten
	release_memory(*m,*A,*b,(void*)0);
	release_memory(*m,*delta,(void*)0,(void*)0);
	*delta = _Rdelta;
	*A = _RA;
	*b = _Rb;
	--*n;
	*m = _Z->length+_N->length*_P->length;
	free(_Z->data);
	free(_Z);
	free(_P->data);
	free(_P);
	free(_N->data);
	free(_N);
	return EXIT_SUCCESS;
}
