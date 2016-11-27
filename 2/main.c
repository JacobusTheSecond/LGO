#include <stdio.h>
#include <stdlib.h>
#include "read_lp.h"

int writeMatrix(double** A, int m, int n){
	//OUTPUTTEST A_B TODO remove
	for(int i=0;i<m;++i){
		printf("[");
		for(int j=0;j<n;++j){
			printf("%3g",A[i][j]);
			if(j != n-1){
				printf(",");
			}
		}
		printf("]\n");
	}
	return printf("\n");
}

int writeVector(double* b, int m){
	for(int i=0;i<m;++i){
		printf("[%3g]\n",b[i]);
	}
	return printf("\n");
}

int  Gauss_solve(double** A, double* b, int m, int n){
	//A ist Matrix mit m Zeilen und n Spalten
	
	//Diagonalisieren
	for(int j=0;j<m;++j){
		int isnonzero = 1;
		if(A[j][j] == 0){
			isnonzero = 0;
			for(int i = j;i<m; ++i){
				if(A[i][j]!=0){
					double* swap = A[j];
					A[j] = A[i];
					A[i] = swap;
					isnonzero = 1;
					break;
				}
			}
		}
		if(isnonzero == 1){
			//Nach unten abziehen
			for(int k=0;k<m;++k){
				if(k!=j){
					double lambda = A[k][j]/A[j][j];
					//Matrix
					for(int l=j;l<n;++l){
						A[k][l] -= lambda*A[j][l];
					}
					//Vektor
					b[k] -= lambda*b[j];

					printf("reduced to:\n");
					writeMatrix(A,m,n);
					printf("=\n");
					writeVector(b,m);
				}
			}
		}
	}

	//normaliseren
	for(int i=0;i<m;++i){
		if(A[i][i]!=0){
			b[i]/= A[i][i];
			A[i][i]=1;
		}else{
			if(b[i]!=0){
				return 0;
			}
		}
	}
	writeMatrix(A,m,n);
	writeVector(b,m);
	return 1;
}


int SimplexAlgorithm(double** A, double* b, double* c,int m,int n){
	//A ist eine Matrix mit m Zeilen und n Spalten
	if(n==0){
		return 0;
	}

	//Gauss if Simplex not needed
	if(n<=m){
		int solvable = Gauss_solve(A,b,m,n);
		if(solvable){
			printf("result = \n");
			writeVector(b,m);
		}else{
			printf("INFEASIBLE (Gauss)\n");
		}
		return 0;
	}

	int* B = malloc(m*sizeof(int));
	//2
	//Compute feasible Basis
	B[0]=0;
	double** A_B;
	int Bsize = 1;
	printf("Basis is beeing generated\n");

	for(int add=1;add<n;++add){
		if(Bsize==m){
			break;
		}
		//Check for linear independence
		
		//fill AB with A_B and A_i
		A_B = malloc(m*sizeof(double*));
		for(int b=0;b<m;++b){
			A_B[b]=malloc((Bsize+1)*sizeof(double));
			for(int lines=0;lines<Bsize+1;++lines){
				if(lines == Bsize){
					A_B[b][lines]=A[b][add];
				}else{
					A_B[b][lines]=A[b][B[lines]];
				}
			}
		}
		
		printf("trying Basis:\n");
		writeMatrix(A_B,m,Bsize+1);

		//Gauss Elimination auf A_B
		for(int i=0;i<Bsize+1;++i){
			for(int j=0;j<m;++j){
				if(A_B[j][i]!=0){
					double lambda = A_B[j][i];

					//Spalte normieren
					for(int k = j;k<m;++k){
						A_B[k][i] = A_B[k][i] / lambda;
					}
					//von Zeilen unter j abziehen
					if(j!=m-1){
						for(int k=j+1;k<m;++k){
							lambda = A_B[k][i];
							for(int l=i;l<Bsize+1;++l){
								A_B[k][l] -= lambda*A_B[j][l];
							}
						}
					}
					//von Spalten rechts von j azbziehen
					if(i!=Bsize){
						for(int k=i+1;k<Bsize+1;++k){
							A_B[j][k]=0;
						}
					}
				}
				printf("reduced to:\n");
				writeMatrix(A_B,m,Bsize+1);
			}
		}

		int addi = 1;
		//leere Spalte -> i wird nicht zu B hinzugefügt
		for(int i=0;i<Bsize+1;++i){
			int nonnull = 0;
			for(int j=0;j<m;++j){
				if(A_B[j][i] != 0){
					nonnull = 1;
				}
			}
			if(nonnull==0){
				addi = 0;
				break;
			}
		}
		if(addi==1){
			B[Bsize]=add;
			++Bsize;
		}
	}
	for(int i=0;i<m;++i){
		free(A_B[i]);
	}
	free(A_B);

	//OUTPUT B
	printf("Basis is [");
	for(int i=0;i<Bsize;++i){
		printf("%d",B[i]);
		if(i!=Bsize-1)printf(",");
	}
	printf("]\n");
	
	//2
	if(Bsize != m){
		printf("INFEASIBLE, because Basis too small\n");
		free(B);
		return 0;
	}
	while(1){	
		//3
		//set N
		int Bindex = 0;
		int Nindex = 0;
		int* N = malloc(sizeof(int)*(n-m));
		for(int i=0;i<n;++i){
			if(B[Bindex]==i){
				++Bindex;
			}else{
				N[Nindex]=i;
				++Nindex;
			}
		}
	
		//OUTPUT N
		printf("Therefore N is [");
		for(int i=0;i<n-m;++i){
			printf("%d",N[i]);
			if(i!=n-m-1)printf(",");
		}
		printf("]\n");
	
		double* p = malloc(m*sizeof(double));
		//3.5 & 4
		//Compute Basic solution / Q
		
		double** Q = malloc(m*sizeof(double*));
		for(int i=0;i<m;++i){
			Q[i] = malloc(n*sizeof(double));
			for(int j=0; j<m;++j){
				Q[i][j] = A[i][B[j]];
			}
			for(int j = m; j<n; ++j){
				Q[i][j] = A[i][N[m-j]];
			}
			p[i] = b[i];
		}
		//p ist zulässige Basislösung
		//Q[m] bis Q[n-1] ist Q aus dem Alg
		Gauss_solve(Q,p,m,n);
		
		//calculate r
		double* rstrich = malloc(sizeof(double)*n);
		for(int i=0;i<n;++i){
			rstrich[i]=c[i];
		}
	
		int Biterator = 0;
		for(int i=0;i<n;++i){
			if(i == B[Biterator]){
				for(int j=m;j<n;++j){
					//TODO faktor überprüfen
					rstrich[N[j-m]]+=rstrich[B[Biterator]]*Q[Biterator][j];
				}
				
				++Biterator;
			}
		}
		writeMatrix(&rstrich,1,n);
		double* r = malloc(sizeof(double)*(n-m));
	
		int allnegative = 1;
		for(int i=0;i<n-m;++i){
			r[i] = rstrich[N[i]];
			if(r[i]>0){
				allnegative = 0;
			}
		}
		free(rstrich);
		printf("r ist :\n");
		writeVector(r,n-m);
		
		//5
		if(allnegative == 1){
			double * x = calloc(n,sizeof(double));
			for(int i=0;i<m;++i){
				x[B[i]] = p[i];
				free(Q[i]);
			}	
			printf("result = \n");
			writeVector(x,n);
			free(B);
			free(N);
			free(p);
			free(r);
			free(Q);
			free(x);
			return 0;
		}
		//6
		//Blands rule pt. 1 
		int alpha;
		for(int i=0;i<n-m;++i){
			if(r[i]>0){
				alpha = N[i];
				break;
			}
		}
	
		//7 & 8
		//BLands rule pt. 2
		int beta=-1;
		for(int i=0;i<m;++i){
			if(Q[i][m+alpha]<0){
				if(beta==-1){
					beta = B[i];
				}else{
					if(p[beta] / Q[beta][m+alpha] < p[i] / Q[i][m+alpha]){
						beta = B[i];
					}
				}
			}
		}

		if(beta==-1){
			printf("UNBOUNDED, because no negative entry in r\n");
			return 0;
		}
	
		//9
		int* newB = malloc(m*sizeof(int));
		int newBiterator = 0;
		printf("new Basis is [");
		for(int i = 0;i<m;++i){
			if(B[i] != beta){
				newB[newBiterator]=B[i];
				++newBiterator;
			}
			if(i==alpha){
				newB[newBiterator]=alpha;
				++newBiterator;
			}
		}
		for(int i=0;i<m;++i){
			printf("%d",newB[i]);
			if(i!=m-1)printf(",");
		}
		printf("]\n");
		free(B);
		B=newB;
	}
}



int main(){
	//Eingabe für unterschiedliche Instanzen
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

	int result = SimplexAlgorithm(A,b,c,m,n);
	free(b);
	free(c);
	for(int i = 0;i<m;++i){
		free(A[i]);
	}
	free(A);
	return EXIT_SUCCESS;
}
