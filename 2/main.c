#include <stdio.h>
#include <stdlib.h>
#include "read_lp.h"

int writeMatrix(double** A, int m, int n){
	//OUTPUTTEST A_B TODO remove
	for(int i=0;i<m;++i){
		for(int j=0;j<n;++j){
			printf("%lf ",A[i][j]);
		}
		printf("\n");
	}
	return printf("\n\n");
}

int  Gauss_solve(double** A, double* b, int m, int n){
	//A ist Matrix mit m Zeilen und n Spalten
	if(m>n){
		return 0;
	}
	//Diagonalisieren
	printf("%d %d\n",m,n);
	for(int j=0;j<m;++j){
		printf("j = %d\n",j);
		int isnonzero = 1;
		printf(" %lf\n",A[j][j]);
		if(A[j][j] == 0){
			isnonzero = 0;
			for(int i = j;i<m; ++i){
				printf(" i = %d\n",i);
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
					printf(" k = %d\n",k);
					double lambda = A[k][j]/A[j][j];
					//Matrix
					for(int l=j;l<n;++l){
						A[k][l] -= lambda*A[j][l];
					}
					//Vektor
					b[k] -= lambda*b[j];
					writeMatrix(A,m,n);
					writeMatrix(&b,1,m);
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
	writeMatrix(&b,1,m);
	return 1;
}


int SimplexAlgorithm(double** A, double* b, double* c,int m,int n){
	//A ist eine Matrix mit m Zeilen und n Spalten
	if(n==0){
		return 0;
	}
	if(n<=m){
		//TODO Gauss
		return 0;
	}
	int* B = malloc(m*sizeof(int));
	
	//2
	//Compute feasible Basis
	B[0]=0;
	double** A_B;
	int Bsize = 1;
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

		writeMatrix(A_B,m,Bsize+1);

		printf("m = %d, n = %d\n",m,n);
		//Gauss Elimination auf A_B
		for(int i=0;i<Bsize+1;++i){
			printf("i = %d\n",i);
			for(int j=0;j<m;++j){
				printf(" j = %d\n",j);
				if(A_B[j][i]!=0){
					printf("  !=0\n");
					double lambda = A_B[j][i];
					//Spalte normieren
					for(int k = j;k<m;++k){
						printf("   k = %d\n",k);
						A_B[k][i] = A_B[k][i] / lambda;
						writeMatrix(A_B,m,Bsize+1);
					}
					//von Zeilen unter j abziehen
					if(j!=m-1){
						printf("  != m-1\n");
						for(int k=j+1;k<m;++k){
							printf("   k = %d\n",k);
							lambda = A_B[k][i];
							for(int l=i;l<Bsize+1;++l){
								printf("    l = %d\n",l);
								A_B[k][l] -= lambda*A_B[j][l];
								writeMatrix(A_B,m,Bsize+1);

							}
						}
					}
					//von Spalten rechts von j azbziehen
					if(i!=Bsize){
						printf("  != Bsize\n");
						for(int k=i+1;k<Bsize+1;++k){
							printf("   k = %d\n",k);
							A_B[j][k]=0;
							writeMatrix(A_B,m,Bsize+1);
						}
					}
				}
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
	printf("B = ");
	for(int i=0;i<Bsize;++i){
		printf("%d ;",B[i]);
	}
	printf("\n");
	
	//2
	if(Bsize != m){
		printf("INFEASIBLE\n");
		free(B);
		return 0;
	}
	
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
	printf("N = ");
	for(int i=0;i<n-m;++i){
		printf("%d ;",N[i]);
	}
	printf("\n");

	double* p = malloc(m*sizeof(double));
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
	
	//calculate r?
	


	//TODO free Q


	free(B);
	free(N);
	return 0;
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
	return EXIT_SUCCESS;
}
