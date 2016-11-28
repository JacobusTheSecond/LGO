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
	return 0;
}

int writeVector(double* b, int m){
	for(int i=0;i<m;++i){
		printf("[%3g]\n",b[i]);
	}
	return 0;
}

int contains(int* Array,int size, int check){
	for(int i=0;i<size;++i){
		if(Array[i]==check)return 1;
	}
	return 0;
}

int  Gauss_solve(double** A, double* b, int m, int n){
	//A ist Matrix mit m Zeilen und n Spalten
	
	//Diagonalisieren
	for(int j=0;j<((m<n)?m:n);++j){
		int isnonzero = 1;
		if(A[j][j] == 0){
			isnonzero = 0;
			for(int i = j;i<m; ++i){
				if(A[i][j]!=0){
					double* swap = A[j];
					double bswap = b[j];
					A[j] = A[i];
					b[j] = b[i];
					A[i] = swap;
					b[i] = bswap;
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

					printf("\nreduced to:\n");
					writeMatrix(A,m,n);
					printf("=\n");
					writeVector(b,m);
				}
			}
		}
	}

	//normaliseren
	int unbound = 0;
	for(int i=0;i<m;++i){
		if(i<n||A[i][i]!=0){
			double lambda = A[i][i];
			b[i]/=lambda;
			for(int j=0;j<n;++j){
				if(A[i][j] != 0)A[i][j]/=lambda;
				if(j!=i){
					if(A[i][j]!=0){
						unbound = 1;
					}
				}
			}
		}else{
			//not feasible
			if(b[i]!=0){
				return 1;
			}
		}
	}
	printf("\nnormalised to:\n");
	writeMatrix(A,m,n);
	printf("=\n");
	writeVector(b,m);

	//slackvariable with no boundries exists -> unbound
	if(unbound ){
		return 2;
	}
	return 0;
}


int SimplexAlgorithm(double** A, double* b, double* c,int m,int n){


	//A ist eine Matrix mit m Zeilen und n Spalten
	
	//leere Matrix wird nicht akzeptiert
	if(n==0){
		printf("empty Matrix!\n");
	}
	
	//Gauss if Simplex not needed
	if(n<=m){
		int solvable = Gauss_solve(A,b,m,n);
		if(solvable==0){
			printf("result = \n");
			writeVector(b,m);
		}else if(solvable==1){
			printf("\nINFEASIBLE (Gauss)\n");
		}else{
			printf("\nUNBOUNDED (Gauss)\n");
		}
		return 0;
	}

	int* B = malloc(m*sizeof(int));
	//2
	//Compute feasible Basis
	B[0]=0;
	double** A_B;
	int Bsize = 1;
	printf("----------\nBasis is beeing generated\n----------\n");

	//durch {1,...,n} iterieren um Kandidaten für die Basis zu erhalten
	for(int add=1;add<n;++add){

		//falls m linear unabhängige Vektoren gefunden wurden, ist B eine feasible Basis
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
		
		printf("\ntrying Basis:\n");
		writeMatrix(A_B,m,Bsize+1);

		//modifizierte Gauss Elimination auf A_B um auf lineare unabhängigkeit zu testen
		for(int i=0;i<Bsize+1;++i){
			for(int j=0;j<m;++j){
				if(A_B[j][i]!=0){
					double lambda = A_B[j][i];

					//Spalte normieren 
					//ist das selbe wie v -> lambda*v
					for(int k = j;k<m;++k){
						A_B[k][i] = A_B[k][i] / lambda;
					}
					//von Zeilen unter j abziehen
					//Gauseliminierung
					if(j!=m-1){
						for(int k=j+1;k<m;++k){
							lambda = A_B[k][i];
							for(int l=i;l<Bsize+1;++l){
								A_B[k][l] -= lambda*A_B[j][l];
							}
						}
					}
					//von Spalten rechts von j azbziehen
					//w = w-v. geht weil pivot spalten
					if(i!=Bsize){
						for(int k=i+1;k<Bsize+1;++k){
							A_B[j][k]=0;
						}
					}
				}
				printf("\nreduced to:\n");
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

		//falls die Vektoren linearely independent sind wird der Index zu B hinzugefügt
		if(addi==1){
			B[Bsize]=add;
			++Bsize;
		}
		for(int i=0;i<m;++i){
			free(A_B[i]);
		}
		free(A_B);
	}

	//OUTPUT B
	printf("\nBasis B is [");
	for(int i=0;i<Bsize;++i){
		printf("%d",B[i]);
		if(i!=Bsize-1)printf(",");
	}
	printf("]\n");
	
	//2
	if(Bsize != m){
		printf("\nINFEASIBLE, because Basis too small\n");
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
		printf("\nTherefore N is [");
		for(int i=0;i<n-m;++i){
			printf("%d",N[i]);
			if(i!=n-m-1)printf(",");
		}
		printf("]\n");
		double* p = malloc(m*sizeof(double));
		
		printf("-----------\nCalculating Basic solution and Q\n----------\n");
		//3.5 & 4
		//Compute Basic solution / Q
		double** Q = malloc(m*sizeof(double*));
		for(int i=0;i<m;++i){
			Q[i] = malloc(n*sizeof(double));
			for(int j=0; j<m;++j){
				Q[i][j] = A[i][B[j]];
			}
			for(int j = m; j<n; ++j){
				Q[i][j] = A[i][N[j-m]];
			}
			p[i] = b[i];
		}

		//Q ist A_B in den ersten m Spalten, und A_N in den anderen n-m
		//Gaußeliminierung auf Q,b macht daraus ( 1_m | -Q_x_N),p, sodass x_B = p + Q_x_N gilt
		writeMatrix(Q,m,n);
		printf("=\n");
		writeVector(p,m);

		Gauss_solve(Q,p,m,n);
		
		//calculate r
		double z_0 = 0;
		double* rstrich = malloc(sizeof(double)*n);
		for(int i=0;i<n;++i){
			rstrich[i]=c[i];
		}

		//macht aus rstrich (r_i in rstrich_N_i) und (0 in rstrich_B_i)
		for(int i=0;i<m;++i){
			if(rstrich[B[i]]!=0){
				for(int j=m;j<n;++j){
					// "-", weil Q anderes Vorzeichen hat
					rstrich[N[j-m]]-=rstrich[B[i]]*Q[i][j];	
				}
				z_0 += rstrich[B[i]]*p[i];
				rstrich[B[i]]=0;
			}
		}

		//liest r aus rstrich
		double* r = malloc(sizeof(double)*(n-m));
		int allnegative = 1;
		for(int i=0;i<n-m;++i){
			r[i] = rstrich[N[i]];
			if(r[i]>0){
				allnegative = 0;
			}
		}
		free(rstrich);	
		
		//T(B) printen
		printf("\n----------\nSimplex Tableau:\n----------\n\n");
		for(int i=0;i<m;++i){
			printf("x_%d = %3g + ",B[i],p[i]);
			for(int j=0;j<n-m;++j){
				//Q[m] bis Q[n] ist -Q_x_N, deswegen -1*...
				printf("%3g*x_%d",(Q[i][m+j]!=0?-1:1)*Q[i][m+j],N[j]);
				if(j != n-m-1)printf(" + ");
			}
			printf("\n");
		}
		printf("-----------------------------------\n");
		//z0 ist nicht wichtig für meine Regel, deswegen nicht berechnet
		printf("z   = %3g + ",z_0);
		for(int i=0;i<n-m;++i){
			printf("%3g*x_%d",r[i],N[i]);
			if(i != n-m-1)printf(" + ");
		}
		printf("\n\n");


		//5
		if(allnegative == 1){
			double * x = calloc(n,sizeof(double));
			for(int i=0;i<m;++i){
				x[B[i]] = p[i];
				free(Q[i]);
			}	
			printf("\nresult = \n");
			writeVector(x,n);
			free(B);
			free(N);
			free(p);
			free(r);
			free(Q);
			free(x);
			return 0;
		}

		printf("\n----------\nChoosing Indices to swap\n----------\n");
		//6
		//Blands rule pt. 1, wählt alpha minimal, sodass r_alpha > 0
		int alpha;
		for(int i=0;i<n-m;++i){
			if(r[i]>0){
				alpha = i;
				break;
			}
		}
		printf("alpha = %d nach Bland's Regel\n",N[alpha]);	
		//7 & 8
		//BLands rule pt. 2, wählt beta minimal, sodass Q_i,alpha < 0 und Q_beta_alpha = max(Q_i,alpha | Q_i,alpha < 0)
		int beta=-1;
		for(int i=0;i<m;++i){
			if(Q[i][m+alpha]>0){
				if(beta==-1){
					beta = i;
				}else{
					//Q[m] bis Q[n] ist -Q_x_N
					if(p[beta]/(-1*Q[beta][m+alpha]) < p[i]/(-1*Q[i][m+alpha])){
						beta = i;
					}
				}
			}
		}

		//wurde kein beta gefunden ist das LP unbounded
		if(beta==-1){
			printf("\nUNBOUNDED, because no negative entry in Q_x_%d for all i\n",N[alpha]);
			free(B);
			free(N);
			for(int i=0;i<m;++i){
				free(Q[i]);
			}
			free(Q);
			free(p);
			free(r);
			return 0;
		}
		printf("beta = %d nach Bland's Regel\n",B[beta]);
		
		printf("\n----------\nGenerating new Base\n----------\n");
		//9
		//B = B\beta v alpha
		int* newB = malloc(m*sizeof(int));
		int newBiterator = 0;
		printf("\nnew Basis is [");
		for(int i = 0;i<n;++i){
			if((contains(B,m,i)&&i!=B[beta]) || i == N[alpha]){
				newB[newBiterator]=i;
				++newBiterator;
			}
		}
		//OUTPUT B
		for(int i=0;i<m;++i){
			printf("%d",newB[i]);
			if(i!=m-1)printf(",");
		}

		//free everything
		printf("]\n");
		free(B);
		B=newB;
		free(N);
		for(int i=0;i<m;++i){
			free(Q[i]);
		}
		free(Q);
		free(p);
		free(r);
		//goto 3
	}
}



int main(int argc, char* argv[]){
	char* in;
	if(argc==1){
		//Eingabe für unterschiedliche Instanzen
		in = malloc(256);
		printf("Programm zum Lösen eines Linearen Programmes mit dem Simplex Algorithmuses.\nDatei-Pfad der zu analysierenden Instanzen eingeben:\n");
		scanf(" %s",in);
		//Deklaration aller notweniden Felder
	}else{
		in = argv[1];
	}
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
	if(argc==1)free(in);
	return EXIT_SUCCESS;
}
