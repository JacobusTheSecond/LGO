#include "read_lp.h"
void release_memory(int       m,
                    double ** A,
                    double *  b,
                    double *  c)
{
   if (A) {
      int i;
      for (i = 0; i < m; i++) {
         if (A[i]) {
            free(A[i]);
         } else {
            /* If A[i] is NULL, then the next entries won't have been
               allocated and their pointer values will be undefined. */
            break;
         }
      }
      free (A);
   }

   if (b) {
      free(b);
   }
   if (c) {
      free(c);
   }
}

/* The function reads an LP instance from filename. The file
   format is expected to be exactly as in the problem specification.
   On return, *m and *n will be the number of rows and columns,
   respectively; A, b and c will be the the matrix, right-hand side
   vector and objective function vector, respectively. Note that
   the order of the data in the input file is c, then b, then A.
   Memory will be allocated for A, b and c; it is the caller's
   responsibility to release the memory later.
*/
int read_LP(const char * filename,
            int *        m,
            int *        n,
            double ***   A,
            double **    b,
            double **    c)
{
   FILE *fp;

   if (!(fp = fopen(filename, "r")))
   {
      fprintf(stderr, "Could not open \"%s\".\n", filename);
      return EXIT_FAILURE;
   } else {
      int i;
      int j;

      fscanf(fp, "%d %d\n", m, n);

      /* Memory allocation. */
      if ((*A = (double**) malloc(*m * sizeof(double*))) &&
          (*b = (double*)  malloc(*m * sizeof(double)))  &&
          (*c = (double*)  malloc(*n * sizeof(double)))) {
         for (i = 0; i < *m; i++) {
            if (!((*A)[i] = (double*) malloc(*n * sizeof(double)))) {
               fprintf(stderr, "Memory allocation failure.\n");
               release_memory(i - 1, *A, *b, *c);
               fclose(fp);
               return EXIT_FAILURE;
            }
         }
      } else {
         fprintf(stderr, "Memory allocation failure.\n");
         release_memory(0, *A, *b, *c);
         fclose(fp);
         return EXIT_FAILURE;
      }

      /* Copying the values into A, b and c. */
      for (j = 0; j < *n; j++) {
         fscanf(fp, "%lf", *c + j);
      }
      for (i = 0; i < *m; i++) {
         fscanf(fp, "%lf", *b + i);
      }
      for (i = 0; i < *m; i++) {
         for (j = 0; j < *n; j++) {
         /* fscanf(fp, "%lf", (*A)[i]+j); */
            fscanf(fp, "%lf", *(*A+i)+j);
         }
      }

      fclose(fp);
      return EXIT_SUCCESS;
   }
}


void test_it(const char * lp_file)
{
   int       m = 0;
   int       n = 0;
   double ** A;
   double *  b;
   double *  c;

   if (read_LP(lp_file, &m, &n, &A, &b, &c) != EXIT_SUCCESS) {

      fprintf(stderr, "Error parsing \"%s\".\n", lp_file);

   } else {

      int i;
      int j;

      printf("A has %d row%s and %d column%s.\n",
             m, (m == 1 ? "" : "s"),
             n, (n == 1 ? "" : "s"));
      printf("c = [");
      for (j = 0; j < n; j++) {
         printf("%.1lf%s", c[j], (j == n-1 ? "]\n" : ", "));
      }
      printf("transpose(b) = [");
      for (i = 0; i < m; i++) {
         printf("%.1lf%s", b[i], (i == m-1 ? "]\n" : ", "));
      }
      for (i = 0; i < m; i++) {
         printf("A[%d] = [", i);
         for (j = 0; j < n; j++) {
            printf("%5.1lf%s", A[i][j], (j == n-1 ? "]\n" : ", "));
         }
      }

      release_memory(m, A, b, c);
   }
}

/*int main(int argc, const char * argv[])
{
   if (argc < 2) {
      fprintf(stderr, "Usage:  %s  <lp file>\n", argv[0]);
      return EXIT_FAILURE;
   }
   test_it(argv[1]);
   return EXIT_SUCCESS;
}*/
