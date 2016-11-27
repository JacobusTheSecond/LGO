#pragma once
#include <stdio.h>
#include <stdlib.h>

void release_memory(int       m,
                    double ** A,
                    double *  b,
                    double *  c);

int read_LP(const char * filename,
            int *        m,
            int *        n,
            double ***   A,
            double **    b,
            double **    c);

void test_it(const char * lp_file);
