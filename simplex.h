#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <stdlib.h>
#include <stdio.h>
#define N 4
#define ERROR 1.0e-6
#define PRINT
#define ABS(a) ({a < 0 ? -a : a;})

float ** simplex_procedure(float * X, int ** B, int n);
void print_vector(char * name, float * vector, int len);
void print_matrix(char * name, float ** matrix, int nb_line, int nb_column);

#endif