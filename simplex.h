#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <time.h>
#include "matrix.h"
#define ERROR 1.0e-6
// #define PRINT

typedef struct Simplex {
	int nb_var;
	int nb_row;
	int nb_column;

	float * A;
	float * b;
	float * c;

	float obj;
	int * basis;
} Simplex;

void simplex_procedure(float * X, int ** B, int n);

#endif