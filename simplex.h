#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <time.h>
#include <math.h>
#include "matrix.h"
#define ERROR 1.0e-2
// #define PRINT
#define MEASURE

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

float ** simplex_procedure(float * X, int ** B, int n);

#endif