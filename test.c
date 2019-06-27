#include <time.h>
#include "io.h"
#include "simplex.h"

/* PARAMETERS */
#define NUMBER_OF_PARTITION 4
#define NUMBER_OF_NEIGHBOUR 4 // Multiple of 2
#define MEASURE
#define PRINT
/* ---------- */

int ** init_B(int n) {
	int i, j;
	int ** B = (int **) malloc(n * sizeof(int *));
	for (i = 0 ; i < n ; i++) {
		B[i] = (int *) calloc(n, sizeof(int));
	}
	for (i = 0 ; i < n-1 ; i++) {
		for (j = i+1 ; j < n ; j++) {
			if (j < i+1+NUMBER_OF_NEIGHBOUR/2) {
				B[i][j] = 1;
				B[j][i] = 1;
			}
		}
	}
	return B;
}

double average(double * vector, int len) {
	int i;
	double sum = 0;
	for (i = 0 ; i < len ; i++) {
		sum += vector[i];
	}
	return sum/(double)len;
}

double * init_X(int n) {
	int i;
	double * X = (double *) malloc(n * sizeof(double));
	for (i = 0 ; i < n ; i++) {
		X[i] = ((double)rand())/((double)RAND_MAX)*100.0;
		// X[i] = (double) i + 1;
	}
	double avg = average(X, n);
	for (i = 0 ; i < n ; i++) {
		X[i] -= avg;
	}
	return X;
}

int main() {
	time_t t;
	srand((unsigned) time(&t));	
	/* init example */
	int ** B = init_B(NUMBER_OF_PARTITION);
	double * X = init_X(NUMBER_OF_PARTITION);

	/* print example */
	#ifdef PRINT
	print_B(B, NUMBER_OF_PARTITION);
	print_vector("X", X, NUMBER_OF_PARTITION);
	#endif

/* ---------------------------------------- */
	#ifdef MEASURE
	clock_t begin = clock();
	#endif

	double ** out = simplex_procedure(X, B, NUMBER_OF_PARTITION);

	#ifdef MEASURE
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("\ntime spent is : %.9lf\n\n", time_spent);
	#endif
/* ---------------------------------------- */

	#ifdef PRINT
	print_matrix("output", out, NUMBER_OF_PARTITION, NUMBER_OF_PARTITION);
	#endif

	int i;
	for (i = 0 ; i < NUMBER_OF_PARTITION ; i++) {
		free(B[i]); free(out[i]);
	} free(X); free(B); free(out);

	return 0;
}
