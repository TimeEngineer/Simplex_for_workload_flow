#include "simplex.h"
#include <time.h>
// #define PRINT
#define N 1024

float average(float * vector, int len) {
	int i;
	float sum = 0;
	for (i = 0 ; i < len ; i++) {
		sum += vector[i];
	}
	return sum/(float)len;
}

int ** init_B(int n) {
	int i, j;
	int ** B = (int **) malloc(n * sizeof(int *));
	for (i = 0 ; i < n ; i++) {
		B[i] = (int *) calloc(n, sizeof(int));
	}
	for (i = 0 ; i < n-1 ; i++) {
		for (j = i+1 ; j < n ; j++) {
			if (j < i+11) {
				B[i][j] = 1;
				B[j][i] = 1;
			}
		}
	}
	return B;
}

float * init_X(int n) {
	int i;
	float * X = (float *) malloc(n * sizeof(float));
	for (i = 0 ; i < n ; i++) {
		X[i] = ((float)rand())/((float)RAND_MAX)*100.0;
		// X[i] = (double) i + 1;
	}
	float avg = average(X, n);
	for (i = 0 ; i < n ; i++) {
		X[i] -= avg;
	}
	return X;
}

int main() {
	time_t t;
	srand(time(&t));
	int i, j;
	// Initialization
	float * X = init_X(N);
	int ** B = init_B(N);

	// for (i = 0 ; i < N ; i++) {
	//  	for (j = 0 ; j < N ; j++) {
	//  		printf("%d", B[i][j]);
	//  	} printf("\n");
	// } printf("\n");

	// print_vector("X", X, N);

/* -------------------------------------------------- */

	float ** alpha = simplex_procedure(X, B, N);
	
/* -------------------------------------------------- */

	// print_matrix("a", alpha, N, N);

	free(X);
	for (i = 0 ; i < N ; i++) {
		free(B[i]);
	} free(B);
	for (i = 0 ; i < N ; i++) {
		free(alpha[i]);
	} free(alpha);
	
	return 0;
}
