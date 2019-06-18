#include "simplex.h"
#include <time.h>
#define N 60

float average(float * vector, int len) {
	int i;
	float sum = 0;
	for (i = 0 ; i < len ; i++) {
		sum += vector[i];
	}
	return sum/(float)len;
}

int main() {
	time_t t;
	srand(time(&t));
	int i, j;
	// Initialization
	float * X = (float *) malloc(N * sizeof(float));
	int ** B = (int **) malloc(N * sizeof(int *));
	printf("B =\n");
	for (i = 0 ; i < N ; i++) {
		printf("[ ");
		X[i] = ((float)rand())/(float)(RAND_MAX) * 100.0;
		// X[i] = (float) i + 1;
		B[i] = (int *) malloc(N * sizeof(int));
		for (j = 0 ; j < N ; j++) {
			if ((i+j) % 2) {
				B[i][j] = 1;
				printf("1 ");
			} else {
				B[i][j] = 0;
				printf("0 ");
			}
		}
		printf("]\n");
	}
	print_vector("N", X, N);
	float avg = average(X, N);
	for (i = 0 ; i < N ; i++) {
		X[i] -= avg;
	}
	print_vector("X", X, N);
	float ** alpha = simplex_procedure(X, B, N);

	printf("\nThe solution is :\n");
	print_matrix("a", alpha, N, N);

	free(X);
	for (i = 0 ; i < N ; i++) {
		free(B[i]);
	} free(B);
	for (i = 0 ; i < N ; i++) {
		free(alpha[i]);
	} free(alpha);
	
	return 0;
}