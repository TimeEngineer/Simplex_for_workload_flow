#include "simplex.h"
#include <time.h>
// #define PRINT
#define MEASURE
#define N 1024

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

	#ifdef PRINT
	printf("B =\n");
	#endif

	for (i = 0 ; i < N ; i++) {
		#ifdef PRINT
		printf("[ ");
		#endif
		// X[i] = ((float)rand())/(float)(RAND_MAX) * 100.0;
		X[i] = (float) i + 1;
		B[i] = (int *) malloc(N * sizeof(int));
		for (j = 0 ; j < N ; j++) {
			if ((i+j) % 256 == 255) {
				B[i][j] = 1;
				#ifdef PRINT
				printf("1 ");
				#endif
			} else {
				B[i][j] = 0;
				#ifdef PRINT
				printf("0 ");
				#endif
			}
		}
		#ifdef PRINT
		printf("]\n");
		#endif
	}

	#ifdef PRINT
	print_vector("N", X, N);
	#endif

	float avg = average(X, N);
	for (i = 0 ; i < N ; i++) {
		X[i] -= avg;
	}
	#ifdef PRINT
	print_vector("X", X, N);
	#endif

/* -------------------------------------------------- */
	#ifdef MEASURE
	printf("start of the procedure\n");
	clock_t begin = clock();
	#endif

	float ** alpha = simplex_procedure(X, B, N);
	
	#ifdef MEASURE
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time spent is : %.9lf\n", time_spent);
	#endif
/* -------------------------------------------------- */

	#ifdef PRINT
	printf("\nThe solution is :\n");
	print_matrix("a", alpha, N, N);
	#endif

	free(X);
	for (i = 0 ; i < N ; i++) {
		free(B[i]);
	} free(B);
	for (i = 0 ; i < N ; i++) {
		free(alpha[i]);
	} free(alpha);
	
	return 0;
}
