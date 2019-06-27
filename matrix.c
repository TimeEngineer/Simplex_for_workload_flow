#include "matrix.h"

void print_vector(char * name, float * vector, int len) {
	int i;
	printf("%s = [ ", name);
	for (i = 0 ; i < len ; i++) {
		if (vector[i] == 0.0) {
			printf("0.00 ");
		}
		else if (vector[i] > 0.0) {
			printf("%.2f ", vector[i]);
		}
		else {
			printf("%.1f ", vector[i]);	
		}
	} printf("]\n");
}

void print_matrix(char * name, float ** matrix, int nb_line, int nb_column) {
	int i, j;
	printf("%s =\n", name);
	for (i = 0 ; i < nb_line ; i++) {
		printf("[ ");
		for (j = 0 ; j < nb_column ; j++) {
			if (matrix[i][j] == 0.0) {
				printf("0.00 ");
			}
			else if (matrix[i][j] > 0) {
				printf("%.2f ", matrix[i][j]);
			}
			else {
				printf("%.1f ", matrix[i][j]);	
			}
		} printf("]\n");
	}
}