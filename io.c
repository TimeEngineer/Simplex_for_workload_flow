#include "io.h"

void print_matrix(char * name, double ** matrix, int nb_row, int nb_column) {
	int i, j;
	printf("%s =\n", name);
	for (i = 0 ; i < nb_row ; i++) {
		printf("[ ");
		for (j = 0 ; j < nb_column ; j++) {
			if (matrix[i][j] == 0.0) {
				printf("0.00 ");
			}
			else if (matrix[i][j] > 0) {
				printf("%.2lf ", matrix[i][j]);
			}
			else {
				printf("%.1lf ", matrix[i][j]);	
			}
		} printf("]\n");
	}
}

void print_vector(char * name, double * vector, int len) {
	int i;
	printf("%s = [ ", name);
	for (i = 0 ; i < len ; i++) {
		if (vector[i] == 0.0) {
			printf("0.00 ");
		}
		else if (vector[i] > 0.0) {
			printf("%.2lf ", vector[i]);
		}
		else {
			printf("%.1lf ", vector[i]);	
		}
	} printf("]\n");
}

void print_B(int ** B, int len) {
	int i, j;
	printf("B =\n");
	for (i = 0 ; i < len ; i++) {
	 	for (j = 0 ; j < len ; j++) {
	 		printf("%d", B[i][j]);
	 	} printf("\n");
	} printf("\n");
}