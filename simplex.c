#include "simplex.h"

float absolute(float a) {
	return a < 0 ? -a : a;
}

void print_simplex(Simplex * ptr_simplex) {
	int i, j;
	printf("\nnb_var = %d nb_row = %d nb_column = %d\n", ptr_simplex->nb_var, ptr_simplex->nb_row, ptr_simplex->nb_column);
	for (i = 0 ; i < ptr_simplex->nb_row ; i++) {
		for (j = 0 ; j < ptr_simplex->nb_column ; j++) {
			if (ptr_simplex->A[i * ptr_simplex->nb_column + j] < 0)
				printf("%.1f ", ptr_simplex->A[i * ptr_simplex->nb_column + j]);
			else if (ptr_simplex->A[i * ptr_simplex->nb_column + j] > 0)
				printf("%.2f ", ptr_simplex->A[i * ptr_simplex->nb_column + j]);
			else
				printf("0.00 ");
		}
		if (ptr_simplex->b[i] < 0)
			printf("| %.1f\n", ptr_simplex->b[i]);
		else if (ptr_simplex->b[i] > 0)
			printf("| %.2f\n", ptr_simplex->b[i]);
		else
			printf("| 0.00\n");
	}
	printf("-----\n");
	for (i = 0 ; i < ptr_simplex->nb_column ; i++) {
		if (ptr_simplex->c[i] < 0)
			printf("%.1f ", ptr_simplex->c[i]);
		else if (ptr_simplex->c[i] > 0)
			printf("%.2f ", ptr_simplex->c[i]);
		else
			printf("0.00 ");
	}
	if (ptr_simplex->obj < 0)
		printf("| %.1f\n", ptr_simplex->obj);
	else if (ptr_simplex->obj > 0)
		printf("| %.2f\n", ptr_simplex->obj);
	else
		printf("| 0.00\n");

	printf("\nbasis = [ ");
	for (i = 0 ; i < ptr_simplex->nb_row ; i++) {
		printf("%d ", ptr_simplex->basis[i]);
	} printf("]\n");
}

Simplex * init_simplex(float * X, int ** B, int n) {
	int i, j, k = 0;
	Simplex * ptr_simplex = (Simplex *) malloc(sizeof(Simplex));
	
	// Count the number of variables
	for (i = 0 ; i < n ; i ++)
		for (j = 0 ; j < n ; j++)
			ptr_simplex->nb_var += B[i][j];
	ptr_simplex->nb_row = n-1;
	ptr_simplex->nb_column = ptr_simplex->nb_var + ptr_simplex->nb_row;
	
	ptr_simplex->A = (float *) calloc(ptr_simplex->nb_row * ptr_simplex->nb_column, sizeof(float));
	// Main Block
	for (i = 0 ; i < n ; i++) {
		for (j = 0 ; j < n ; j++) {
			if (B[i][j]) {
				if (i < n-1)
					ptr_simplex->A[i * ptr_simplex->nb_column + k] = 1.0;
				if (j < n-1)
					ptr_simplex->A[j * ptr_simplex->nb_column + k] = -1.0;
				k++;
			}
		}
	}
	// Eye
	for (i = 0 ; i < ptr_simplex->nb_row ; i++)
		ptr_simplex->A[i * ptr_simplex->nb_column + ptr_simplex->nb_var + i] = 1.0;
	
	// Bounds
	ptr_simplex->b = (float *) malloc(ptr_simplex->nb_row * sizeof(float));
	for (i = 0 ; i < ptr_simplex->nb_row ; i++) {
		if (X[i] < 0) {
			for (j = 0 ; j < ptr_simplex->nb_var ; j++)
				ptr_simplex->A[i * ptr_simplex->nb_column + j] *= -1.0;
			ptr_simplex->b[i] = -X[i];
		}
		else
			ptr_simplex->b[i] = X[i];
	}

	// Prices
	ptr_simplex->c = (float *) calloc(ptr_simplex->nb_var * ptr_simplex->nb_row, sizeof(float));
	ptr_simplex->obj = 0;
	for (i = 0 ; i < ptr_simplex->nb_row ; i++) {
		for (j = 0 ; j < ptr_simplex->nb_var ; j++)
			ptr_simplex->c[j] += ptr_simplex->A[i * ptr_simplex->nb_column + j];
		ptr_simplex->obj += ptr_simplex->b[i];
	}

	// Basis
	ptr_simplex->basis = (int *) calloc(ptr_simplex->nb_row, sizeof(int));
	for (i = 0 ; i < ptr_simplex->nb_row ; i++)
		ptr_simplex->basis[i] = i + ptr_simplex->nb_var;

	return ptr_simplex;
}

// Multiply the row by a factor
void multiply_row_by_factor(Simplex * ptr_simplex, int row, float factor) {
	int j;
	for (j = 0 ; j < ptr_simplex->nb_column ; j++)
		ptr_simplex->A[row * ptr_simplex->nb_column + j] *= factor;
	ptr_simplex->b[row] *= factor;
}

// Substract a row by a row multiplied by a factor
void sub_row_by_row_by_factor(Simplex * ptr_simplex, int row0, int row1, float factor) {
	int j;
	for (j = 0 ; j < ptr_simplex->nb_column ; j++)
		ptr_simplex->A[row0 * ptr_simplex->nb_column + j] -= ptr_simplex->A[row1 * ptr_simplex->nb_column + j] * factor;
	ptr_simplex->b[row0] -= ptr_simplex->b[row1] * factor;
}

// Natural pivot for simplex method
void pivot(Simplex * ptr_simplex, int * row, int * column) {
	*row = -1;
	*column = -1;
	int i, j;
	float max = ERROR, min = 1e38;
	// Find the best column
	for (j = 0 ; j < ptr_simplex->nb_column ; j++) {
		float cur = ptr_simplex->c[j];
		if (cur > max) {
			max = cur;
			*column = j;
		}
	}
	if (*column == -1) { return; }
	// Find the best row
	for (i = 0 ; i < ptr_simplex->nb_row ; i++) {
		float cur = ptr_simplex->A[i * ptr_simplex->nb_column + *column];
		float ratio = ptr_simplex->b[i]/cur;
		if (cur > ERROR && ratio < min) {
			min = ratio;
			*row = i;
		}
	}
}

void simplex_method(Simplex * ptr_simplex) {
	int i, j, k;
	while (1) {
		pivot(ptr_simplex, &i, &j );
		if (i == -1) { break; }

		printf("\npivot is : row %d column %d\n", i, j);

		// Column j enter in the basis
		ptr_simplex->basis[i] = j;

		// Normalize the row, then reduce column
		float factor = 1/ptr_simplex->A[i * ptr_simplex->nb_column + j];
		multiply_row_by_factor(ptr_simplex, i, factor);
		for (k = 0 ; k < ptr_simplex->nb_row ; k++) {
			if (k != i) {
				float factor = ptr_simplex->A[k * ptr_simplex->nb_column + j];
				sub_row_by_row_by_factor(ptr_simplex, k, i, factor);
			}
		}
		factor = ptr_simplex->c[j];
		for (k = 0 ; k < ptr_simplex->nb_column ; k++)
			ptr_simplex->c[k] -= ptr_simplex->A[i * ptr_simplex->nb_column + k] * factor;
		ptr_simplex->obj -= ptr_simplex->b[i] * factor;

		print_simplex(ptr_simplex);
	}
}

void simplex_procedure(float * X, int ** B, int n) {
	Simplex * ptr_simplex = init_simplex(X, B, n);
	print_simplex(ptr_simplex);
	simplex_method(ptr_simplex);
}