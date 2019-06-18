#include "simplex.h"

float absolute(float a) {
	return a < 0 ? -a : a;
}

int index_min_by_column(float ** table, int nb_line, int column) {
	int index = 0, i = 0;
	float min = table[i][column];
	for (i = 1 ; i < nb_line ; i++) {
		if (table[i][column] < min) {
			min = table[i][column];
			index = i;
		}
	}
	return index;
}

void multiply_line_by_factor(float ** table, int nb_column, int line, float factor) {
	int j;
	for (j = 0 ; j < nb_column ; j++) {
		table[line][j] *= factor;
	}
}

void add_line_by_line(float ** table, int nb_column, int line0, int line1) {
	int j;
	for (j = 0 ; j < nb_column ; j++) {
		table[line0][j] += table[line1][j];
	}
}

void sub_line_by_line_by_factor(float ** table, int nb_column, int line0, int line1, float factor) {
	int j;
	for (j = 0 ; j < nb_column ; j++) {
		table[line0][j] -= table[line1][j] * factor;
	}
}

void pivot(float ** table, int nb_line, int nb_column, int * line, int * column, int species) {
	*line = -1;
	*column = -1;
	int i, j;
	float max = ERROR, min = 1e38;
	for (j = 0 ; j < nb_column+species-3 ; j++) {
		if (table[nb_line-species][j] > max) {
			max = table[nb_line-species][j];
			*column = j;
		}
	}
	if (*column == -1) { return; }
	for (i = 0 ; i < nb_line-2 ; i++) {
		float ratio = table[i][nb_column-1]/table[i][*column];
		if (table[i][*column] > ERROR && ratio < min) {
			min = ratio;
			*line = i;
		}
	}
}

int find_basis_variable(float ** table, int nb_line, int column) {
	int i, index = -1;
	for(i = 0 ; i < nb_line-2 ; i++) {
		if (absolute(table[i][column] - 1.0) < ERROR) {
			if (index == -1) {
				index = i;   // found first '1', save this row number.
			}
			else {
				return -1; // found second '1'.
			}
		}
		else if (absolute(table[i][column]) > ERROR) {
			return -1;
		}
	}
	return index;
}

#ifdef PRINT
void print_table(float ** table, int nb_line, int nb_column, int species, int pivot_line, int pivot_column) {
	int i, j;
	printf("\nThe table is\n");
	for (i = 0 ; i < nb_line ; i++) {
		if (species == 1 && i == nb_line-2) { continue; }
		if (i > nb_line - 3) {
			printf("------\n");
		}
		for (j = 0 ; j < nb_column ; j++) {
			if (i == pivot_line && j == pivot_column) { printf("\033[0;31m"); }
			if (species == 1 && j == nb_column-2) { continue; }
			if (j > nb_column - 3) {
				printf("| ");			
			}
			if (table[i][j] == 0.0) {
				printf("0.00000000 ");
			}
			else if (table[i][j] > 0.0) {
				printf("%.8f ", table[i][j]);
			}
			else {
				printf("%.7f ", table[i][j]);	
			}
			if (i == pivot_line && j == pivot_column) { printf("\033[0m"); }
		} printf("\n");
	}
	printf("\n");
}
#endif

void print_vector(char * name, float * vector, int len) {
	int i;
	printf("%s = [ ", name);
	for (i = 0 ; i < len ; i++) {
		if (vector[i] == 0.0) {
			printf("0.00000000 ");
		}
		else if (vector[i] > 0.0) {
			printf("%.8f ", vector[i]);
		}
		else {
			printf("%.7f ", vector[i]);	
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
				printf("0.00000000 ");
			}
			else if (matrix[i][j] > 0) {
				printf("%.8f ", matrix[i][j]);
			}
			else {
				printf("%.7f ", matrix[i][j]);	
			}
		} printf("]\n");
	}
}

void print_basis(char * name, int * vector, int len) {
	int i;
	printf("%s = [ ", name);
	for (i = 0 ; i < len ; i++) {
		printf("%d ", vector[i]);
	} printf("]\n");
}

// search in the basis which contains the line number
int search_basis_line(int * basis, int len, int line) {
	int i;
	for (i = 0 ; i < len ; i++) {
		if (basis[i] == line) {
			return i;
		}
	}
}

float ** simplex_procedure(float * X, int ** B, int n) {
	int i, j, k = 0, nb_var = 0;
	// Initialization of the table

	for (i = 0 ; i < n ; i ++) {
		for (j = 0 ; j < n ; j++) {
			// count the number of variables
			nb_var += B[i][j];
		}
	}
	// nb_var -= n; // if we put 1 for B[i][i]
	printf("nb_var = %d\n", nb_var);

	int nb_line = 2 * n;
	int nb_column = nb_var + nb_line;
	float ** table = (float **) malloc(nb_line * sizeof(float *));
	for (i = 0 ; i < nb_line ; i ++) {
		table[i] = (float *) calloc(nb_column, sizeof(float));
		for (j = 0 ; j < nb_column ; j++) {
			if (j >= nb_var && j < nb_column - 2 && i < nb_line - 2) { // Identity
				if (i == j - nb_var) {
					table[i][j] = 1.0;
				}
			}
			else if (j == nb_column - 1) { // Last column
				if (i < nb_line - 2) {
					if (i % 2) { // odd
						table[i][j] = -X[i/2];
					}
					else { // even
						table[i][j] = X[i/2];
					}
				}
			}
			else if (j < nb_var && i == nb_line - 1) { // Last line
				table[i][j] = -1.0;
			}
			else if (j == nb_column - 2) { // Temporary column
				table[i][j] = -1.0;
			}
		}
	}
	for (i = 0 ; i < n ; i++) { // Main block
		for (j = 0 ; j < n ; j++) {
			if (B[i][j]) {
				if (i < n-1) {
					table[i*2][k] = X[i];
					table[i*2+1][k] = -X[i];
				}
				if (j < n-1) {
					table[j*2][k] = -X[i];
					table[j*2+1][k] = X[i];
				}
				k++;
			}
		}
	}
	int * basis = (int *) calloc(nb_line-2, sizeof(int));
	for (i = 0 ; i < nb_line-2 ; i++) {
		basis[i] = i + nb_line-2;
	}
	#ifdef PRINT
	print_table(table, nb_line, nb_column, 2, -1, -1);
	print_basis("basis", basis, nb_line-2);
	#endif
	// Begin of the procedure

	// STEP 1 : PRE-INITIALIZATION :
	#ifdef PRINT
	printf("\n------------------------------\n"
			 "----- PRE-INITIALIZATION -----\n"
			 "------------------------------\n");
	#endif
	k = index_min_by_column(table, nb_line, nb_column-1);
	// column k out, column line_column-2 in
	basis[k] = nb_column-2;
	multiply_line_by_factor(table, nb_column, k, -1);
	for (i = 0 ; i < nb_line - 1 ; i++) {
		if (i != k) {
			add_line_by_line(table, nb_column, i, k);
		}
	}
	#ifdef PRINT
	print_table(table, nb_line, nb_column, 2, -1, -1);
	print_basis("basis", basis, nb_line-2);
	#endif
	// STEP 2 : SIMPLEX
	#ifdef PRINT
	printf("\n-------------------\n"
			 "----- SIMPLEX -----\n"
			 "-------------------\n");
	#endif
	while (1) {
		pivot(table, nb_line, nb_column, &i, &j, 2);
		if (i == -1) {
			break;
		}
		#ifdef PRINT
		printf("pivot is : line %d column %d\n", i, j);
		#endif
		// column j in
		basis[i] = j;
		multiply_line_by_factor(table, nb_column, i, 1/table[i][j]);

		for (k = 0 ; k < nb_line-1 ; k++) {
			if (k != i) {
				sub_line_by_line_by_factor(table, nb_column, k, i, table[k][j]);
			}
		}
		#ifdef PRINT
		print_table(table, nb_line, nb_column, 2, i, j);
		print_basis("basis", basis, nb_line-2);
		#endif
	}

	// STEP 3 : INITIALIZATION
	#ifdef PRINT
	printf("\n--------------------------\n"
			 "----- INITIALIZATION -----\n"
			 "--------------------------\n");
	#endif
	i = find_basis_variable(table, nb_line, nb_column-2);
	if (i != -1) {
		// Searching for the first correct pivot
		for (j = 0 ; j < nb_column-2 ; j++) {
			if (absolute(table[i][j]) > ERROR) {
				break;
			}
		}
		#ifdef PRINT
		printf("pivot is : line %d column %d\n", i, j);
		#endif
		// column j in
		basis[i] = j;
		multiply_line_by_factor(table, nb_column, i, 1/table[i][j]);

		for (k = 0 ; k < nb_line-1 ; k++) {
			if (k != i) {
				sub_line_by_line_by_factor(table, nb_column, k, i, table[k][j]);
			}
		}
		#ifdef PRINT
		print_table(table, nb_line, nb_column, 2, i, j);
		print_basis("basis", basis, nb_line-2);
		#endif
	}
	for (i = 0 ; i < nb_line-2 ; i++) {
		j = basis[i];
		for (k = 0 ; k < nb_line ; k++) {
			if (k != i && k != nb_line-2) {
				sub_line_by_line_by_factor(table, nb_column, k, i, table[k][j]);
			}
		}
	}
	#ifdef PRINT
	print_table(table, nb_line, nb_column, 1, -1, -1);
	#endif
	// STEP 4 : SIMPLEX
	#ifdef PRINT
	printf("\n-------------------\n"
			 "----- SIMPLEX -----\n"
			 "-------------------\n");
	#endif
	while(1) {
		pivot(table, nb_line, nb_column, &i, &j, 1);
		if (i == -1) { break; }
		#ifdef PRINT
		printf("pivot is : line %d column %d\n", i, j);
		#endif
		// column j in
		basis[i] = j;
		multiply_line_by_factor(table, nb_column, i, 1/table[i][j]);

		for (k = 0 ; k < nb_line ; k++) {
			if (k != i && k != nb_line-2) {
				sub_line_by_line_by_factor(table, nb_column, k, i, table[k][j]);
			}
		}
		#ifdef PRINT
		print_table(table, nb_line, nb_column, 1, i, j);
		print_basis("basis", basis, nb_line-2);
		#endif
	}

	// STEP 5 : SOLUTION
	#ifdef PRINT
	printf("\n--------------------\n"
			 "----- SOLUTION -----\n"
			 "--------------------\n");
	#endif
	float * x = (float *) calloc(nb_var, sizeof(float));
	for (j = 0 ; j < nb_var ; j++) {
		i = search_basis_line(basis, nb_line-2, j);
		x[j] = table[i][nb_column-1];
	}
	#ifdef PRINT
	printf("The solution is :\n");
	print_vector("x", x, nb_var);
	#endif

	// Free table
	free(basis);
	for (i = 0 ; i < nb_line ; i++) {
		free(table[i]);
	} free(table);

	// Initialization of the output
	float ** alpha = (float **) malloc(n * sizeof(float *));
	k = 0;
	for (i = 0 ; i < n ; i++) {
		alpha[i] = (float *) calloc(n, sizeof(float));
		for (j = 0 ; j < n ; j++) {
			if (B[i][j]) {
				alpha[i][j] = x[k];
				k++;
			}
		}
	}
	free(x);
	return alpha;
}