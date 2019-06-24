#include "simplex.h"

float absolute(float a) {
	return a < 0 ? -a : a;
}

// Returns the index of the minimum value of the column
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

// Multiply the line by a factor
void multiply_line_by_factor(float ** table, int nb_column, int line, float factor) {
	int j;
	for (j = 0 ; j < nb_column ; j++) {
		table[line][j] *= factor;
	}
}

// Add a line by a line
void add_line_by_line(float ** table, int nb_column, int line0, int line1) {
	int j;
	for (j = 0 ; j < nb_column ; j++) {
		table[line0][j] += table[line1][j];
	}
}

// Substract a line by a line multiplied by a factor
void sub_line_by_line_by_factor(float ** table, int nb_column, int line0, int line1, float factor) {
	int j;
	for (j = 0 ; j < nb_column ; j++) {
		table[line0][j] -= table[line1][j] * factor;
	}
}

// Natural pivot for simplex method
void pivot(float ** table, int nb_line, int nb_column, int * line, int * column, int species) {
	*line = -1;
	*column = -1;
	int i, j;
	float max = ERROR, min = 1e38;
	// Find the best column
	for (j = 0 ; j < nb_column+species-3 ; j++) {
		if (table[nb_line-species][j] > max) {
			max = table[nb_line-species][j];
			*column = j;
		}
	}
	if (*column == -1) { return; }
	// Find the best line
	for (i = 0 ; i < nb_line-2 ; i++) {
		float ratio = table[i][nb_column-1]/table[i][*column];
		if (table[i][*column] > ERROR && ratio < min) {
			min = ratio;
			*line = i;
		}
	}
}

void pivot_bland(float ** table, int nb_line, int nb_column, int * line, int * column, int species) {
	*line = -1;
	*column = -1;
	int i, j;
	float max = ERROR, min = 1e38;
	// Find the best column
	for (j = 0 ; j < nb_column+species-3 ; j++) {
		if (table[nb_line-species][j] > max) {
			*column = j;
			break;
		}
	}
	if (*column == -1) { return; }
	// Find the best line
	for (i = 0 ; i < nb_line-2 ; i++) {
		float ratio = table[i][nb_column-1]/table[i][*column];
		if (table[i][*column] > ERROR && ratio < min) {
			min = ratio;
			*line = i;
		}
	}
}

// Find in the basis which index contains the column
int find_basis_line(int * basis, int len, int column) {
	int i;
	for (i = 0 ; i < len ; i++) {
		if (basis[i] == column) {
			return i;
		}
	}
	return -1;
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

// Initialize the table for the problem
float ** initialize_table(int n, int nb_line, int nb_column, int nb_var, float * X, int ** B) {
	printf("start of the procedure\n");
	clock_t begin = clock();

	int i, j, k = 0;
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
					table[i*2][k] = 1.0;
					table[i*2+1][k] = -1.0;
				}
				if (j < n-1) {
					table[j*2][k] = -1.0;
					table[j*2+1][k] = 1.0;
				}
				k++;
			}
		}
	}

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time0 spent is : %.9lf\n", time_spent);

	return table;
}

int * initialize_basis(int nb_line) {
	int i;
	int * basis = (int *) calloc(nb_line-2, sizeof(int));
	for (i = 0 ; i < nb_line-2 ; i++) {
		basis[i] = i + nb_line-2;
	}
	return basis;
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
	
	#ifdef PRINT
	printf("nb_var = %d\n", nb_var);
	#endif

	int nb_line = 2 * n;
	int nb_column = nb_var + nb_line;
	float ** table = initialize_table(n, nb_line, nb_column, nb_var, X, B);

	int * basis = initialize_basis(nb_line);

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

	printf("start of the procedure\n");
	clock_t begin = clock();

	// Find the index which contains the minimum value in the last column 
	k = index_min_by_column(table, nb_line, nb_column-1);

	// column k leave the basis, column line_column-2 enter in the basis
	basis[k] = nb_column-2;

	// Eliminate the negative values ​​from the last column
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

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time1 spent is : %.9lf\n", time_spent);
	
	// STEP 2 : SIMPLEX
	#ifdef PRINT
	printf("\n-------------------\n"
			 "----- SIMPLEX -----\n"
			 "-------------------\n");
	#endif

	printf("start of the procedure\n");
	begin = clock();

	while (1) {
		pivot_bland(table, nb_line, nb_column, &i, &j, 2);
		if (i == -1) { break; }

		#ifdef PRINT
		printf("pivot is : line %d column %d\n", i, j);
		#endif

		// column j enter in the basis
		basis[i] = j;

		// Normalize the line, then reduce column
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

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time2 spent is : %.9lf\n", time_spent);

	// STEP 3 : INITIALIZATION
	#ifdef PRINT
	printf("\n--------------------------\n"
			 "----- INITIALIZATION -----\n"
			 "--------------------------\n");
	#endif

	printf("start of the procedure\n");
	begin = clock();

	// Determine if the column nb_column-2 is in the basis
	i = find_basis_line(basis, nb_line-2, nb_column-2);
	// If it's in the basis, you have to extract it 
	if (i != -1) {
		// Searching for the first correct pivot
		for (j = 0 ; j < nb_column-2 ; j++) {
			if (absolute(table[i][j]) > ERROR) { break; }
		}

		if (j == nb_column-1) { printf("Error\n"); exit(1); }

		#ifdef PRINT
		printf("pivot is : line %d column %d\n", i, j);
		#endif

		// Column j enter in the basis
		basis[i] = j;

		// Normalize the line, then reduce column
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

	// Reduce the table, transition to a problem of the first species
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

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time3 spent is : %.9lf\n", time_spent);

	// STEP 4 : SIMPLEX
	#ifdef PRINT
	printf("\n-------------------\n"
			 "----- SIMPLEX -----\n"
			 "-------------------\n");
	#endif

	printf("start of the procedure\n");
	begin = clock();

	while (1) {
		pivot_bland(table, nb_line, nb_column, &i, &j, 1);
		if (i == -1) { break; }

		#ifdef PRINT
		printf("pivot is : line %d column %d\n", i, j);
		#endif

		// Column j enter in the basis
		basis[i] = j;

		// Normalize the line, then reduce column
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

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time4 spent is : %.9lf\n", time_spent);

	// STEP 5 : SOLUTION
	#ifdef PRINT
	printf("\n--------------------\n"
			 "----- SOLUTION -----\n"
			 "--------------------\n");
	#endif

	printf("start of the procedure\n");
	begin = clock();

	// Get the solution of the problem
	float * x = (float *) calloc(nb_var, sizeof(float));
	for (j = 0 ; j < nb_var ; j++) {
		i = find_basis_line(basis, nb_line-2, j);
		if (i != -1) {
			x[j] = table[i][nb_column-1];
		}
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

	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time5 spent is : %.9lf\n", time_spent);
	return alpha;
}