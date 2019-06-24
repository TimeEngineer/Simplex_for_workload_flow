#include <stdio.h>            /* C input/output                       */
#include <stdlib.h>           /* C standard library                   */
#include <glpk.h>             /* GNU GLPK linear/mixed integer solver */

#define N 256

void print_matrix(char * name, double ** matrix, int nb_line, int nb_column) {
	int i, j;
	printf("%s =\n", name);
	for (i = 0 ; i < nb_line ; i++) {
		printf("[ ");
		for (j = 0 ; j < nb_column ; j++) {
			if (matrix[i][j] == 0.0) {
				printf("0.0000 ");
			}
			else if (matrix[i][j] > 0) {
				printf("%.4lf ", matrix[i][j]);
			}
			else {
				printf("%.3lf ", matrix[i][j]);	
			}
		} printf("]\n");
	}
}

void print_vector(char * name, double * vector, int len) {
	int i;
	printf("%s = [ ", name);
	for (i = 0 ; i < len ; i++) {
		if (vector[i] == 0.0) {
			printf("0.0000 ");
		}
		else if (vector[i] > 0.0) {
			printf("%.4lf ", vector[i]);
		}
		else {
			printf("%.3lf ", vector[i]);	
		}
	} printf("]\n");
}

int ** init_B(int n, int * nb_var, int * count) {
	int i, j;
	*nb_var = 0, *count = 0;
	int ** B = (int **) malloc(n * sizeof(int *));
	for (i = 0 ; i < n ; i++) {
		B[i] = (int *) malloc(n * sizeof(int));
		for (j = 0 ; j < n ; j++) {
			if ((i+j) % 2) {
				B[i][j] = 1;
				if (i < n-1) {
					*count += 1;
				} else {
					*nb_var += 1;
				}
			} else {
				B[i][j] = 0;
			}
		}
	}
	*nb_var += *count;
	*count <<= 2;
	return B;
}

double ** init_table(int ** B, int n, int nb_line, int nb_var) {
	int i, j, k = 0;
	double ** table = (double **) malloc(nb_line * sizeof(double *));
	for (i = 0 ; i < nb_line ; i ++) {
		table[i] = (double *) calloc(nb_var, sizeof(double));
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
	return table;
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
		// X[i] = ((float)rand())/(float)(RAND_MAX) * 100.0;
		X[i] = (double) i + 1;
	}
	double avg = average(X, n);
	for (i = 0 ; i < n ; i++) {
		X[i] -= avg;
	}
	return X;
}

#include <time.h>

int main(void) {
	int i, j, k = 0;
	int nb_var, count;
	int nb_line = 2*(N-1);

	int ** B = init_B(N, &nb_var, &count);
	double * X = init_X(N);
	printf("nb_var = %d, count = %d\n", nb_var, count);

	printf("start of the procedure\n");
	clock_t begin = clock();

	double ** table = init_table(B, N, nb_line, nb_var);

	// print_matrix("table", table, nb_line, nb_var);
	// print_vector("X", X, N);

	/* declare variables */
	glp_prob *lp;
	int ia[1+131000], ja[1+131000];
	double ar[1+131000], z;
	/* create problem */
	lp = glp_create_prob();
	glp_set_prob_name(lp, "partition optimizer");
	glp_set_obj_dir(lp, GLP_MAX);
	/* fill problem */
	glp_add_rows(lp, 2*(N-1));

	char name[8];
	for (i = 1 ; i <= 2*(N-1) ; i++) {
		sprintf(name, "%d", i);
		glp_set_row_name(lp, i, name);
		if (i%2) {
			glp_set_row_bnds(lp, i, GLP_UP, 0.0, X[(i-1)/2]);
		} else {
			glp_set_row_bnds(lp, i, GLP_UP, 0.0, -X[(i-1)/2]);
		}
	}
	
	glp_add_cols(lp, nb_var);
	for (i = 1 ; i <= nb_var ; i++) {
		sprintf(name, "x%d", i);
		glp_set_col_name(lp, i, name);
		glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0);
		glp_set_obj_coef(lp, i, -1.0);
	}
	j = 0;
	for (i = 1 ; i <= count ; i++) {
		while (table[j/nb_var][j%nb_var] == 0) {
			j++;
			if (j == nb_line * nb_var) { break; } 
		}
		if (j == nb_line * nb_var) { break; }
		ia[i] = j/nb_var+1, ja[i] = j%nb_var+1, ar[i] = table[j/nb_var][j%nb_var];
		j++;
	}

	for (i = 0 ; i < nb_line ; i++) {
		free(table[i]);
	} free(table);

	glp_load_matrix(lp, count, ia, ja, ar);
	/* solve problem */
	glp_simplex(lp, NULL);
	/* recover and display results */
	z = glp_get_obj_val(lp);
	double * x = (double *) calloc(nb_var, sizeof(double));
	for (i = 0 ; i < nb_var ; i++) {
		x[i] = glp_get_col_prim(lp, i+1);
	}
	// printf("z = %g\n", z);
	// print_vector("x", x, nb_var);

	double ** alpha = (double **) malloc(N * sizeof(double *));
	k = 0;
	for (i = 0 ; i < N ; i++) {
		alpha[i] = (double *) calloc(N, sizeof(double));
		for (j = 0 ; j < N ; j++) {
			if (B[i][j]) {
				alpha[i][j] = x[k];
				k++;
			}
		}
	}
	
	for (i = 0 ; i < N ; i++) {
		free(alpha[i]);
		free(B[i]);
	} free(alpha); free(B); free(X);
	// print_matrix("a", alpha, N, N);
	
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time spent is : %.9lf\n", time_spent);
	/* housekeeping */
	glp_delete_prob(lp);
	glp_free_env();
	return 0;
}