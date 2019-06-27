#include "simplex.h"

double ** init_table(int ** B, int n, int nb_row, int nb_var) {
	int i, j, k = 0;
	double ** table = (double **) malloc(nb_row * sizeof(double *));
	for (i = 0 ; i < nb_row ; i ++) {
		table[i] = (double *) calloc(nb_var, sizeof(double));
	}
	for (i = 0 ; i < n ; i++) { // Main block
		for (j = 0 ; j < n ; j++) {
			if (B[i][j]) {
				if (i < n-1) {
					table[i][k] = 1.0;
				}
				if (j < n-1) {
					table[j][k] = -1.0;
				}
				k++;
			}
		}
	}
	return table;
}

double ** simplex_procedure(double * X, int ** B, int n) {
	int i, j, k = 0;

	/* count the number of variables and the number of non-zeros */
	int nb_var = 0, nb_no0 = 0;
	for (i = 0 ; i < n-1 ; i++) {
		for (j = i+1 ; j < n ; j++) {
			if (B[i][j]) {
				if (j < n-1) {
					nb_no0 += 2;
				} else {
					nb_no0++;
					nb_var++;
				}
			}
		}
	}
	nb_var += nb_no0;
	nb_no0 <<= 1;

	/* number of row */
	int nb_row = n-1;
	
	/* initialize the table */
	double ** table = init_table(B, n, nb_row, nb_var);

	/* declare variables */
	glp_prob *lp;
	int ia[1+nb_no0], ja[1+nb_no0];
	double ar[1+nb_no0];
	
	/* create problem */
	lp = glp_create_prob();
	glp_set_prob_name(lp, "partition optimizer");
	glp_set_obj_dir(lp, GLP_MAX); // maximize
	
	/* fill problem */
	char name[8];
	glp_add_rows(lp, nb_row);
	glp_add_cols(lp, nb_var);

	for (i = 1 ; i <= nb_row ; i++) {
		sprintf(name, "%d", i);
		glp_set_row_name(lp, i, name);
		glp_set_row_bnds(lp, i, GLP_FX, X[i-1], X[i-1]);
	}

	for (i = 1 ; i <= nb_var ; i++) {
		sprintf(name, "x%d", i);
		glp_set_col_name(lp, i, name);
		glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0.0);
		glp_set_obj_coef(lp, i, -1.0);
	}

	j = 0;
	for (i = 1 ; i <= nb_no0 ; i++) {
		while (table[j/nb_var][j%nb_var] == 0) {
			j++;
			if (j == nb_row * nb_var) { break; } 
		}
		if (j == nb_row * nb_var) { break; }
		ia[i] = j/nb_var+1, ja[i] = j%nb_var+1, ar[i] = table[j/nb_var][j%nb_var];
		j++;
	}

	for (i = 0 ; i < nb_row ; i++) {
		free(table[i]);
	} free(table);

	/* load the problem */
	glp_load_matrix(lp, nb_no0, ia, ja, ar);
	
	/* config parameters */
	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_ERR;

	/* solve problem */
	glp_simplex(lp,&parm);
	
	/* recover and display results */
	// double z = -glp_get_obj_val(lp);
	// printf("minimum number of swap %.2f\n", z);
	
	double * x = (double *) calloc(nb_var, sizeof(double));
	for (i = 0 ; i < nb_var ; i++) {
		x[i] = glp_get_col_prim(lp, i+1);
	}

	double ** out = (double **) malloc(n * sizeof(double *));
	k = 0;
	for (i = 0 ; i < n ; i++) {
		out[i] = (double *) calloc(n, sizeof(double));
		for (j = 0 ; j < n ; j++) {
			if (B[i][j]) {
				out[i][j] = x[k];
				k++;
			}
		}
	}
	
	free(x);

	glp_delete_prob(lp);
	glp_free_env();

	return out;
}