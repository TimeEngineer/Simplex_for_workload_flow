# simplex

How to use this library for partitionning :

First install `glpk` library, in linux you can type :

```sudo apt-get install glpk-utils libglpk-dev glpk-doc```

io.h contains basic input/output for matrix and vector :

- void print_matrix(char * name, double ** matrix, int nb_row, int nb_column);

- void print_vector(char * name, double * vector, int len);

- void print_B(int ** B, int len);

simplex.h contains the solution of your problem :

- double ** simplex_procedure(double * X, int ** B, int n);

Where :

- X is a vector of normalized value of the partition

- B is a matrix of boundary to know where are the neighbours of each partition

- n is the number of partition (lenght of X and B)

Also there is a test included where you can choose your parameters.

In order to run the test, you can write :

```make && ./main```.