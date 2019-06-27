#ifndef _IO_H_
#define _IO_H_

#include <stdio.h>

void print_matrix(char * name, double ** matrix, int nb_row, int nb_column);
void print_vector(char * name, double * vector, int len);
void print_B(int ** B, int len);

#endif