#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <stdio.h>            /* C input/output                       */
#include <stdlib.h>           /* C standard library                   */
#include <glpk.h>             /* GNU GLPK linear/mixed integer solver */

double ** simplex_procedure(double * X, int ** B, int n);

#endif