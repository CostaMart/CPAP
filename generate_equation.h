#ifndef GENERATE_EQUATION
#define GENERATE_EQUATION

#include <stdbool.h>

bool isPositiveSymmetric(double **matrix, int n);
double **mediumMatrix(double **m1, double **m2, int order);
double *generateRandomLinearSystem(int order);
double *solutionVector(int lunghezza);

#endif //