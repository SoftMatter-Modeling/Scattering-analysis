#ifndef EIGENV3BY3
#define EIGENV3BY3

/*
CHANGE LOG
   Feb 01, 2019: Started this code
OBJECTIVE
   3x3 Eigenvalue calculation from numerical recipes
*/

void nrerror(char error_text[]);
double *vectort(long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
void free_vector(double *v, long nl, long nh);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void jacobi(double **a, int n, double d[], double **v, int nrot);
void eigsrt(double d[], double **v, int n);

#endif
