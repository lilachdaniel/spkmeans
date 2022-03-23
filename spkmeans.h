#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

double **wam(double** vectors, int n, int d);
double calc_weight(double *x1, double *x2, int d);
double calc_l2(double *x1, double *x2, int d);
double **lnorm(double **d_mat, double **w_mat, int n);
double **mult_diag_mat_diag(double **diag, double **mat, int n);
double **ddg(double **w_mat, int n);
double **pow_diag_mat(double **diag, int n);
void free_mat(double **mat, int n);
double **i_minus_mat(double **mat, int n);