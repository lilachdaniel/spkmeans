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
void free_mat(double **mat, int rows);


double **i_minus_mat(double **mat, int n);

// void display_mat(double **matrix, int num_rows, int num_cols);

/* Declarations for Jacobi) */
static double **Jac(double **A, int num_cols, int num_rows);
void find_Rotation_Matrix(double **A, double **P, int num_rows, int *i, int *j, double *c, double *s);
double construct_A_tag(double **A, double **A_tag, int i, int j, double c, double s, int num_rows, double prev_off);
void fast_Mult(double **V, int num_rows, int i, int j, double c, double s);
int is_diagonal_matrix(double **A, int n);

