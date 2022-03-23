#include "spkmeans.h"


/* Receives n observations of size d
 * Returns a nxn weighted adjacency matrix */
double **wam(double **vectors, int n, int d) {
    int i,j;
    double **w_mat = (double **)malloc(sizeof(double *) * n);
    assert(w_mat != NULL && "wam: error in memory allocation");

    for (i = 0; i < n; ++i){
        w_mat[i] = (double *) malloc(sizeof(double) * n);
        assert(w_mat[i] != NULL && "wam: error in memory allocation");
    }

    for (i = 0; i < n; ++i){
        for (j = i; j < n; ++j) {
            if (i == j) { /* self loops not allowed */
                w_mat[i][j] = 0;
            }
            else {
                w_mat[i][j] = calc_weight(vectors[i], vectors[j], d);
                w_mat[j][i] = w_mat[i][j]; /* supposed to be symmetric */
            }
        }
    }
    return w_mat;
}

/* Receives two observations of size d
 * Returns the weight of the connection between them */
double calc_weight(double *x1, double *x2, int d) {
    double l2, weight;

    l2 = calc_l2(x1, x2, d);
    weight = exp(-l2/2);

    return weight;
}

/* Receives two vectors of size d
 * Returns the Euclidean norm of their difference */
double calc_l2(double *x1, double *x2, int d){
    int i;
    double l2 = 0;
    double x;

    for (i = 0; i < d; ++i){
        x = pow(x1[i] - x2[i], 2);
        l2 += x;
    }

    return sqrt(l2);
}

/* Receives a nxn weighted adjacency matrix
 * Returns a nxn diagonal degree matrix */
double **ddg(double **w_mat, int n) {
    int i, j, sum;
    double **d_mat = (double **)malloc(sizeof(double *) * n);
    assert(d_mat != NULL && "ddg: error in memory allocation");

    for (i = 0; i < n; ++i) {
        d_mat[i] = (double *)malloc(sizeof(double) * n);
        assert(d_mat[i] != NULL && "ddg: error in memory allocation");

        d_mat[i][i] = 0;
        for (j = 0; j < n; ++j) {
            d_mat[i][i] += w_mat[i][j]; /* d_mat[i][i] = sum(w_mat[i]) */
            if (i != j) { /* diagonal matrix */
                d_mat[i][j] = 0;
            }
        }
    }

    return d_mat;
}

/* Receives a nxn diagonal degree matrix and a nxn weighted adjacency matrix
 * Returns a nxn normalized graph laplacian matrix */
double **lnorm(double **d_mat, double **w_mat, int n){
    int i, j;
    double **ln_mat, **d_pow;

    /*initialize Lnorm matrix*/
    ln_mat = (double **)malloc(sizeof(double *) * n);
    assert(ln_mat != NULL && "lnorm: error in memory allocation");

    for (i = 0; i < n; ++i){
        ln_mat[i] = (double *) malloc(sizeof(double) * n);
        assert(ln_mat[i] != NULL && "lnorm: error in memory allocation");
    }

    /* create I - D^-0.5 * W * D^-0.5 */
    d_pow = pow_diag_mat(d_mat, n); /* D^-0.5 */
    ln_mat = mult_diag_mat_diag(d_pow, w_mat, n); /* D^-0.5 * W * D^-0.5 */
    ln_mat = i_minus_mat(ln_mat, n); /* I - D^-0.5 * W * D^-0.5 */

    /* free temporary mem */
    free_mat(d_pow, n);

    return ln_mat;

}

/* Receives a nxn matrix
 * Returns a nxn matrix (I - mat) */
double **i_minus_mat(double **mat, int n) {
    int i, j;

    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++) {
            if(i == j) { /* I[i][i] == 1 */
                mat[i][j] = 1 - mat[i][j];
            }
            else { /* I[i][j] == 0 (i != j) */
                mat[i][j] = - mat[i][j];
            }
        }
    }
    return mat;
}

/* Receives a nxn diagonal matrix and a nxn matrix
 * Returns the multiplication (diag * mat * diag) */
double **mult_diag_mat_diag(double **diag, double **mat, int n){
    int i, j, k;
    double **mult;

    /* initialize mult matrix */
    mult = (double **)malloc(sizeof(double *) * n);
    assert(mult != NULL && "mult_diag_mat_diag: error in memory allocation");

    for (i = 0; i < n; ++i){
        mult[i] = (double *) malloc(sizeof(double) * n);
        assert(mult[i] != NULL && "mult_diag_mat_diag: error in memory allocation");
    }

    /* multiply */
    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j) {
            mult[i][j] = diag[i][i] * mat[i][j] * diag[j][j];
        }
    }

    return mult;
}

/* Receives a nxn diagonal matrix
 * Returns a nxn diagonal matrix (diag^-0.5) */
double **pow_diag_mat(double **diag, int n) {
    int i;
    double **mat = (double **)malloc(sizeof(double *) * n);
    assert(mat != NULL && "pow_diag_mat: error in memory allocation");

    for (i = 0; i < n; ++i) {
        mat[i] = (double *)calloc(n, sizeof(double));
        assert(mat[i] != NULL && "pow_diag_mat: error in memory allocation");

        mat[i][i] = pow(diag[i][i], -0.5);
    }

    return mat;
}

/* Receives a nxn matrix
 * Frees memory of that matrix */
void free_mat(double **mat, int n){
    int i;

    for (i = 0; i < n; ++i) {
        free(mat[i]);
    }

    free(mat);
}

/* Receives a matrix and prints matrix to stdout */
void display_mat(double **matrix, int num_rows, int num_cols) {
    int i, j;
    for (i = 0; i < num_rows; ++i) {
        for(j = 0; j < num_cols; ++j) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}