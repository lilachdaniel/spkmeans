#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

void print_debug(double **T, int N, int k){
    int i, j;
    for(i = 0; i < N; i++){
        for(j = 0; j < k; j++)
            printf("%f,", T[i][j]);
        printf("\n");
    }
}


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
    int i, j;
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
    double **ln_mat, **d_pow;

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
    int i, j;
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
void free_mat(double **mat, int rows){
    int i;

    for (i = 0; i < rows; ++i) {
        free(mat[i]);
    }

    free(mat);
}


/*************
 * Jacobi
 * */

/* Recieves Matrix A and dimensions and pointer to 1D eigenvalues Array
   Returns Array of eigenvectors and places in return_eigvals the eigenvalues */
double **Jac(double **A, int num_cols, int num_rows){
    int idx, sub_idx, i, j, loop = 0;
    double ** P, ** V, **return_array, ** A_tag, **temp_mat, c, s;
    double convergence = 1, prev_off = 0, new_off;
    
    /* Initiating V as Identity Matrix*/
    V = (double**)calloc(num_rows, sizeof(double*));
    assert(V);
    for(idx = 0; idx < num_rows; idx++){
        V[idx] = (double*)calloc(num_cols, sizeof(double));
        assert(V[idx]);
        V[idx][idx] = 1;
    }
    
    /* Initiating P */
    P = (double**)calloc(num_rows, sizeof(double*));
    assert(P);
    for(idx = 0; idx < num_rows; idx++){
        P[idx] = (double*)calloc(num_cols, sizeof(double));
        assert(P[idx]);
    }
    /* Initiating A_tag */
    A_tag = (double**)calloc(num_rows, sizeof(double*));
    assert(A_tag);
    for(idx = 0; idx < num_rows; idx++){
        A_tag[idx] = (double*)calloc(num_cols, sizeof(double));
        assert(A_tag[idx]);
    }
    /* Initiating temp_mat */
    temp_mat = (double**)calloc(num_rows, sizeof(double*));
    assert(temp_mat);
    for(idx = 0; idx < num_rows; idx++){
        temp_mat[idx] = (double*)calloc(num_cols, sizeof(double));
        assert(temp_mat[idx]);
    }
    /* Calulating initial off value of A */
    for(idx = 0; idx < num_rows; idx++){
        for(sub_idx = 0; sub_idx < num_cols; sub_idx++){
            if(idx != sub_idx)
                prev_off += A[idx][sub_idx]*A[idx][sub_idx];
        }
    }
    /*
    // printf("A:\n");
    // print_debug(A, num_rows, num_cols);
    // printf("V:\n");
    // print_debug(V, num_rows, num_cols);
    */

    /* Main loop of pivot. finding A_tag and changing A. 
       Finding P and multiplying to find V.
       Stopping after convergence is smaller than EPS or max rotations*/
    while ((convergence > EPS && loop++ < MAX_JAC_IT)){
        /* printf("in line %d, loop = %d\n", __LINE__, loop); */

        find_Rotation_Matrix(A, P, num_rows, &i, &j, &c, &s);
        
        
        /* fast_Mult(V, num_rows, i, j, c, s); */
        square_mat_mult(V, P, temp_mat, num_rows);

        new_off = construct_A_tag(A, A_tag, i, j, c, s, num_rows);

        convergence = prev_off - new_off;
        prev_off = new_off;

        if(is_diagonal_matrix(A, num_rows)){
            convergence = 0;
        }
        /*
        // printf("in loop %d convergence = %f\n", loop, convergence);
        // scanf("continue?%d", &idx);
        // printf("A:\n");
        // print_debug(A, num_rows, num_cols);
        // printf("V:\n");
        // print_debug(V, num_rows, num_cols);
        */

    }
    
    /* Place EigenValues in return_eigvals, free memory and return V */    
    /* NEED TO SORT EIGENVALUES!! */
    /* Initiating return_array */
    return_array = (double**)calloc(num_rows + 1, sizeof(double*));
    assert(return_array);
    for(idx = 0; idx < num_rows + 1; idx++){
        return_array[idx] = (double*)calloc(num_cols, sizeof(double));
        assert(return_array[idx]);
        for(sub_idx = 0; sub_idx < num_cols; sub_idx++)
            return_array[idx][sub_idx] = (idx == 0) ? A[sub_idx][sub_idx] : V[idx - 1][sub_idx];
    }
    free_mat(temp_mat, num_rows);
    free_mat(P, num_rows); 
    free_mat(A_tag, num_rows); 
    free_mat(V, num_rows); 

    return return_array;
}

/* Recieves pointer to A and P and constructs Rotation Matrix in P and
   values of c and s and values i and j coordinates of abs max of A */
void find_Rotation_Matrix(double **A, double **P, int num_rows, int *i, int *j, double *c, double *s){
    int idx_i, idx_j;
    double t, sign, theta;
    double abs_max = log(0);

    /* Finding absolute max value of A off diagonal */
    for(idx_i = 0; idx_i < num_rows; idx_i++){
        for(idx_j = idx_i+1; idx_j < num_rows; idx_j++){
            if (fabs(A[idx_i][idx_j]) > abs_max){
                abs_max = fabs(A[idx_i][idx_j]);
                *i = idx_i;
                *j = idx_j;
            }
        }
    }

    /* Finding c and s */
    theta = (A[*j][*j] - A[*i][*i]) / (2 * A[*i][*j]);
    sign = (theta == 0) ? 1 : fabs(theta)/theta;
    t = sign / (fabs(theta) + sqrt(theta*theta + 1));
    *c = 1 / (sqrt(t*t + 1));
    *s = t*(*c);

    /* Constructing P */
    for(idx_i = 0; idx_i < num_rows; idx_i++){
        for(idx_j = 0; idx_j < num_rows; idx_j++){
            P[idx_i][idx_j] = (idx_i == idx_j) ? 1 : 0; 
        }
    }
    P[*i][*i] = *c;
    P[*j][*j] = *c;
    P[*i][*j] = *s;
    P[*j][*i] = -*s;
}

/* Calculates A' using A_tag as temporary space and places in A using i, j cords of abs max off-
   diag element of A and c, s calculated in Find_Rotation_Matrix and returns convergance off value */
double construct_A_tag(double **A, double **A_tag, int i, int j, double c, double s, int num_rows){
    int r, ind, sub_ind;
    double A_tag_off = 0;
    for(r = 0; r < num_rows; r++){
        A_tag[r][i] = (r != i && r != j) ? c*A[r][i] - s*A[r][j] : A[r][i];
        A_tag[r][j] = (r != i && r != j) ? c*A[r][j] + s*A[r][i] : A[r][j];

        A_tag[i][r] = (r != i && r != j) ? c*A[r][i] - s*A[r][j] : A[r][i];
        A_tag[j][r] = (r != i && r != j) ? c*A[r][j] + s*A[r][i] : A[r][j];
    }
    A_tag[i][i] = c*c*A[i][i] + s*s*A[j][j] - 2*s*c*A[i][j];
    A_tag[j][j] = s*s*A[i][i] + c*c*A[j][j] + 2*s*c*A[i][j]; 
    A_tag[i][j] = 0; /* (c*c - s*s)*A[i][j] + s*c*(A[i][i]-A[j][j]); */
    A_tag[j][i] = 0;

    /* Copying A_tag to A */
    for(r = 0; r < num_rows; r++){
        A[r][i] =  A_tag[r][i];
        A[r][j] = A_tag[r][j];

        A[i][r] = A_tag[i][r];
        A[j][r] = A_tag[j][r];
    }

    /* Calulating convergence off value of new A */
    for(ind = 0; ind < num_rows; ind++){
        for(sub_ind = 0; sub_ind < num_rows; sub_ind++){
            if(ind != sub_ind)
                A_tag_off += A[ind][sub_ind]*A[ind][sub_ind];
        }
    }
    return A_tag_off;
}

/* function to multiply two square matrices then copy in to first */
void square_mat_mult(double **first, double **second, double **result, int N){
    int i, j, k;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            for (k = 0; k < N; k++)
                result[i][j] += first[i][k] * second[k][j];

    for (i = 0; i < N; i++) 
        for (j = 0; j < N; j++) 
            first[i][j] = result[i][j];
}

/* Multiply */
void fast_Mult(double **V, int num_rows, int i, int j, double c, double s){
    double i_elem_in_row, j_elem_in_row;
    int index;
    for(index = 0; index < num_rows; index++){        
        i_elem_in_row = c*V[index][i] - s*V[index][j];
        j_elem_in_row = s*V[index][i] + c*V[index][j];
        V[index][i] = i_elem_in_row;
        V[index][j] = j_elem_in_row;
    }
}

/* Recieves square matrix A and dimension n 
   Returns 1 if A is diagonal otherwise 0 */
int is_diagonal_matrix(double **A, int n){
    int d_idx, d_sub_idx;
    for(d_idx = 0; d_idx < n; d_idx++){
        for(d_sub_idx = 0; d_sub_idx < n; d_sub_idx++){
            if(d_idx != d_sub_idx && A[d_idx][d_sub_idx] != 0){
                return 0;
            }
        }
    }
    return 1;
}