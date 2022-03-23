#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "spkmeans.h"

int check_wam(double **vectors, int n, int d);
int check_calc_l2(double *x1, double *x2, int d, double expected_ans);
int check_lnorm(double **d, double **w, int n);
int check_mult_mat(double **a, double **b, int n, double **expected_ans);
void display_mat(double **matrix, int numberOfLines, int numberOfColumns);
double **init_mat(int n, int d, int max_val);

int main() {
    int max_val = 20;
    int max_size = 6;
    int max_iter = 1;
    int i, n, d;
    double *a1, *a2, *b, *b1, *c, *c2, **y1, **y2, **y, **mat, **w;

    //    /* check calc_l2 */
    //
    //    a1 = (double *)malloc(sizeof(double)*5);
    //    a1[0] = 4;
    //    a1[1] = 5;
    //    a1[2] = 20;
    //    a1[3] = -3;
    //    a1[4] = -9;
    //
    //    a2 = (double *)malloc(sizeof(double)*5);
    //    a2[0] = 1;
    //    a2[1] = 1;
    //    a2[2] = 10;
    //    a2[3] = 1;
    //    a2[4] = 0;
    //
    //    b = (double *)malloc(sizeof(double)*7);
    //    b[0] = -7;
    //    b[1] = -5;
    //    b[2] = -2;
    //    b[3] = -4;
    //    b[4] = -202;
    //    b[5] = 0;
    //    b[6] = 0;
    //
    //    b1 = (double *) calloc(7, sizeof(double));
    //
    //    c = (double *) calloc(3, sizeof(double));
    //    c2 = (double *) calloc(3, sizeof(double));
    //
    //
    ////    double a[] = {3, 4, 10, -4, -9}; /* 14.899664425751 */
    ////    double b[] = {-7, -5, -2 ,-4, 202, 0 , 0}; /* 202.23253941935 */
    ////    double c[] = {0, 0, 0}; /* 0 */
    ////
    ////    double a1[] = {4, 5, 20, -3, -9};
    ////    double a2[] = {1, 1, 10, 1, 0};
    ////
    ////    double b1[7] = {0};
    ////
    ////    double c2[3] = {0};
    //
    ////    check_calc_l2(a1, a2, 5, 14.899664425751);
    ////    check_calc_l2(b1, b, 7, 202.23253941935);
    //
    //    check_calc_l2(c, c2, 3, 0);
    //
    //    free(a1);
    //    free(a2);
    //    free(b);
    //    free(b1);
    //    free(c);
    //    free(c2);
    //
    //
    //    /* check mult_diag_mat_diag */
    //
    //    y = (double **) malloc(sizeof(double *) * 3);
    //    y1 = (double **) malloc(sizeof(double *) * 3);
    //    y2 = (double **) malloc(sizeof(double *) * 3);
    //
    //    for (i = 0; i < 3; ++i) {
    //        y[i] = (double *) malloc(sizeof(double *) * 3);
    //        y1[i] = (double *) malloc(sizeof(double *) * 3);
    //        y2[i] = (double *) malloc(sizeof(double *) * 3);
    //    }
    //
    //    y[0][0] = 5;
    //    y[0][1] = 6;
    //    y[0][2] = 1;
    //    y[1][0] = 9;
    //    y[1][1] = 14;
    //    y[1][2] = 6;
    //    y[2][0] = -16;
    //    y[2][1] = -16;
    //    y[2][2] = -6;
    //
    //
    //    y1[0][0] = 1;
    //    y1[0][1] = 2;
    //    y1[0][2] = -1;
    //    y1[1][0] = 3;
    //    y1[1][1] = 2;
    //    y1[1][2] = 0;
    //    y1[2][0] = -4;
    //    y1[2][1] = 0;
    //    y1[2][2] = 2;
    //
    //    y2[0][0] = 3;
    //    y2[0][1] = 4;
    //    y2[0][2] = 2;
    //    y2[1][0] = 0;
    //    y2[1][1] = 1;
    //    y2[1][2] = 0;
    //    y2[2][0] = -2;
    //    y2[2][1] = 0;
    //    y2[2][2] = 1;
    //
    //
    //    //    double x1[4][3] = {{4, 6, -2}, {-1, 1, 3}, {-3, 2, 0}, {5, 7, 8}};
    //    //    double x2[3][2] = {{4, 6}, {-1, 1}, {-3, 2}};
    //    //    double x[4][2] = {{16, 26}, {-14, 1}, {-14, -16}, {-11, 53}}; /* result */
    //
    ////    double y1[3][3] = {{1, 2, -1}, {3, 2, 0}, {-4, 0, 2}};
    ////    double y2[3][3] = {{3, 4, 2}, {0, 1, 0}, {-2, 0, 1}};
    ////    double y[3][3] = {{5, 6, 1}, {9, 14, 6}, {-16, -16, -6}}; /* result */
    //
    //    //    double z1[2][4] = {{1, 6, 3, 18}, {7, 12, 0, 23}};
    //    //    double z2[4][3] = {{2, 12, 0}, {3, 4, 5}, {23, 15, 1}, {0, 0, 3}};
    //    //    double z[2][3] = {{89, 81, 87}, {50, 132, 129}}; /* result */
    //
    //    check_mult_mat(y1, y2, 3, y);
    //
    //    free_mat(y, 3);
    //    free_mat(y1, 3);
    //    free_mat(y2, 3);
    //

    /* random matrix initialization */
    srand(time(0));

    for (i = 0; i < max_iter; ++i) {
        printf("here i= %d\n", i);

        n = rand() % max_size;
        d = rand() % max_size;
        printf("d = %d, n = %d\n", d, n);
        mat = init_mat(n, d, max_val);
        display_mat(mat, n, d);
        check_wam(mat, n, d);
        w = wam(mat, n, d);
        check_lnorm(ddg(w, n), w, n); /* free ddg */
        free_mat(w, n);
        free_mat(mat, n);
        printf("here after free i= %d\n", i);

    }

    return 0;
}

double **init_mat(int n, int d, int max_val){
    int i, j;
    double **mat;

    mat = (double **)malloc(sizeof(double*) * n);

    for (i = 0; i < n; ++i) {
        //            printf("i = %d\n", i);
        mat[i] = (double *) malloc(sizeof(double) * d);
        for(j = 0; j < d; ++j) {
            mat[i][j] = (double)(rand() % max_val);
            //                printf("%d%d\n", i, j);
            if (rand() % 1 == 1) {
                mat[i][j] = - mat[i][j];
            }
        }
    }

    return mat;
}


int check_wam(double** vectors, int n, int d) {
    int i, j;
    double weight;
    double **w = wam(vectors, n, d);
    printf("hi\n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            /* check non-negativity */
            if (w[i][j] < 0) {
                printf("Error in wam: w%d%d is negative\n", i, j);
                display_mat(w, n, n);
            }
            /* check symmetry */
            if (w[i][j] != w[j][i]) {
                printf("Error in wam: asymmetry in %d%d\n", i, j);
                display_mat(w, n, n);
            }
            /* check self loops */
            if (w[i][i] != 0) {
                printf("Error in wam: self loop in %d\n", i);
                display_mat(w, n, n);
            }
            /* check weight */
            weight = calc_weight(vectors[i], vectors[j], d);
            if (i != j && w[i][j] != weight) { /* self loops!!!!!!!! */
                printf("Error in wam: incorrect weight in w%d%d. Expected weight: %f\n", i, j, weight);
                display_mat(w, n, n);
            }
            printf("in check wam (%d,%d)\n", i, j);
        }
    }
    free_mat(w, n);
    return 0; /* Everything went fine! */
}

int check_calc_l2(double *x1, double *x2, int d, double expected_ans) {
    double actual_ans = calc_l2(x1, x2, d);
    if (actual_ans != expected_ans) {
        printf("Error in calc_l2: d == %d, expected value: %f, actual value: %f\n", d, expected_ans, actual_ans);
    }

    return 0; /* Everything went fine! */

}

int check_lnorm(double **d, double **w, int n) { /* test based on wikipedia "Laplacian matrix normalization" */
    int i, j;
    double **ln = lnorm(d, w, n);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j && d[i][i] != 0 && ln[i][j] != 1) {
                printf("Error in lnorm: i == j == %d, ln[%d][%d] expected value: 1\n", i, i, i, j);
            }
            else if (i != j && w[i][j] != 0 && ln[i][j] != -1/sqrt(d[i][i]*d[j][j])) {
                printf("Error in lnorm: i == %d, j == %d, ln[%d][%d] = %f, expected value: "
                       "-1/sqrt(deg(v%d)*deg(v%d))\n", i, j, i, j, ln[i][j], i, j);
            }
            else if (ln[i][j] != 0) {
                printf("Error in lnorm: i == %d, j == %d, ln[%d][%d] expected value: 0\n", i, j, i, j);
            }
        }
    }

    free_mat(ln, n);

    return 0; /* Everything went fine! */
}

int check_mult_mat(double **a, double **b, int n, double **expected_ans) {
    int i, j;
    double **c = mult_diag_mat_diag(a, b, n);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (c[i][j] != expected_ans[i][j]) {
                printf("Error in mult_diag_mat_diag: n == %d , i == %d, j == %d\n", n, i, j);
                display_mat(c, n, n);
            }
        }
    }

    free_mat(c, n);

    return 0; /* Everything went fine! */
}

void display_mat(double **matrix, int numberOfLines, int numberOfColumns) {
    int row, columns;
    for (row=0; row<numberOfLines; row++)
    {
        for(columns=0; columns<numberOfColumns; columns++)
        {
            printf("%f ", matrix[row][columns]);
        }
        printf("\n");
    }
}







