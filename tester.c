#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "spkmeans.h"

int check_wam(double **vectors, int n, int d);
int check_calc_l2(double *x1, double *x2, int d, double expected_ans);
int check_lnorm(double **d, double **w, int n);
double **init_mat(int n, int d, int max_val);
int random_test(int max_val, int max_size, int max_iter);

int main() {
    int max_val = 50;
    int max_size = 7;
    int max_iter = 1;
    int i, n, d;
    double *a1, *a2, *b, *b1, *c, *c2, **y1, **y2, **y;

        /* check calc_l2 */

        a1 = (double *)malloc(sizeof(double)*5);
        assert(a1 != NULL && "tester.main: error in memory allocation");
        a1[0] = 4;
        a1[1] = 5;
        a1[2] = 20;
        a1[3] = -3;
        a1[4] = -9;

        a2 = (double *)malloc(sizeof(double)*5);
        assert(a2 != NULL && "tester.main: error in memory allocation");
        a2[0] = 1;
        a2[1] = 1;
        a2[2] = 10;
        a2[3] = 1;
        a2[4] = 0;

        b = (double *)malloc(sizeof(double)*7);
        assert(b != NULL && "tester.main: error in memory allocation");
        b[0] = -7;
        b[1] = -5;
        b[2] = -2;
        b[3] = -4;
        b[4] = -202;
        b[5] = 0;
        b[6] = 0;

        b1 = (double *) calloc(7, sizeof(double));
        assert(b1 != NULL && "tester.main: error in memory allocation");

        c = (double *) calloc(3, sizeof(double));
        assert(c != NULL && "tester.main: error in memory allocation");

        c2 = (double *) calloc(3, sizeof(double));
        assert(c2 != NULL && "tester.main: error in memory allocation");



    //    double a[] = {3, 4, 10, -4, -9}; /* 14.899664425751 */
    //    double b[] = {-7, -5, -2 ,-4, 202, 0 , 0}; /* 202.23253941935 */
    //    double c[] = {0, 0, 0}; /* 0 */
    //
    //    double a1[] = {4, 5, 20, -3, -9};
    //    double a2[] = {1, 1, 10, 1, 0};
    //
    //    double b1[7] = {0};
    //
    //    double c2[3] = {0};

    //    check_calc_l2(a1, a2, 5, 14.899664425751);
    //    check_calc_l2(b1, b, 7, 202.23253941935);

        check_calc_l2(c, c2, 3, 0);

        free(a1);
        free(a2);
        free(b);
        free(b1);
        free(c);
        free(c2);


//        /* check mult_diag_mat_diag */
//
//        y = (double **) malloc(sizeof(double *) * 3);
//        assert(y != NULL && "tester.main: error in memory allocation");
//        y1 = (double **) malloc(sizeof(double *) * 3);
//        assert(y1 != NULL && "tester.main: error in memory allocation");
//        y2 = (double **) malloc(sizeof(double *) * 3);
//        assert(y2!= NULL && "tester.main: error in memory allocation");
//
//        for (i = 0; i < 3; ++i) {
//            y[i] = (double *) malloc(sizeof(double *) * 3);
//            assert(y[i] != NULL && "tester.main: error in memory allocation");
//            y1[i] = (double *) malloc(sizeof(double *) * 3);
//            assert(y1[i] != NULL && "tester.main: error in memory allocation");
//            y2[i] = (double *) malloc(sizeof(double *) * 3);
//            assert(y2[i] != NULL && "tester.main: error in memory allocation");
//        }
//
//        y[0][0] = 5;
//        y[0][1] = 6;
//        y[0][2] = 1;
//        y[1][0] = 9;
//        y[1][1] = 14;
//        y[1][2] = 6;
//        y[2][0] = -16;
//        y[2][1] = -16;
//        y[2][2] = -6;
//
//
//        y1[0][0] = 1;
//        y1[0][1] = 2;
//        y1[0][2] = -1;
//        y1[1][0] = 3;
//        y1[1][1] = 2;
//        y1[1][2] = 0;
//        y1[2][0] = -4;
//        y1[2][1] = 0;
//        y1[2][2] = 2;
//
//        y2[0][0] = 3;
//        y2[0][1] = 4;
//        y2[0][2] = 2;
//        y2[1][0] = 0;
//        y2[1][1] = 1;
//        y2[1][2] = 0;
//        y2[2][0] = -2;
//        y2[2][1] = 0;
//        y2[2][2] = 1;
//
//
//        //    double x1[4][3] = {{4, 6, -2}, {-1, 1, 3}, {-3, 2, 0}, {5, 7, 8}};
//        //    double x2[3][2] = {{4, 6}, {-1, 1}, {-3, 2}};
//        //    double x[4][2] = {{16, 26}, {-14, 1}, {-14, -16}, {-11, 53}}; /* result */
//
//    //    double y1[3][3] = {{1, 2, -1}, {3, 2, 0}, {-4, 0, 2}};
//    //    double y2[3][3] = {{3, 4, 2}, {0, 1, 0}, {-2, 0, 1}};
//    //    double y[3][3] = {{5, 6, 1}, {9, 14, 6}, {-16, -16, -6}}; /* result */
//
//        //    double z1[2][4] = {{1, 6, 3, 18}, {7, 12, 0, 23}};
//        //    double z2[4][3] = {{2, 12, 0}, {3, 4, 5}, {23, 15, 1}, {0, 0, 3}};
//        //    double z[2][3] = {{89, 81, 87}, {50, 132, 129}}; /* result */
//
//        check_mult_mat(y1, y2, 3, y);
//
//        free_mat(y, 3);
//        free_mat(y1, 3);
//        free_mat(y2, 3);

        /* random matrix test for wam and lnorm */
        srand(time(0));
        random_test(max_val, max_size, max_iter);

    return 0;
}

/* Tests wam and lnorm functions on random matrices */
int random_test(int max_val, int max_size, int max_iter) {
    double **obs_mat, **w_mat, **d_mat;
    int n, d, i;

    for (i = 0; i < max_iter; ++i) {
        /* choose n and d */
        n = rand() % max_size;
        d = rand() % max_size;

        /* init matrix */
        obs_mat = init_mat(n, d, max_val);
        display_mat(obs_mat, n, d);

        /* check wam */
        check_wam(obs_mat, n, d);

        /* check lnorm */
        w_mat = wam(obs_mat, n, d);
        d_mat = ddg(w_mat, n);
        check_lnorm(d_mat, w_mat, n);

        /* free memory */
        free_mat(w_mat, n);
        free_mat(obs_mat, n);
        free_mat(d_mat, n);

    }
    return 0; /* everything went fine! */
}

/* Initializes n random observations of size d */
double **init_mat(int n, int d, int max_val){
    int i, j;
    double **obs_mat;

    obs_mat = (double **)malloc(sizeof(double*) * n);
    assert(obs_mat != NULL && "tester.init_mat: error in memory allocation");

    for (i = 0; i < n; ++i) {
        obs_mat[i] = (double *)malloc(sizeof(double) * d);
        assert(obs_mat[i] != NULL && "tester.init_mat: error in memory allocation");

        for(j = 0; j < d; ++j) {
            obs_mat[i][j] = (double)(rand() % max_val);
            if (rand() % 1 == 1) {
                obs_mat[i][j] = - obs_mat[i][j];
            }
        }
    }

    return obs_mat;
}

/* Tests wam function */
int check_wam(double** vectors, int n, int d) {
    int i, j;
    double weight;
    double **w_mat = wam(vectors, n, d);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            /* check non-negativity */
            if (w_mat[i][j] < 0) {
                printf("Error in wam: w_mat%d%d is negative\n", i, j);
                display_mat(w_mat, n, n);
            }
            /* check symmetry */
            if (w_mat[i][j] != w_mat[j][i]) {
                printf("Error in wam: asymmetry in %d%d\n", i, j);
                display_mat(w_mat, n, n);
            }
            /* check self loops */
            if (w_mat[i][i] != 0) {
                printf("Error in wam: self loop in %d\n", i);
                display_mat(w_mat, n, n);
            }
            /* check weight */
            weight = calc_weight(vectors[i], vectors[j], d);
            if (i != j && w_mat[i][j] != weight) {
                printf("Error in wam: incorrect weight in w_mat%d%d. Expected weight: %f\n", i, j, weight);
                display_mat(w_mat, n, n);
            }
        }
    }

    free_mat(w_mat, n);

    return 0; /* Everything went fine! */
}

/* Tests calc_l2 function */
int check_calc_l2(double *x1, double *x2, int d, double expected_ans) {
    double actual_ans = calc_l2(x1, x2, d);
    if (actual_ans != expected_ans) {
        printf("Error in calc_l2: d == %d, expected value: %f, actual value: %f\n", d, expected_ans, actual_ans);
    }

    return 0; /* Everything went fine! */

}

/* Tests lnorm function
 * Note: this test is based on wikipedia "Laplacian matrix normalization" */
int check_lnorm(double **d, double **w, int n) {
    int i, j;
    double **ln = lnorm(d, w, n);

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j && d[i][i] != 0.0) {
                if (ln[i][j] != 1) {
                    printf("Error in lnorm: i == j == %d, ln[%d][%d] = %f, expected value: 1\n", i, i, j, ln[i][j]);
                }
            }
            else if (i != j && w[i][j] != 0.0) {
                if (ln[i][j] != ((-1/sqrt(d[i][i]*d[j][j])) * w[i][j])) {
                    printf("Error in lnorm: i == %d, j == %d, ln[%d][%d] = %f\n", i, j, i, j, ln[i][j]);
                }
            }
            else if (ln[i][j] != 0.0) {
                printf("Error in lnorm: i == %d, j == %d, ln[%d][%d] = %f expected value: 0\n", i, j, i, j, ln[i][j]);
            }
        }
    }

    free_mat(ln, n);

    return 0; /* Everything went fine! */
}
//
///* Receives a matrix and prints matrix to stdout */
//void display_mat(double **matrix, int num_rows, int num_cols) {
//    int i, j;
//    for (i = 0; i < num_rows; ++i) {
//        for(j = 0; j < num_cols; ++j) {
//            printf("%f ", matrix[i][j]);
//        }
//        printf("\n");
//    }
//}







