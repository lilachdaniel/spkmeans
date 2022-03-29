#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "spkmeans.h"

#define GOAL_SPK 0
#define GOAL_WAM 1
#define GOAL_DDG 2
#define GOAL_LNORM 3
#define GOAL_JACOBI 4

static PyObject* general_capi(PyObject *self, PyObject *args);
double **form_T(double **vectors, int *k, int N, int d);
double** renormalize(double **U, int N, int k);
int heuristic(int k, double **J, int N);
int cmp(const void * a, const void * b);

static double **fit_c(int k, int max_it, double **centroids, double **vectors, int vec_size, int num_vecs, double eps);
int find_cent_ind(double *v, double **centroids, int k, int vec_size);
double *add_vecs(double *v1, double *v2, int vec_size);
int update_centroids(int k, int vec_size, double **centroids, double **sums, int *counts, double eps);
double delta_norm_pow2(double *v1, double *v2, int vec_size);

/********************
 * Receives PyObject which are the vectors
 * Transfers to C and and runs form_T
 * Transfers back to PyObject and 
 * Returns T
 * ***************/
static PyObject* general_capi(PyObject *self, PyObject *args){
    int k, N, d, i, j, goal, return_cols;
    double **result, **W, **D, **vectors;
    PyObject *py_vectors;
    PyObject *item;
    PyObject *sub_item;
    PyObject *pyresult;
    PyObject *sub_pyresult;

    if(!PyArg_ParseTuple(args, "Oiiii", &py_vectors, &k, &N, &d, &goal)){
        return NULL;
    }
    /* Allocating memory for vectors */
    vectors = (double**)malloc(N*sizeof(double*));
    if(vectors==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }

    for(j = 0; j < N; j++){
        vectors[j] = (double*)malloc(d*sizeof(double));
        if(vectors[j]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }

    /* Transfering list of lists from PyObject py_vectors to double ** vectors */
    for (i = 0; i < N; i++){
        item = PySequence_GetItem(py_vectors, i);
        for(j = 0; j < d; j++){
            sub_item = PySequence_GetItem(item, j);
            vectors[i][j] = PyFloat_AsDouble(sub_item);
        }
    }

    /* for each goal the columns of the return matrix is different. */
    return_cols = N; 
    
    /* Run the C function */
    switch (goal)
    {
    case GOAL_SPK:
        result = form_T(vectors, &k, N, d); 
        return_cols = k;
        break;
    
    case GOAL_WAM:
        result = wam(vectors, N, d);
        break;
    case GOAL_DDG:
        W = wam(vectors, N, d);
        result = ddg(W, N);

        free_mat(W, N);
        break;

    case GOAL_LNORM:
        W = wam(vectors, N, d);
        D = ddg(W, N);
        result = lnorm(D, W, N);
        
        free_mat(W, N);
        free_mat(D, N);
        break;

    case GOAL_JACOBI:
        result = Jac(vectors, N, d);
        /* trying lilach's way where result is one dimension bigger so I must increase N */
        N += 1;
        break;

    default:
        printf("INVALID GOAL FROM PYTHON IN C API");
        exit(1);
    }
    
    
    pyresult = PyList_New(N); /* prepare list of lists to return */

    for (i = 0; i < N; i++) {
        sub_pyresult = PyList_New(return_cols);
        for(j = 0; j < return_cols; j ++){
            PyList_SetItem(sub_pyresult, j, Py_BuildValue("d", result[i][j]));
        }
        PyList_SetItem(pyresult,i, Py_BuildValue("O", sub_pyresult));
    }
 
    
    free_mat(result, N);

    /* Freeing vectors 
    * (if goal was Jacobi I have increased N so I must decrease before free vectors)*/
    if(goal == GOAL_JACOBI)
        N -=1;

    free_mat(vectors, N);

    return pyresult;
}

/****************
 * Kmeans Methods from HW2
 * *************/

static double **fit_c(int k, int max_it, double **centroids, double **vectors, int vec_size, int num_vecs, double eps){
    int i, j, l, cent_ind, smaller_than_e;
    double **sums;
    int *counts;

    /* Initializing sums and counts and allocating memory with starting values 0*/
    counts = (int*)malloc(k*sizeof(int));
    if(counts==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    sums = (double**)malloc(k * sizeof (double *));
    if(sums==NULL){
        free(counts);
        printf("An Error Has Occurred");
        exit(1);
    }

    for(j = 0; j < k; j++){
        sums[j] = (double*)malloc(vec_size*sizeof(double));
        if(sums[j]==NULL){
            free(counts);
            free(sums);
            printf("An Error Has Occurred");
            exit(1);
        }
    }


    /* K-means Cluster Dividing*/
    for(i = 0; i < max_it; i++){
        /* Initialize lists to calculate new centroids later*/
        for(j = 0; j < k; j++){
            counts[j] = 0;
            for(l = 0; l < vec_size; l++){
                sums[j][l] = 0;
            }
        }
        
        for(j = 0; j < num_vecs; j++){
            cent_ind = find_cent_ind(vectors[j], centroids, k, vec_size);

            /* Update sums and counts*/
            sums[cent_ind] = add_vecs(sums[cent_ind], vectors[j], vec_size);
            counts[cent_ind]++;
        }

        smaller_than_e = update_centroids(k, vec_size, centroids, sums, counts, eps);
        if (smaller_than_e){
            break;
        } 
    }

    free(counts);

    for(i = 0; i < k; i++){
        free(sums[i]);
    }
    free(sums); 

    return centroids;
}


double delta_norm_pow2(double *v1, double *v2, int vec_size){
    double sum = 0;
    int p;
    for (p = 0; p < vec_size; p++){
        sum += (v1[p]-v2[p])*(v1[p]-v2[p]);
    }
    return sum;
}


int find_cent_ind(double *v, double **centroids, int k, int vec_size){
    int min_ind = 0;
    int ind;
    double delt;
    double min_delta = delta_norm_pow2(v, centroids[0], vec_size);
    for (ind = 0; ind < k; ind++){
        delt = delta_norm_pow2(v, centroids[ind], vec_size);
        if (delt < min_delta){
            min_delta = delt;
            min_ind = ind;
        }
    }
    return min_ind;
}

int update_centroids(int k, int vec_size, double **centroids, double **sums, int *counts, double eps) {
    int smaller_than_e = 1; /*true = 1, false = 0*/
    int i,j,index;
    double *prev_cent;
    double delta_norm;
    prev_cent = (double*)malloc(vec_size*sizeof(double)); /* Allocating space for 1d array of a vector*/
    if(prev_cent==NULL){
            free(counts);
            for(i = 0; i < k; i++){
                free(sums[i]);
            }
            free(sums);
            printf("An Error Has Occurred");
            exit(1);
        }

    /*update all k centroids*/
    for (i = 0; i < k; i++) {
        for (index = 0; index < vec_size; index++){/* want to make copy*/
            prev_cent[index] = centroids[i][index];
        } 
        for (j = 0; j < vec_size; j++) {
            centroids[i][j] = sums[i][j] / counts[i];
        }

        /*update smaller_than_e*/
        delta_norm = pow(delta_norm_pow2(prev_cent, centroids[i], vec_size), 0.5);

        if (delta_norm >= eps) {
            smaller_than_e = 0;
        }
    }
    free(prev_cent);
    return smaller_than_e;
}

double *add_vecs(double *v1, double *v2, int vec_size) {
    int i;
    double *res = (double*)malloc(vec_size*sizeof(double ));

    for (i = 0; i < vec_size; i++) {
        res[i] = v1[i] + v2[i];
    }

    return res;
}

//receive python object (=the initial centroids) and returns python object (=the final centroids)
static PyObject* fit_capi(PyObject *self, PyObject *args) {
    
    int i, j, k, max_it, vec_size, num_vecs;
    PyObject *pyresult;
    PyObject *sub_pyresult;
    double **fit_result;
    PyObject *py_centroids;
    PyObject *py_vectors;
    double **centroids;
    double **vectors;
    double eps;
    PyObject *item;
    PyObject *sub_item;


 
    if(!PyArg_ParseTuple(args, "iiOOiid", &k, &max_it, &py_centroids, &py_vectors, &vec_size, &num_vecs, &eps)) {
        return NULL;
    }

    
    
   /* Allocating memory for vectors and centroids */
    vectors = (double**)malloc(num_vecs*sizeof(double*));
    if(vectors==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(j = 0; j < num_vecs; j++){
        vectors[j] = (double*)malloc(vec_size*sizeof(double));
        if(vectors[j]==NULL){

            printf("An Error Has Occurred");
            exit(1);
        }
    }
    centroids = (double**)malloc(k*sizeof(double*));
    if(centroids==NULL){
        printf("An Error Has Occurred");
        exit(1);
    }
    for(j = 0; j < k; j++){
        centroids[j] = (double*)malloc(vec_size*sizeof(double));
        if(centroids[j]==NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }
    /* Transfering list of lists from PyObject py_vectors to double ** vectors */
    for (i = 0; i < num_vecs; i++){
        item = PySequence_GetItem(py_vectors, i);
        for(j = 0; j < vec_size; j++){
            sub_item = PySequence_GetItem(item, j);
            vectors[i][j] = PyFloat_AsDouble(sub_item);
        }
    }
    /* Transfering list of lists from PyObject py_centroids to double ** centriods */
    for (i = 0; i < k; i++){
        item = PySequence_GetItem(py_centroids, i);
        for(j = 0; j < vec_size; j++){
            sub_item = PySequence_GetItem(item, j);
            centroids[i][j] = PyFloat_AsDouble(sub_item);
        }
    }

    pyresult = PyList_New(k);
   
    fit_result = fit_c(k, max_it, centroids, vectors, vec_size, num_vecs, eps);
   
    for (i = 0; i < k; i++) {
        sub_pyresult = PyList_New(vec_size);
        for(j = 0; j < vec_size; j ++){
            PyList_SetItem(sub_pyresult, j, Py_BuildValue("d", fit_result[i][j]));
        }
        PyList_SetItem(pyresult,i, Py_BuildValue("O", sub_pyresult));
    }

    for(i = 0; i < k; i++){
        free(centroids[i]);
    }
    free(centroids);

    for ( i = 0; i < num_vecs; i++)
    {
        free(vectors[i]);
    }
    free(vectors);
    

    return pyresult;
}


/****************
 * C - Python API Methods
 * *************/


static PyMethodDef capiMethods[] = {
        {"general_capi",
         (PyCFunction) general_capi,
         METH_VARARGS,
         PyDoc_STR("Function according to goal")},
         {"fit",
         (PyCFunction) fit_capi,
         METH_VARARGS,
         PyDoc_STR("k-means algorithm with randomized centroids initialization")},
         {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef= {
        PyModuleDef_HEAD_INIT,
        "spkmeansmodule",
        NULL,
        -1,
        capiMethods
};

PyMODINIT_FUNC
PyInit_spkmeansmodule(void) {
    PyObject *m;
    m= PyModule_Create(&moduledef);

    if(!m) {
        return NULL;
    }

    return m;
}

double **form_T(double **vectors, int *k, int N, int d){
    double **W, **L, **D, **J, **U, **T;
    int i, j;

    W = wam(vectors, N, d);
    D = ddg(W, N);
    L = lnorm(D, W, N);
    J = Jac(L, N, N);
    *k = heuristic(*k, J, N); 

    /* Allocating space for U and putting the eigenvectors in
    assuming that if they needed to be sorted we will sort them in Jac */
    U = (double**)malloc(N*sizeof(double*));
    assert(U);
    for(i = 0; i < N; i++){
        U[i] = (double*)malloc((*k)*sizeof(double));
        assert(U[i]);
        for(j = 0; j < *k; j++){
            U[i][j] = J[i + 1][j];
        }
    }

    T = renormalize(U, N, *k); /* function returns normalized U */
    
    free_mat(W, N);
    free_mat(D, N);
    free_mat(L, N);
    free_mat(J, N + 1);

   

    return T;

}

/***************
 * receives U and normalizes it
 * *************/
double** renormalize(double **U, int N, int k){
    int i, j;
    double sum = 0;
    for(i = 0; i < N; i++){
        for(j = 0; j < k; j++){
            sum += pow(U[i][j], 2);
        }
        for(j = 0; j < k; j++){
            U[i][j] /= pow(sum, 0.5);
        }
    }

    return U;
}


/* meant to return k if k!=0 otherwise calculate k */
int heuristic(int k, double **J, int N){
    int i, max_ind = 0;
    double delta, delta_max = 0;
    if (k)
        return k;
    
    qsort(J[0], N, sizeof(double), cmp);

    for (i = 0; i < (int)(N/2); i++){
        delta = fabs(J[0][i]-J[0][i + 1]);
        
        if(delta > delta_max){
            delta_max = delta;
            max_ind = i;
        }
    }
    return max_ind + 1;
}

int cmp(const void * a, const void * b) {
    
    if ( ( *(double*)a - *(double*)b ) > 0 ) return 1;
    else if ( ( *(double*)a - *(double*)b ) < 0 ) return -1;
    else return 0;
    
}