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

//static PyObject* general_capi(PyObject *self, PyObject *args);

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
    return_cols = d; 
    
    /* Run the C function */
    switch (goal)
    {
    case GOAL_SPK:
        result = form_T(vectors, &k, N, d); 
        return_cols = k;
        break;
    case GOAL_WAM:
        //result = wam(vectors, N, d);
        break;
    case GOAL_DDG:
        //W = wam(vectors, N, d);
        //result = ddg(W, N);
        /* Freeing W */
        for ( i = 0; i < N; i++)
        {
            free(W[i]);
        }
        free(W);
        break;

    case GOAL_LNORM:
        // W = wam(vectors, N, d);
        // D = ddg(W, N);
        // result = lnorm()

        /* Freeing W */
        for ( i = 0; i < N; i++)
        {
            free(W[i]);
        }
        free(W);
        /* Freeing D */
        for ( i = 0; i < N; i++)
        {
            free(D[i]);
        }
        free(D);
        break;

    case GOAL_JACOBI:
        // eigenvalues = (double*)malloc(N*sizeof(double));
        // assert(eigenvalues);  TRYING NEW WAY OF JAC
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
 
    /* Freeing result */
    for ( i = 0; i < N; i++)
    {
        free(result[i]);
    }
    free(result);

    /* Freeing vectors 
    * (if goal was Jacobi I have increased N so I must decrease before free vectors)*/
    if(goal == GOAL_JACOBI)
        N -=1;

    for ( i = 0; i < N; i++)
    {
        free(vectors[i]);
    }
    free(vectors);
   

    return pyresult;
}

static PyMethodDef capiMethods[] = {
        {"general_capi",
         (PyCFunction) general_capi,
         METH_VARARGS,
         PyDoc_STR("Function according to goal")},
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