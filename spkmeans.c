#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


enum Goal {
    WAM = 0,
    DDG = 1,
    LNORM = 2,
    JACOBI = 3
};

int main(int argc, char *argv[]){
    int vec_size, sum_vecs;
    double **vectors, **W, **D, **L, **J;
    enum Goal gl;
    int size_vec_amount_vecs[2];

    if (argc != 3) err(TRUE);
    if (!strcmp(argv[1], "wam")) gl = WAM;
    else if (!strcmp(argv[1], "ddg")) gl = DDG;
    else if (!strcmp(argv[1], "lnorm")) gl = LNORM;
    else if (!strcmp(argv[1], "jacobi")) gl = JACOBI;
    else err(TRUE);

    /* Finding length of vectors and amount of vectors*/
    find_lengths_and_amount(argv[2], size_vec_amount_vecs);
    vec_size = size_vec_amount_vecs[0];
    sum_vecs =  size_vec_amount_vecs[1];

    /* Reading vectors from file*/
    vectors = collect(argv[2], sum_vecs, vec_size);
    
    switch(gl){
        case(WAM):
            W = wam(vectors, sum_vecs, vec_size);
            print_mat(W, sum_vecs);
            free_mat(W, sum_vecs); 
            return 0;
        
        case(DDG):
            W = wam(vectors, sum_vecs, vec_size);
            D = ddg(W, sum_vecs);
        
            print_mat(D, sum_vecs);
        
            free_mat(W, sum_vecs);
            free_mat(D, sum_vecs);
            return 0;
        
        case(LNORM):
            W = wam(vectors, sum_vecs, vec_size);
            D = ddg(W, sum_vecs);
            L = lnorm(D, W, sum_vecs);

            print_mat(L, sum_vecs);
        
            free_mat(W, sum_vecs);
            free_mat(D, sum_vecs);
            free_mat(L, sum_vecs);
            return 0;
        
        case(JACOBI):
            J = Jac(vectors, sum_vecs, vec_size);         
            
            print_Jac(J, sum_vecs);

            free_mat(J, sum_vecs + 1);
            return 0;

        default:
            err(TRUE);
    }
    return 1;
}

/* Creates a matrix of vectors from given file */
double **collect(char *input_filename, int sum_vectors, int vec_size){
    int i = 0, j = 0;
    double tmp;
    double **vectors;

    FILE *ifp = NULL;
    ifp = fopen(input_filename, "r");
    if(!ifp){
        err(FALSE);
    }

    /* Allocating memory for array of vectors*/
    vectors = (double**)malloc(sum_vectors*sizeof(double*));
    if(!vectors){
        err(FALSE);
    }
    for (i = 0; i < sum_vectors; i++){
        vectors[i] = (double*)malloc(vec_size*sizeof(double));
        if(!vectors[i]){
            err(TRUE);
        }
    }

    /* Load file as doubles in to array vectors*/
    for (i = 0; i < sum_vectors; i++){
        for (j = 0; j < vec_size; j++){
            fscanf(ifp, "%lf", &tmp);
            vectors[i][j] = tmp;
            fgetc(ifp);
        }
    }

    fclose(ifp);

    return vectors;
}

/* Receives square matrix and prints it */
void print_mat(double ** mat, int size){
    int i, j;
    for (i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            printf("%.4f", mat[i][j]);
            if(j != size - 1) printf(",");
            else printf("\n");
        } 
    }
}

/* Reads file to find amount of vectors and their length */
void find_lengths_and_amount(char *input_filename, int *size_vec_amount_vecs){
    int sum_vectors = 0, sum_cords = 0;
    char c;

    FILE *ifp = NULL;
    ifp = fopen(input_filename, "r");
    if(!ifp){
        err(FALSE);
    }

    /* Finding size of vector and amount of vectors*/
    while ((c = fgetc(ifp)) != EOF) {
        if ( c == '\n' ){
            sum_vectors++;
            sum_cords++;
        }
        else if (c == ','){
            sum_cords++;
        }
    }
    fclose(ifp);

    size_vec_amount_vecs[0] = sum_cords/sum_vectors;
    size_vec_amount_vecs[1] = sum_vectors;
}

/* Receives Jacobi matrix and prints eigenvalues and eigenvectors */
void print_Jac(double ** mat, int size){
    int i, j;
    /* print eigenvalues */
    for (i = 0; i < size; i++){
        /* Changing -0.0000 to 0.0000 */
        if(mat[0][i] > -0.0001)
            printf("%.4f", fabs(mat[0][i]));
        else
            printf("%.4f", mat[0][i]);
        
        if(i != size - 1) printf(",");
        else printf("\n");
    }

    /* print eigenvectors */
    for (i = 1; i < size + 1; i++){
        for(j = 0; j < size; j++){
            printf("%.4f", mat[i][j]);
            if(j != size - 1) printf(",");
            else printf("\n");
        }
    }
}