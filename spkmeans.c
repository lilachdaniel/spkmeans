#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

double **collect(char *input_filename, int sum_vectors, int vec_size);
void find_lengths_and_amount(char *input_filename, int *size_vec_amount_vecs);
void print_mat(double ** mat, int size); 
void print_Jac(double ** mat, int size); 

enum Goal {
    WAM = 0,
    DDG = 1,
    LNORM = 2,
    JACOBI = 3
};

int main(int argc, char *argv[]){
    int vec_size, sum_vecs, *size_vec_amount_vecs;
    double **vectors, **W, **D, **L, **J;
    enum Goal gl;

    size_vec_amount_vecs = (int*)malloc(2*sizeof(int));
    assert(size_vec_amount_vecs);

    if(argc != 3) goto term_input;
    if(!strcmp(argv[1], "wam")) gl = WAM;
    else if(!strcmp(argv[1], "ddg")) gl = DDG;
    else if(!strcmp(argv[1], "lnorm")) gl = LNORM;
    else if(!strcmp(argv[1], "jacobi")) gl = JACOBI;
    else goto term_input;

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
            exit(0);
        
        case(DDG):
            W = wam(vectors, sum_vecs, vec_size);
            D = ddg(W, sum_vecs);
        
            print_mat(D, sum_vecs);
        
            free_mat(W, sum_vecs);
            free_mat(D, sum_vecs);
            exit(0);
        
        case(LNORM):
            W = wam(vectors, sum_vecs, vec_size);
            D = ddg(W, sum_vecs);
            L = lnorm(D, W, sum_vecs);

            print_mat(L, sum_vecs);
        
            free_mat(W, sum_vecs);
            free_mat(D, sum_vecs);
            free_mat(L, sum_vecs);
            exit(0);
        
        case(JACOBI):
            J = Jac(vectors, sum_vecs, vec_size);         
            
            print_Jac(J, sum_vecs);

            free_mat(J, sum_vecs + 1);
            exit(0);
        default:
            goto term_input;
            break;
    }



term_input:
    printf("Invalid Input!");
    exit(1);
}


double **collect(char *input_filename, int sum_vectors, int vec_size){
    int i = 0, j = 0;
    double tmp;
    double **vectors;

    FILE *ifp = NULL;
    ifp = fopen(input_filename, "r");
    if(ifp==NULL){
        printf("An Error Has Occurred");
        exit(1);
        }

    /* Allocating memory for array of vectors*/
    vectors = (double**)malloc(sum_vectors*sizeof(double*));
    if(vectors==NULL){
        printf("An Error Has Occurred");
        exit(1);
        }
    for (i = 0; i < sum_vectors; i++){
        vectors[i] = (double*)malloc(vec_size*sizeof(double));
        if(vectors[i]==NULL){
        printf("An Error Has Occurred");
        exit(1);
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

void print_mat(double ** mat, int size){
    int i, j;
    for (i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            printf("%.4f", mat[i][j]);
            if(j != size - 1) printf(",");
            else printf("\n");
        } 
    }
    printf("\n");
}

void find_lengths_and_amount(char *input_filename, int *size_vec_amount_vecs){

    int sum_vectors = 0, sum_cords = 0;
    char c;

    FILE *ifp = NULL;
    ifp = fopen(input_filename, "r");
    if(ifp==NULL){
        printf("An Error Has Occurred");
        exit(1);
        }

    /* Finding size of vector and amount of vectors*/
    while ( ( c = fgetc( ifp ) ) != EOF ) {
        if ( c == '\n' ){
            sum_vectors++;
            sum_cords++;
        }
        else if (c == ',')
        {
            sum_cords++;
        }
    }
    fclose(ifp);

    size_vec_amount_vecs[0] = sum_cords/sum_vectors;
    size_vec_amount_vecs[1] = sum_vectors;
}


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
    printf("\n");
}