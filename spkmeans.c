#include "spkmeans.h"


double **collect(char *input_filename, int sum_vectors, int vec_size);
void find_lengths_and_amount(char *input_filename, int *size_vec_amount_vecs);
void print_mat(double ** mat, int size); 
void print_Jac(double ** mat, int size); 

enum Goal {
    wam = 0,
    ddg = 1,
    lnorm = 2,
    jacobi = 3
};

int main(int argc, char *argv[]){
    int vec_size, sum_vecs, *size_vec_amount_vecs;

    size_vec_amount_vecs = (int*)malloc(2*sizeof(int));
    assert(size_vec_amount_vecs);
    double **vectors, **W, **D, **L, **J;

    if(argc != 3) goto term_input;
    enum Goal gl;
    if(!strcomp(argv[1], "wam")) gl = wam;
    else if(!strcomp(argv[1], "ddg")) gl = ddg;
    else if(!strcomp(argv[1], "lnorm")) gl = lnorm;
    else if(!strcomp(argv[1], "jacobi")) gl = jacobi;
    else goto term_input;

    /* Finding length of vectors and amount of vectors*/
    find_lengths_and_amount(argv[2], size_vec_amount_vecs);
    vec_size = size_vec_amount_vecs[0];
    sum_vecs =  size_vec_amount_vecs[1];

    /* Reading vectors from file*/
    vectors = collect(argv[2], sum_vecs, vec_size);
    
    switch(gl){
        case(wam):
            W = wam(vectors, sum_vecs, vec_size);
            print_mat(W, sum_vecs);
            free_mat(W, sum_vecs); 
            exit(0);
        
        case(ddg):
            W = wam(vectors, sum_vecs, vec_size);
            D = ddg(W, sum_vecs);
        
            print_mat(D, sum_vecs);
        
            free_mat(W, sum_vecs);
            free_mat(D, sum_vecs);
            exit(0);
        
        case(lnorm):
            W = wam(vectors, sum_vecs, vec_size);
            D = ddg(W, sum_vecs);
            L = lnorm(D, W, sum_vecs);

            print_mat(D, sum_vecs);
        
            free_mat(W, sum_vecs);
            free_mat(D, sum_vecs);
            free_mat(L, sum_vecs);
            exit(0);
        
        case(jacobi):
            J = Jac(vectors, sum_vecs, vec_size);         
            
            print_Jac(J, sum_vecs + 1);

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
    printf("[");
    for (i = 0; i < size; i++){
        printf("[");
        for(j = 0; j < size; j++){
            printf("%.3f", mat[i][j]);
            if(j != size - 1) printf(", ");
        }
        printf("]\n");
    }
    printf("]\n");
}

void print_Jac(double ** mat, int size){
    int i, j;

    /* print eigenvalues */
    printf("[");
    for (i = 0; i < size; i++){
        printf("%.4f", mat[0][i]);
        if(i != size - 1) printf(", ");
    }
    printf("]\n");
    /* print eigenvectors */
    printf("[");
    for (i = 1; i < size + 1; i++){
        printf("[");
        for(j = 0; j < size; j++){
            printf("%.4f", mat[i][j]);
            if(j != size - 1) printf(", ");
        }
        printf("]\n");
    }
    printf("]\n");
}