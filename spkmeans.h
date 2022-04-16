#define MAX_JAC_IT 100
#define EPS 0.00001
#define FALSE 0
#define TRUE 1

typedef int bool;

double **wam(double** vectors, int n, int d);
double calc_weight(double *x1, double *x2, int d);
double calc_l2(double *x1, double *x2, int d);
double **lnorm(double **d_mat, double **w_mat, int n);
double **mult_diag_mat_diag(double **diag, double **mat, int n);
double **ddg(double **w_mat, int n);
double **pow_diag_mat(double **diag, int n);
void free_mat(double **mat, int rows);
void err(int is_invalid_input);
double **i_minus_mat(double **mat, int n);
double **Jac(double **A, int num_cols, int num_rows);
void find_Rotation_Matrix(double **A, double **P, int num_rows, int *i, int *j,
                          double *c, double *s);
double construct_A_tag(double **A, double **A_tag, int i, int j, double c,
                       double s, int num_rows);
void fast_Mult(double **V, int num_rows, int i, int j, double c, double s);
int is_diagonal_matrix(double **A, int n);
double **collect(char *input_filename, int sum_vectors, int vec_size);
void find_lengths_and_amount(char *input_filename, int *size_vec_amount_vecs);
void print_mat(double ** mat, int size);
void print_Jac(double ** mat, int size);