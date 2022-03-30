#ifndef MATRIX_CALC_C_SPK_ALGORITHM_H
#define MATRIX_CALC_C_SPK_ALGORITHM_H
#define verifyNotNULL(var) if((var)==NULL) {printf("An Error Has Occured"); exit(-1);}

/* structs: */
typedef struct ROTATION_MATRIX{
    double c;
    double s;
    int i;
    int j;
} rotation_mat;

typedef struct Jacobi_output{
    double *eigenValues;
    double **V;
} Jacobi_output;

/* functions: */
void free_contiguous_mat(double **mat);
double** calc_weighted_adjacency_matrix(double **dots, int d, int n);
double** calc_diagonal_degree_matrix(double **W, int n);
double **calc_L_norm(double **D, double **W, int n);
Jacobi_output *jacobi(double **A, int n);
double **calc_T(Jacobi_output *jacobiOutput, int n, int *k_pointer);
void print_2d_array(double **array, int row_num, int col_num);
void print_2d_array_transpose(double **array, int row_num, int col_num);
double **deepCopy2DArray(double **A, int row_num, int col_num);
void print_list(double *, int);
double **kmeans(double **dots,double **centroids, int k, int d, int n, int max_iter);

#endif
