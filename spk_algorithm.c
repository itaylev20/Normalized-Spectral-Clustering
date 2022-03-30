#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define EPSILON (1E-15)
#define MAX_JACOBI_ITERS 100
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
double calc_weight(const double *dot1, const double *dot2, int d);
double** alloc_nXm_matrix(int n, int m);
double** alloc_nXn_matrix(int n);
void free_contiguous_mat(double **mat);
double** create_identity_matrix(int n);
double** calc_weighted_adjacency_matrix(double **dots, int d, int n);
double** calc_diagonal_degree_matrix(double **W, int n);
double **calc_sqrt_minus_half_of_diagonal_mat(double **D, int n);
double **calc_L_norm(double **D, double **W, int n);
int *find_off_diag_max_abs_val(double **A, int n);
double **calc_A_tag(double **A_tag, double **A, int n, rotation_mat *P);
double calc_of_f_square(double **A, int n);
void multiply_rotation_matrix(double **V, rotation_mat *P,int n);
Jacobi_output *jacobi(double **A, int n);
rotation_mat *calc_rotation_mat(rotation_mat *P,double **A, int n);
int compare_eigenStruct(const void *a, const void *b);
double **calc_T(Jacobi_output *jacobiOutput, int n, int *k_pointer);

/* Kmeans functions: */
double calc_distance(const double *dot, const double *centroid, int d);
int find_nearest_centroid(double *dot, double **centroids, int k, int d);
void update_centroids(double **dots,const int *dots_location, double **new_centroids, int n, int d, int k);
int check_equals_2d_list(double **list1, double **list2, int row_num, int col_num);
void print_2d_array(double **array, int row_num, int col_num);
void print_2d_array_transpose(double **array, int row_num, int col_num);
double **deepCopy2DArray(double **A, int row_num, int col_num);
void print_list(double *, int);
double **kmeans(double **dots,double **centroids, int k, int d, int n, int max_iter);

/*-------------------------------- calculate weight: ------------------------------------------------------*/

double calc_weight(const double *dot1, const double *dot2, int d){
    double sum = 0;
    int i;
    for(i=0 ; i<d ; i++)
        sum += ((dot1[i] - dot2[i]) * (dot1[i] - dot2[i]));
    return exp(-((sqrt(sum))/2));
}

/*--------------------- allocating contiguous memory for a nXm matrix: ----------------------------------*/

double** alloc_nXm_matrix(int n, int m){
    int i;
    double **output_matrix = (double**) calloc(n ,sizeof(double*));
    double *matrix_data = (double*)calloc(n * m , sizeof(double));

    verifyNotNULL(output_matrix)
    verifyNotNULL(matrix_data)

    for (i = 0; i < n; i++)
        output_matrix[i] = matrix_data + i * m;

    return output_matrix;
}
/*------------------------- allocation memory for a nXn matrix: ----------------------------*/
double** alloc_nXn_matrix(int n){
    return alloc_nXm_matrix(n,n);
}

/*------------------------- free contiguous memory of a nXm matrix: ---------------------------*/
void free_contiguous_mat(double **mat){
    free (mat[0]);
    free(mat);
}

/*----------------------------- creating Identity matrix n: --------------------------------------*/
double** create_identity_matrix(int n){
    double** I=alloc_nXn_matrix(n);
    int i,j;

    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            I[i][j]= i==j ? 1 :0 ;
    return I;
}

/*-------------------------------  calculate weighted adjacency matrix: ----------------------------------*/

double** calc_weighted_adjacency_matrix(double **dots, int d, int n){
    /* n is the number of dots. d is the dimension of each dot. */
    int i,j;
    double  weight_ij;
    double **weighted_adj_mat;

    /* allocation memory for the output matrix: */
    weighted_adj_mat = alloc_nXn_matrix(n);

    /* filling the primary diagonal with zeros: */
    for (i=0 ; i<n ; i++){
        weighted_adj_mat[i][i] = 0;
    }

    /* filling the rest of the matrix with the weights: */
    for (i=0 ; i<n ; i++){
        for (j=i+1 ; j<n ; j++){  /* only upper triangular */
            weight_ij = calc_weight(dots[i], dots[j], d);
            weighted_adj_mat[i][j] = weight_ij;
            weighted_adj_mat[j][i] = weight_ij;
        }
    }


    return weighted_adj_mat;
}

/*------------------------------ calc_diagonal_degree_matrix --------------------------------------------------------*/

double** calc_diagonal_degree_matrix(double **W, int n){
    /* W is symmetric. n is the dimensions of W. */
    int i,j;
    double line_sum;
    double **diagonal_degree_mat;

    /* allocation memory for the output matrix: */
    diagonal_degree_mat = alloc_nXn_matrix(n);

    /* filling the output matrix: */
    for (i=0 ; i<n ; i++){
        line_sum = 0;
        for (j=0 ; j<n ; j++){
            line_sum += W[i][j];
        }
        diagonal_degree_mat[i][i] = line_sum;
    }

    return diagonal_degree_mat;
}

/* ---------------------------- calc_sqrt_minus_half_of_diagonal_mat: ----------------------------------------*/

double **calc_sqrt_minus_half_of_diagonal_mat(double **D, int n){
    /* D is symmetric. n is the dimensions of D. */
    int i;
    double **output_mat;

    /* allocation memory for the output matrix: */
    output_mat = alloc_nXn_matrix(n);

    /* filling the output matrix: */
    for (i=0 ; i<n ; i++){
        output_mat[i][i] = 1/(sqrt(D[i][i]));
    }

    return output_mat;
}

/* ------------------------------------- calc_L_norm: -----------------------------------------------------------*/

double **calc_L_norm(double **D, double **W, int n){
    /* W is symmetric and D is diagonal */
    int i,j;
    double Dii, Djj, res_ij;
    double **D_sqrt_minus_half, **L_norm;
    D_sqrt_minus_half = calc_sqrt_minus_half_of_diagonal_mat(D, n);
    L_norm = alloc_nXn_matrix(n);

    /* L_norm = I - (D_sqrt_minus_half * W * D_sqrt_minus_half): */
    for(i=0 ; i<n ; i++){
        Dii = D_sqrt_minus_half[i][i];
        for (j=i+1 ; j<n ; j++){
            Djj = D_sqrt_minus_half[j][j];
            res_ij = -(Dii * W[i][j] * Djj); /* the minus is because of the "I-" */
            L_norm[i][j] = res_ij;
            L_norm[j][i] = res_ij;
        }
        L_norm[i][i] = 1; /* because of the "I-" */
    }
    free_contiguous_mat(D_sqrt_minus_half);



    return L_norm;
}

/*//////////////////////////////////////// Jacobi Algorithm: //////////////////////////////////////////////////////*/
/* ----------------------- Jacobi main function: ---------------------------------*/

Jacobi_output *jacobi(double **A, int n) {
    double **A_tag, **V, **temp;
    double *eigenValues;
    Jacobi_output *output_struct;
    rotation_mat *P;
    int i = 1;

    A = deepCopy2DArray(A, n, n);
    V = create_identity_matrix(n);

    P = (rotation_mat *) malloc(sizeof(rotation_mat));
    verifyNotNULL(P)
    calc_rotation_mat(P, A, n);

    A_tag = calc_A_tag(alloc_nXn_matrix(n), A, n, P);




    multiply_rotation_matrix(V, P, n);
    while (calc_of_f_square(A, n) - calc_of_f_square(A_tag, n) > EPSILON && i <= MAX_JACOBI_ITERS) {
        temp = A;
        A = A_tag;
        calc_rotation_mat(P, A, n);
        A_tag=calc_A_tag(temp, A, n, P);
        multiply_rotation_matrix(V, P, n);
        i++;
    }

    output_struct = (Jacobi_output *) malloc(sizeof(Jacobi_output));
    eigenValues = (double *) calloc(n,sizeof(double) );
    verifyNotNULL(output_struct)
    verifyNotNULL(eigenValues)


    for (i = 0; i < n; i++)
        eigenValues[i] = A_tag[i][i];

    output_struct->eigenValues = eigenValues;
    output_struct->V = V;

    free(P);
    free_contiguous_mat(A_tag);
    free_contiguous_mat(A);


    return output_struct;
}

/* ----------------------- multiply rotation matrix: ---------------------------------*/

void multiply_rotation_matrix(double **V, rotation_mat *P, int n) {
    int r;
    double res_i, res_j;
    for (r = 0; r < n; r++) {
        /* Note: using temp vars to avoid sync issues */
        res_i = V[r][P->i] * P->c + V[r][P->j] * P->s * -1;
        res_j = V[r][P->i] * P->s + V[r][P->j] * P->c;
        V[r][P->i] = res_i;
        V[r][P->j] = res_j;
    }
}

/* ----------------------- Calculate A': ---------------------------------*/

/* Calculate A' by the guidance in Section 6 - Relation between A and A': */
double **calc_A_tag(double **A_tag, double **A, int n, rotation_mat *P) {
    int r, col;

    /* Copying A to A', but for each row, at indices i,j assign the formula as described in Section 6. */
    for (r = 0; r < n; r++) {
        for (col = 0; col < n; col++) {
            if (col == P->i) {
                A_tag[r][col] = P->c * A[r][P->i] - P->s * A[r][P->j];
                A_tag[col][r] = A_tag[r][col];
            } else if (col == P->j) {
                A_tag[r][col] = P->c * A[r][P->j] + P->s * A[r][P->i];
                A_tag[col][r] = A_tag[r][col];
            } else if (r != P->j && r !=
                                    P->i) /* To prevent override bc this is a symmetric matrix - and in the above if we are doing it also to the opposite indices TODO - rewrite this comment. */
                A_tag[r][col] = A[r][col];
        }
    }
    /* Overriding A_tag indices for the last 3 formulas in Section 6 (Indices ii,jj,ij,ji) */
    A_tag[P->i][P->i] = pow(P->c, 2) * A[P->i][P->i] + pow(P->s, 2) * A[P->j][P->j] - 2 * P->s * P->c * A[P->i][P->j];
    A_tag[P->j][P->j] = pow(P->s, 2) * A[P->i][P->i] + pow(P->c, 2) * A[P->j][P->j] + 2 * P->s * P->c * A[P->i][P->j];
    A_tag[P->i][P->j] = 0;
    A_tag[P->j][P->i] = 0;

    return A_tag;
}

/* ------------------------ find_off_diag_max_abs_val -----------------------------*/

int *find_off_diag_max_abs_val(double **A, int n){
    /*A is symmetric. n is the dimensions of A. Returns a list of [i,j]. */
    int i, j;
    double max_val = 0, abs_curr;
    int *ij_list = (int*) calloc(2,sizeof(int));
    verifyNotNULL(ij_list)

    /*if all of-diagonal elements are zeros: */
    ij_list[0] = 0;
    ij_list[1] = 1;

    for (i=0 ; i<n ; i++){
        for (j = i+1 ; j<n ; j++){
            abs_curr = fabs(A[i][j]);
            if (abs_curr > max_val){
                ij_list[0] = i;
                ij_list[1] = j;
                max_val = abs_curr;
            }
        }
    }
    return ij_list;
}
/* ----------------------- calculate of f(A)^2 for Section 5 - Convergence: ---------------------------------*/

double calc_of_f_square(double **A, int n) {
    double diagonal_elementwise_sum_square = 0, matrix_elementwise_sum_square = 0, f_norm;
    int i, j;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            matrix_elementwise_sum_square += pow(A[i][j], 2);

    for (i = 0; i < n; i++)
        diagonal_elementwise_sum_square += pow(A[i][i], 2);

    f_norm = sqrt(matrix_elementwise_sum_square);

    return f_norm - diagonal_elementwise_sum_square;
}

/* ----------------------- calc_rotation_mat: ---------------------------------*/

rotation_mat *calc_rotation_mat(rotation_mat *P, double **A, int n) {
    /* returns a struct that represent the rotation matrix of A. */
    int sign_theta, i, j;
    double Aii, Ajj, Aij, theta, t, max_val = 0, abs_curr;

/* TODO, tests it
    int *ij = find_off_diag_max_abs_val(A, n);
    P->i = ij[0];
    P->j = ij[1];
*/
/*if all of-diagonal elements are zeros: */
    P->i = 0;
    P->j = 1;

    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            abs_curr = fabs(A[i][j]);
            if (abs_curr > max_val) {
                P->i = i;
                P->j = j;
                max_val = abs_curr;
            }
        }
    }

    Aii = A[P->i][P->i];
    Ajj = A[P->j][P->j];
    Aij = A[P->i][P->j];

    theta = (Ajj - Aii) / (2 * Aij);
    sign_theta = (theta == 0) ? 1 : (int) (theta / fabs(theta));
    t = sign_theta / (fabs(theta) + sqrt((theta * theta) + 1));

    P->c = 1 / (sqrt(t * t + 1));
    P->s = t * (P->c);
/* todo
    free(ij);
*/
    return P;
}

/*/////////////////////////////////// after Jacobi: /////////////////////////////////////////////////*/

typedef struct eigenStruct{
    double *Pointer;
    int originalIndex;  /* for a stable sort*/
} eigenStruct;

/* ------------------------ comparator of pointers:: ----------------------------------------*/

int compare_eigenStruct(const void *a, const void *b){
    /* compare function for a stable sort based on the eigenvalue */
    int res;

    eigenStruct A = *(eigenStruct*)a;
    eigenStruct B = *(eigenStruct*)b;


    if (*(A.Pointer) < *(B.Pointer))
        res = -1;
    else if (*(A.Pointer) > *(B.Pointer))
        res = 1;
    else{
        if ((A.originalIndex) < (B.originalIndex)) res = -1;
        else if ((A.originalIndex) > (B.originalIndex)) res = 1;
        else res = 0;
    }
    return res;
}

/* ------------------------ calc_T: ----------------------------------*/

double **calc_T(Jacobi_output *jacobiOutput, int n, int *k_pointer){
    /* if k == 0 it recalculate k and change it in memory. */
    int i,j, original_j, k;
    double line_sum, gap, max_gap = 0 ;
    double *eigenValues = jacobiOutput -> eigenValues;
    double **V = jacobiOutput -> V;
    double **U;
    double **T;
    eigenStruct eigenStructI;
    /* creating a list of pointers to the eigenvalues, that will be used to sort the eigenvectors according to the eigenvalue and original index: */
    eigenStruct *pointers_list = (eigenStruct*) malloc(n*sizeof(eigenStruct));
    verifyNotNULL(pointers_list)

    /* filling a pointers list of the eigenvalues: */
    for (i=0 ; i<n ; i++){
        eigenStructI.Pointer = &(eigenValues[i]);
        eigenStructI.originalIndex = i; /* saving the original index in eigenvalue array to enable a stable sort */
        pointers_list[i] = eigenStructI;
    }

    /* sorting the pointers list: */
    qsort(pointers_list, n, sizeof(eigenStruct), compare_eigenStruct);



    /* finding k using the sorted pointers list: */
    if (*k_pointer == 0){
        k = 1;
        for (i=1 ; i<n/2+1 ; i++){
            gap = fabs((*(pointers_list[i-1].Pointer)) - (*(pointers_list[i].Pointer)));
            if (gap > max_gap){
                max_gap=gap;
                k=i;
            }
        }
        *k_pointer = k;
    /* or extracting k that was given as a parameter: */
    } else{
        k = *k_pointer;
    }



    /* forming the matrix U that contain the k first eigenvectors: */
    U = alloc_nXm_matrix(n,k);
    for(j=0; j<k; j++){
        original_j = (int)(((pointers_list[j].Pointer) - eigenValues)); /* the column index of the eigenvector that fits eigenvalue number j. */
        for (i=0; i<n; i++){
            U[i][j] = V[i][original_j];
        }
    }



    /* forming the matrix T from U: */
    T = alloc_nXm_matrix(n,k);
    for(i=0; i<n; i++){
        line_sum = 0;
        for (j=0; j<k; j++){
            line_sum += U[i][j] * U[i][j];
        }
        for (j=0; j<k; j++){
            if (line_sum != 0)
                T[i][j] = U[i][j]/(sqrt(line_sum));
            else /* according to the answers at the forum it doesnt suppose to happen in a valid input data */
                T[i][j] = U[i][j];
        }
    }

    free(pointers_list);
    free_contiguous_mat(U);

    return T;
}

/*//////////////////////////////////////// kmeans functions: //////////////////////////////////////////////////////*/

/*-------------------------------- distance: ------------------------------------------------------*/

double calc_distance(const double *dot, const double *centroid, int d){
    double sum = 0;
    int i;
    for(i=0 ; i<d ; i++)
        sum += ((dot[i]-centroid[i])*(dot[i]-centroid[i]));
    return sum;
}

/*--------------------------------- print_list: -------------------------------------------  Temp!  -----*/

void print_list(double *array, int len){
    int i;
    double num;
    for(i=0; i < len; i++){
        num = array[i];
        if (num > -0.00005 && num < 0)
            printf("%.4f", 0.0);
        else
            printf("%.4f", array[i]);
        if(i!= len-1)
            printf(",");
    }
    printf("\n");
}

/*------------------------------ print_2d_array: ------------------------------------------------*/

void print_2d_array(double **array, int row_num, int col_num){
    int i, j;
    double num;
    for(i=0; i < row_num; i++)
    {
        for(j = 0; j < col_num; j++) {
            num = array[i][j];
            if (num > -0.00005 && num < 0)
                printf("%.4f", 0.0);
            else
                printf("%.4f", array[i][j]);
            if(j!= col_num-1)
                printf(",");
        }
        printf("\n");
    }
}
void print_2d_array_transpose(double **array, int row_num, int col_num){
    int i, j;
    double num;
    for(i=0; i < col_num; i++)
    {
        for(j = 0; j < row_num; j++) {
            num = array[j][i];
            if (num > -0.00005 && num < 0)
                printf("%.4f", 0.0);
            else
                printf("%.4f", array[j][i]);
            if(j!= row_num-1)
                printf(",");
        }
        printf("\n");
    }
}

/* ----------------------------- deepCopy2DArray: --------------------------------------------*/

double **deepCopy2DArray(double **A, int row_num, int col_num){
    int i, j;
    double **B = (double**) alloc_nXm_matrix(row_num,col_num);
    for (i = 0; i < row_num; i++)
        for (j = 0; j < col_num; j++)
            B[i][j] = A[i][j];

    return B;
}

/*------------------------------- find_nearest_centroid: -----------------------------------------*/
/* returns the index of the nearest_centroid of the dot*/
int find_nearest_centroid(double *dot, double **centroids, int k, int d){
    double min_distance = calc_distance(dot, centroids[0], d);
    double distance;
    int min_index = 0, i;
    for (i = 1 ; i < k ; ++i){
        distance = calc_distance(dot, *(centroids+i), d);
        if (distance < min_distance){
            min_distance = distance;
            min_index = i;
        }
    }
    return min_index;
}

/*------------------------------ check_equals_2d_list: ------------------------------------------ */

/*compare between two 2 dimensions lists of the same dimensions:*/
int check_equals_2d_list(double **list1, double **list2, int row_num, int col_num){
    int i, j;
    for(i = 0; i < row_num; i++) {
        for (j = 0; j < col_num; j++) {
            if (list1[i][j] != list2[i][j])
                return 0;
        }
    }
    return 1;
}
/*------------------------------- update_centroid: -----------------------------------------*/

void update_centroids(double **dots,const int *dots_location, double **centroids, int n, int d, int k) {

    /* space allocation */
    int i, j;
    int *num_of_dots_per_cluster = (int*) malloc(k * sizeof(int));
    double **sum_of_dots_per_cluster = (double **) malloc(k*sizeof(double *));
    verifyNotNULL(num_of_dots_per_cluster)
    verifyNotNULL(sum_of_dots_per_cluster)


    for (i = 0; i < k; i++) {
        sum_of_dots_per_cluster[i] = (double *) calloc(d, sizeof(double));
        verifyNotNULL(sum_of_dots_per_cluster)
    }

    /* initializing num_of_dots_per_cluster */
    for (i = 0; i < k; i++) {
        num_of_dots_per_cluster[i] = 0;
    }

    /* assignments in the sum and num arrays */
    for (i=0; i<n; i++){
        num_of_dots_per_cluster[dots_location[i]] += 1;
        for(j=0; j<d; j++)
            sum_of_dots_per_cluster[dots_location[i]][j] += dots[i][j];
    }

    /* assignments in the centroids array */
    for (i=0; i<k; i++){
        for(j=0; j<d; j++)
            centroids[i][j] = (sum_of_dots_per_cluster[i][j])/num_of_dots_per_cluster[i];
    }
    /* free space */
    free(num_of_dots_per_cluster);
    for(i = 0; i<k; i++)
        free(sum_of_dots_per_cluster[i]);
    free(sum_of_dots_per_cluster);
}

/*---------------------------------- kmeans Algorithm: ---------------------------------------------*/

double **kmeans(double **dots,double **centroids, int k, int d, int n, int max_iter) {

    /* variables declaration: */
    int i, j, iter_counter = 0;
    int *dots_location;
    double **old_centroids;

    /* Initialize a centroids matrix, and a dots cluster location array :*/
    dots_location = (int*) calloc(n, sizeof(int));
    old_centroids = (double**) malloc(k*sizeof(double*));
    verifyNotNULL(dots_location)
    verifyNotNULL(old_centroids)

    for(i=0; i<k; i++) {
        old_centroids[i] = (double *) malloc(d * sizeof(double));
        verifyNotNULL(old_centroids[i])
    }

    for(i=0; i<n; i++)
        dots_location[i] = -1;


    /* the kmeans algorithm flow :*/
    while(iter_counter < max_iter){

        for(i=0 ; i<n ; i++) {
            dots_location[i] = find_nearest_centroid(dots[i], centroids, k, d);
        }

        for(i=0 ; i<k ; i++)
            for(j=0 ; j<d ; j++)
                old_centroids[i][j] = centroids[i][j];

        update_centroids(dots, dots_location, centroids, n, d, k);

        if(check_equals_2d_list(old_centroids, centroids, k, d))
            break;

        iter_counter+=1;
    } /* end of while */

    /* frees */
    free(dots_location);
    for(i=0 ; i<k ; i++) {
        free(old_centroids[i]);
    }
    free(old_centroids);

    return centroids;
}