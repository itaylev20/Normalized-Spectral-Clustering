#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spk_algorithm.h"

#define NJump 150
#define MAX_ITER 300
#define MAX_BUFFER 700

/* goals functions for main: */
double** getWeightedAdjacencyMatrix(double **dots, int d, int n);
double** getDiagonalDegreeMatrix(double **dots, int d, int n);
double** getNormalizedGraphLaplacian(double **dots, int d, int n);
Jacobi_output* getJacobi(double **dots, int d, int n);
void printWeightedAdjacencyMatrix(double **dots, int d, int n );
void printDiagonalDegreeMatrix(double **dots, int d, int n );
void printNormalizedGraphLaplacian(double **dots, int d, int n );
void printJacobi(double **dots, int n );
void printSpectralClustering(double **dots, int d, int n, int k, int max_iter );
void printByGoal(double **points, int d, int n,char *goal );
/* free functions: */
void free2DMalloc(double **points,int n);
/* functions for main: */
double ** scanInput(char* filename, int *n, int *d);
double** getT(double **dots, int d, int n,int *k);


/*------------------------ goals functions for main ----------------------------*/
double** getWeightedAdjacencyMatrix(double **dots, int d, int n){
    return calc_weighted_adjacency_matrix(dots,d,n);
}
double** getDiagonalDegreeMatrix(double **dots, int d, int n){
    double **W = getWeightedAdjacencyMatrix(dots,d,n);
    double **D = calc_diagonal_degree_matrix(W,n);
    free_contiguous_mat(W);
    return D;
}
double** getNormalizedGraphLaplacian(double **dots, int d, int n){
    double **W = getWeightedAdjacencyMatrix(dots,d,n);
    double **D = getDiagonalDegreeMatrix(dots,d,n);
    double **L_norm = calc_L_norm(D,W,n);
    free_contiguous_mat(W);
    free_contiguous_mat(D);
    return L_norm;
}

Jacobi_output *getJacobi(double **dots, int d, int n) {
    /* returns the output of jacobi operated on L_norm, for the SPk algorithm */
    double **L_norm = getNormalizedGraphLaplacian(dots, d, n);
    Jacobi_output *jacobiOutput = jacobi(L_norm, n);
    free_contiguous_mat(L_norm);
    return jacobiOutput;
}
double** getT(double **dots, int d, int n, int *k){
    Jacobi_output *jacobi_output = getJacobi(dots,d,n);
    double **T = calc_T(jacobi_output,n,k);
    free_contiguous_mat(jacobi_output->V);
    free(jacobi_output->eigenValues);
    free(jacobi_output);
    return T;
}
void printWeightedAdjacencyMatrix(double **dots, int d, int n ){
    double **W = getWeightedAdjacencyMatrix(dots, d, n);
    print_2d_array(W, n, n);
    free_contiguous_mat(W);
}
void printDiagonalDegreeMatrix(double **dots, int d, int n ){
    double **D = getDiagonalDegreeMatrix(dots, d, n);
    print_2d_array(D, n, n);
    free_contiguous_mat(D);
}
void printNormalizedGraphLaplacian(double **dots, int d, int n ){
    double **L_norm = getNormalizedGraphLaplacian(dots, d, n);
    print_2d_array(L_norm, n, n);
    free_contiguous_mat(L_norm);
}
void printJacobi(double **dots, int n ){
    Jacobi_output *jacobiOutput = jacobi(dots, n);
    print_list(jacobiOutput -> eigenValues, n);
    print_2d_array_transpose(jacobiOutput -> V, n, n);

    /*frees:*/
    free_contiguous_mat(jacobiOutput->V);
    free(jacobiOutput->eigenValues);
    free(jacobiOutput);
}
void printSpectralClustering(double **dots, int d, int n, int k, int max_iter ){
    double **T,**final_centroids, **initialCentroids;
    T = getT(dots,d,n,&k);


    initialCentroids = deepCopy2DArray(T,k,k); /* defining the first k 'dots' (from T) as the initial  centroids. Like in Ex1. */

    final_centroids = kmeans(T, initialCentroids, k, k, n, max_iter);
    print_2d_array(final_centroids,k,k);

    /* frees: */
    free_contiguous_mat(T);
    free_contiguous_mat(initialCentroids);
    /* no need to free final_centroids because it points to the same address as initialCentroids */

}

void printByGoal(double **points, int d, int n, char *goal ){ /* Print by Goal except for Spk (spk uses the function above) */
    if (strcmp(goal,"wam") == 0){
        printWeightedAdjacencyMatrix(points,d,n);
    }
    else  if (strcmp(goal,"ddg") == 0){
        printDiagonalDegreeMatrix(points,d,n);
    }
    else  if (strcmp(goal,"lnorm") == 0){
        printNormalizedGraphLaplacian(points,d,n);

    }
    else  if (strcmp(goal,"jacobi") == 0){
        printJacobi(points,n);
    }
    else { /* Not suppose to happen - just in case. */
        printf("Invalid input.");
    }
}

/*////////////////////////////////////////////////////// frees: /////////////////////////////////////////////////////*/

void free2DMalloc(double **points,int n) {
    int i;
    for(i=0; i<n; i++) {
        free(points[i]);
    }
    free(points);
}

/*------------------------------------- Scan Matrix-----------------------------------------------*/

double ** scanInput(char* filename, int *n, int *d) {
    double **vectors;
    double *currentVector;
    char *num;
    int i;
    char *line;

    FILE *file = fopen(filename, "r");
    verifyNotNULL(file)

    line = malloc(MAX_BUFFER*sizeof(char));
    verifyNotNULL(line)


    /* Scan to get d: */
    if(fgets(line, MAX_BUFFER, file) && *line != EOF) {
        strtok(line, ",");
        *d = 1;
        while (strtok(NULL, ","))
            *d += 1;
    }
    else
        printf("An Error Has Occured");

    /* Now we will scan all the input, line by line, each time we will allocate a vector of D size to assign the values.
     We will define **double vectors as first NJump size, and if necessary we will increase the size.
     In the necessary, we will resize the vectors to the actual n size. */

    vectors = (double**)malloc(NJump * sizeof(double *));
    verifyNotNULL(vectors)

    rewind(file);
    *n=0;
    while(fgets(line, MAX_BUFFER, file) && *line != EOF){

        if((*n)!=0 && (*n) % NJump == 0){
            vectors = (double **)realloc(vectors, ((((*n) / NJump) + 1) * NJump) * sizeof(double *));
            verifyNotNULL(vectors)
        }

        currentVector = (double *)malloc((*d) * sizeof(double));
        verifyNotNULL(currentVector)

        num = strtok(line, ",");
        i = 0;
        while (num != NULL){
            currentVector[i]= atof(num);
            num = strtok(NULL, ",");
            i++;
        }

        vectors[*n] = currentVector;
        (*n) += 1;
    }


    if (*n % NJump != 0) {
        vectors = (double **) realloc(vectors, (*n) * sizeof(double *));
        verifyNotNULL(vectors)
    }

    fclose(file);
    /*frees: */
    free(line);

    return vectors;
}

/*//////////////////////////////////////// main: ///////////////////////////////////////////////////////////*/

int main(int argc, char* argv[]) {
    int k,d=-1,n=-1;
    char *goal,*filename;
    double **points;

    /* As described in the instructions - we can assume input is always valid */
    assert(argc == 4 ); /* TODO, Make more informative - add error massage */
    k = atoi(argv[1]);
    goal = argv[2];
    filename = argv[3];

    points = scanInput(filename, &n, &d);
    verifyNotNULL(points)

    if (! (k >= 0 && k<n)){  /* TODO Discussion ID 159904, -  No need to check K if goal!=spk. */
        printf("Invalid Input!");
        assert(-1);
    }

    if (strcmp(goal,"spk") == 0){
        printSpectralClustering(points, d, n, k, MAX_ITER);
    }
    else {
        printByGoal(points,d,n,goal);
    }

    free2DMalloc(points,n);

    return 0;
}
