double **kmeans(double **dots,double **centroids, int k, int d, int n, int max_iter);
void printByGoal(double **points, int d, int n,char *goal );
void print_2d_array(double **array, int row_num, int col_num);
void print_2d_array_transpose(double **array, int row_num, int col_num);
void free2DMalloc(double **points,int n);
double** getT(double **dots, int d, int n,int *k);
void free_contiguous_mat(double **mat);
