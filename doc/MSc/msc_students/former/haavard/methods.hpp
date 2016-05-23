/*
 * Function to delete matrices of doubles
 * with N rows (no. of columns not required)
 */
void delete_matrix(double** matrix, int N);
 
/*
 * Function to setup NxM matrices of doubles
 */
double** create_matrix(int N, int M);

// random numbers with gaussian distribution
double gaussian_deviate(long * idum);

//inverting a matrix
void inverse(double **a, int n);