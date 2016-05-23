#include "methods.hpp"
#include "lib.h"

/*
 * Function to delete matrices of doubles
 * with N rows (no. of columns not required)
 */
void delete_matrix(double** matrix, int N) {
  for (int i = 0; i < N; i++) {
    delete[] matrix[i];
  }
  delete[] matrix;
  matrix = NULL;
}
/*
 * Function to setup NxM matrices of doubles
 */
double** create_matrix(int N, int M) {
  double** matrix = new double*[N];
  for (int i = 0; i < N;i++) {
    matrix[i] = new double [M];
  }
  return matrix;
}

// random numbers with gaussian distribution
double gaussian_deviate(long * idum)
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if ( idum < 0) iset =0;
  if (iset == 0) {
    do {
      v1 = 2.*ran1(idum) -1.0;
      v2 = 2.*ran1(idum) -1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.);
    fac = sqrt(-2.*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset =0;
    return gset;
  }
}



/* The function
**                inverse()
** perform a mtx inversion of the input matrix a[][] with
** dimension n. The method is described in Numerical Recipes
** sect. 2.3, page 48.
*/

void inverse(double **a, int n)
{        
  int          i,j, *indx;
  double       d, *col, **y;

  // allocate space in memory
  indx = new int[n];
  col  = new double[n];
  y    = (double **) matrix(n, n, sizeof(double)); 
   
  ludcmp(a, n, indx, &d);   // LU decompose  a[][] 

 
  for(j = 0; j < n; j++) {

    // initialize right-side of linear equations 

    for(i = 0; i < n; i++) col[i] = 0.0;
    col[j] = 1.0;

    lubksb(a, n, indx, col);

    // save result in y[][] 

    for(i = 0; i < n; i++) y[i][j] = col[i];

  }   //j-loop over columns 
   
  // return the inverse matrix in a[][] 

  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) a[i][j] = y[i][j];
  } 
  free_matrix((void **) y);     // release local memory 
  delete [] col;
  delete []indx;

}  // End: function inverse()