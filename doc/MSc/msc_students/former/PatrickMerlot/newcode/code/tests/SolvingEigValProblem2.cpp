/*
  Solves an eigenvalue problem AX=eX for a given Matrix A,
  in output: eigenvalues e and eigenvectors X
*/
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
using namespace  std;
// output file as global variable
ofstream ofile;  

// function declarations 

void initialise(double&, double&, int&, int&) ;
double potential(double);
int comp(const double *, const double *);
void output(double, double, int, double *);

int main(int argc, char* argv[])
{



  char *outfilename;
  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] << 
      " read also output file on same line" << endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }
  ofile.open(outfilename); 


//   //   Read in data 
    double** A;
   int dimA;
   dimA =3;
//   A = new double * [dimA];
//   for(int i=0; i<dimA; i++)
//     {
//       A[i]=new double [dimA];
//     }

//   double *dd, *ee;
//   dd = new double[dimA];
//   ee = new double[dimA];
  
//   A[0][0]=0;
//   A[0][1]=1;
//   A[0][2]=-1;
//   A[1][1]=1;
//   A[1][2]=0;
//   A[2][2]=1;
//   A[1][0]= A[0][1];
//   A[2][0]= A[0][2];
//   A[2][1]= A[1][2];

//   printf("A=\n");
//   for(int i=0; i<dimA; i++)
//     for(int j=0; j<dimA; j++)
//       {
// 	printf("%2.0lf ",A[i][j]);
// 	if(j==dimA-1)
// 	  printf("\n");
//       }

//   tred2(A, dimA, dd, ee);
//    printf("A=\n");
//   for(int i=0; i<dimA; i++)
//     for(int j=0; j<dimA; j++)
//       {
// 	printf("%2.0lf ",A[i][j]);
// 	if(j==dimA-1)
// 	  printf("\n");
//       }

//   for(int i = 0; i < dimA; i++)    
//     for(int j = 0; j < dimA; j++) 
//       {
// 	A[i][j] = 0.0;
// 	A[i][i] = 1.0;
//       }




//   for(int i=0; i<dimA; i++)
//       printf("d[%d]=%lf\n",i,dd[i]);
//   for(int i=0; i<dimA; i++)
//       printf("e[%d]=%lf\n",i,ee[i]);

//   // diagonalize and obtain eigenvalues
//   tqli(dd, ee, dimA, A);      
  
//   for(int i=0; i<dimA; i++)
//     printf("eigVal[%d]=%lf\n",i,dd[i]);

// printf("Z=\n");
//   for(int i=0; i<dimA; i++)
//     for(int j=0; j<dimA; j++)
//       {
// 	printf("%2.0lf ",A[i][j]);
// 	if(j==dimA-1)
// 	  printf("\n");
//       }

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  double       **a, **z, sum, *ddd, *eee;
  
 ddd= new double [dimA]; 
  eee = new double [dimA];
  z = (double **) matrix(dimA, dimA, sizeof(double));   
  a = (double **) matrix(dimA, dimA, sizeof(double));   
  a[0][0] = 0.0;
  a[0][1] = 1.0;
  a[0][2] = -1.0;
  a[1][0] = 1.0;
  a[1][1] = 1.0;
  a[1][2] = 0.0;
  a[2][0] = -1.0;
  a[2][1] = 0.0;
  a[2][2] = 1.0;


printf("A=\n");
  for(int i=0; i<dimA; i++)
    for(int j=0; j<dimA; j++)
      {
	printf("%2.0lf ",a[i][j]);
	if(j==dimA-1)
	  printf("\n");
      }

// set eigenvector matrix
  for(int i = 0; i < dimA; i++)    
    for(int j = 0; j < dimA; j++) 
      {
	z[i][j] = 0.0;
	z[i][i] = 1.0;
      }

  printf("\n Z=\n");
  for(int i=0; i<dimA; i++)
    for(int j=0; j<dimA; j++)
      {
	printf("%12.4E ",z[i][j]);
	if(j==dimA-1)
	  printf("\n");
      }
 


  tred2(a, dimA, ddd, eee);
 for(int i = 0; i < dimA; i++)    
    for(int j = 0; j < dimA; j++) 
	z[i][j] = a[i][j];


  tqli(ddd,eee,dimA, z);

  for(int i = 0; i < dimA; i++) {
        printf("  d[%2d] = %12.4E\n",i, ddd[i]);
  }
  
 printf("\nZ=\n");
  for(int i=0; i<dimA; i++)
    for(int j=0; j<dimA; j++)
      {
	printf("%12.4E ",z[i][j]);
	if(j==dimA-1)
	  printf("\n");
      }

  // compute norm of each eigenvector
  printf("\nnorm of eigenvectors\n");
  double sum2=0.0;
  for(int j=0;j<dimA;j++)
    {
      sum2=0.0;
    for(int i=0;i<dimA;i++)      
      sum2 += pow(z[i][j],2); 
    printf("norm of z(%12.4E)=%12.4E\n",ddd[j],sum2);
    }
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



//   // Sort eigenvalues as an ascending series 
//   qsort(d,(UL) max_step - 1,sizeof(double),
//          (int(*)(const void *,const void *))comp);
//   // send results to ouput file
//   output(r_min , r_max, max_step, d);
//   delete [] r; delete [] w; delete [] e; delete [] d; 
//   free_matrix((void **) z); // free memory
  ofile.close();  // close output file
  return 0;
} // End: function main() 

/*
  The function potential()
  calculates and return the value of the 
  potential for a given argument x.
  The potential here is for the hydrogen atom
*/        

double potential(double x)
{
   return -2./x;

} // End: function potential()  

/*
  The function   int comp()                  
  is a utility function for the library function qsort()
  to sort double numbers after increasing values.
*/       

int comp(const double *val_1, const double *val_2)
{
  if((*val_1) <= (*val_2))       return -1;
  else  if((*val_1) > (*val_2))  return +1;
  else                     return  0; 
} // End: function comp() 


void initialise(double& r_min, double& r_max, int& orb_l, int& max_step) 
{
  // cout << "Min vakues of R = ";
//   cin >> r_min;
//   cout << "Max value of R = ";
//   cin >> r_max;
//   cout << "Orbital momentum = ";
//   cin >> orb_l;
//   cout << "Number of steps = ";
//   cin >> max_step;
}  // end of function initialise   




void output(double r_min , double r_max, int max_step, double *d)
{
  int i;
  ofile << "RESULTS:" << endl;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile <<"R_min = " << setw(15) << setprecision(8) << r_min << endl;  
  ofile <<"R_max = " << setw(15) << setprecision(8) << r_max << endl;  
  ofile <<"Number of steps = " << setw(15) << max_step << endl;  
  ofile << "Five lowest eigenvalues:" << endl;
  for(i = 0; i < 5; i++) {
    ofile << setw(15) << setprecision(8) << d[i] << endl;
  }
}  // end of function output         

