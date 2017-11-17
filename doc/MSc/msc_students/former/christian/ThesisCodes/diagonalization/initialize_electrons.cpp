#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;



void InitializeTwoElectrons(int N, mat &A, double xMin, double xMax, mat &SaveEigenvector, double omega_r)
{
    double h = (xMax-xMin)/N;

    //Set up the vector x and the matrix A:
    vec x = xMin + linspace(0, N, N+1)*h;

    //No repulsion:
    //vec V = omega_r*omega_r*(x%x);

    //With Repulsion:
    vec V = omega_r*omega_r*(x%x) + 1.0/x;

    SaveEigenvector.col(0) = x.subvec(1, N-1);    //Saves the x vector for output.

    double Constant = 1/(h*h);
    A.diag(0)  =  2*Constant + V.subvec(1,N-1);     //Set d_i elements in A
    A.diag(1)  = -Constant*ones(N-2);               //Set e_i elements in A
    A.diag(-1) = A.diag(1);                         //Set e_i elements in A

    return;
}
