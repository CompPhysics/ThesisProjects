//     Simple test case of matrix operations
//     using Blitz++
#include <blitz/array.h> 
#include <iostream>
using namespace std;
using namespace blitz;

int main()
{
  // je joue avec Blitz
  Array<int,2>  PAT1;
  cout << "matrix PAT= " << PAT1 << endl;
  cout << "size PAT1= " << PAT1.size() << endl;

Array<double,1>  PAT0;
  cout << "matrix PAT0= " << PAT0 << endl;
  cout << "size PAT0= " << PAT0.size() << endl;

  Range rPat(0,5);
 cout << "matrix rPAT= " << rPat << endl;

  PAT0.resize(rPat);
 cout << "matrix PAT0= " << PAT0 << endl;
 cout << "size PAT0= " << PAT0.size() << endl;

 PAT0(3)=432;
 cout << "matrix PAT0= " << PAT0 << endl;
 cout << "size PAT0= " << PAT0.size() << endl;
  
double data[] = { 1, 2, 3, 4 };
Array<double,2> PAT2(data, shape(2,2), neverDeleteData);   // Make a 2x2 array
  cout << "size PAT2= " << PAT2.size() << endl;
  cout << "PAT2= " << PAT2 << endl;
Array<double,2> PAT3(data, shape(1,4), neverDeleteData);   
 cout << "PAT3= " << PAT3 << endl;
Array<double,2> PAT4(data, shape(4,1), neverDeleteData);   
  cout << "PAT4= " << PAT4 << endl;
  cout << "size PAT4= " << PAT4.size() << endl;

Array<double,2> PAT5(data, shape(4,4), neverDeleteData);   
 Range aa(0,3);
 Range bb(0,3);
 cout << "PAT3= " << PAT3(0,aa) << endl;
 PAT5 = PAT3(0,2)*PAT3(0,1);
 cout << "PAT3= " << PAT3 << endl;
  cout << "size PAT5=PAT3(i)*PAT3(j)= " << PAT5 << endl;

  PAT2.resizeAndPreserve(shape(2,4));
  cout << "PAT2= " << PAT2 << endl;
  cout << "PAT2(1,2)= " << PAT2(1,1) << endl;


  Array<double,1> PAT6(4),PAT7(4);
  Array<double,2> PAT8(4);
  PAT6= 1,2,3,4;
  PAT7=1,0,0,1;

  firstIndex i;    // Placeholder for the first index
  secondIndex j;   // Placeholder for the second index
  PAT8=PAT6(i)*PAT7(j);
  cout << "PAT6= " << PAT6 << endl;
  cout << "PAT7= " << PAT7 << endl;
  cout << "PAT8=PAT6*PAT7 " << PAT8 << endl;
  cout << "size PAT8= " << PAT8.size() << endl;


  // Create two 4x4 arrays.  We want them to look like matrices, so
  // we'll make the valid index range 1..4 (rather than 0..3 which is
  // the default).

  Range r(1,4);
  cout << "r= " << r << endl;

  Array<float,2> A(r,r), B(r,r);
  A=3;
  cout << "test A(1,1)=" << A(1,1) << endl;
  
  // The first will be a Hilbert matrix:
  //
  // a   =   1
  //  ij   -----
  //       i+j-1
  //
  // Blitz++ provides a set of types { firstIndex, secondIndex, ... }
  // which act as placeholders for indices.  These can be used directly
  // in expressions.  For example, we can fill out the A matrix like this:

//   firstIndex i;    // Placeholder for the first index
//   secondIndex j;   // Placeholder for the second index

  A = 1.0 / (i+j-1); 

  cout << "A = " << A << endl;
  // Now the A matrix has each element equal to a_ij = 1/(i+j-1).

  // The matrix B will be the permutation matrix
  //
  // [ 0 0 0 1 ]
  // [ 0 0 1 0 ]
  // [ 0 1 0 0 ]
  // [ 1 0 0 0 ]
  //
  // Here are two ways of filling out B:

  B = (i == (5-j));         // Using an equation -- a bit cryptic

  cout << "B = " << B << endl;

  B = 0, 0, 0, 1,           // Using an initializer list
    0, 0, 1, 0,           
    0, 1, 0, 0,
    1, 0, 0, 0;

  cout << "B = " << B << endl;

  // Now some examples of tensor-like notation.

  Array<float,3> C(r,r,r);  // A three-dimensional array: 1..4, 1..4, 1..4

  thirdIndex k;             // Placeholder for the third index

  // This expression will set
  //
  // c    = a   * b
  //  ijk    ik    kj

  C = A(i,k) * B(k,j);

  // In real tensor notation, the repeated k index would imply a
  // contraction (or summation) along k.  In Blitz++, you must explicitly
  // indicate contractions using the sum(expr, index) function:

  Array<float,2> D(r,r);

  D = sum(A(i,k) * B(k,j), k);

  // The above expression computes the matrix product of A and B.

  cout << "D = " << D << endl;

  // Now let's fill out a two-dimensional array with a radially symmetric
  // decaying sinusoid.

  int N = 64;                   // Size of array: N x N
  Array<float,2> F(N,N);
  float midpoint = (N-1)/2.;
  int cycles = 3;
  float omega = 2.0 * M_PI * cycles / double(N);
  float tau = - 10.0 / N;

  F = cos(omega * sqrt(pow2(i-midpoint) + pow2(j-midpoint)))
    * exp(tau * sqrt(pow2(i-midpoint) + pow2(j-midpoint)));



  Array<double,2> E (5,5);
  E=0.0;
  for(int i=0;i<5;i++)
    E(i,i)=5.0;
  cout << "E = i: " << E << endl;

  E/=6.0;
  cout << "E= : " << E << endl;

  return 0;
}

