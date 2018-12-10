#ifndef MATRIX_TEST_IS_INCLUDED
#define MATRIX_TEST_IS_INCLUDED
  
#include <blitz/array.h>

#include <iostream>

using namespace  blitz;

    // ******   data declaration  ******* 

typedef   struct {
  int
           matrType,
           dataType,
            matrRow,
            matrCol;
} INPUTDATA;

    // *****  function declaration *******

void inputData(INPUTDATA &data);

  /*
  ** The template function 
  **      realMatrixTest()
  ** creates three matrices Array<typename T,2>
  ** A, B, C and fills all elements with values
  ** according to the input data and calulateS
  ** the sum C = A + B and product C = A x B and 
  ** writes the result to standard output
  */

template <typename T>
void realMatrixTest(INPUTDATA input)
{
  int    idum = 10;

  Array<T,2> A(input.matrRow, input.matrCol), 
             B(input.matrRow, input.matrCol), 
             C(input.matrRow, input.matrCol); 

  switch(input.matrType) {
    case 1: randomRealVal(A, B, input);
            break;
    case 2: readRealVal(A, B, input);
            break;
  } // end switch

  C = A + B;        // calculate the sum of matrix C

  cout << endl << "A = " << A  << endl;
  cout << endl << "B = " << B  << endl;
  cout << endl << "Test sum matrix C = A + B";
  cout << endl << "C = " << C  << endl;

  MatrixMult(A,B,C);   // calculate the matrix product C = A x B

  cout << endl << "Test matrix product C = A x B";
  cout << endl << "C = " << C  << endl;

} // End: template function realMatrixTest()

  /*
  ** The template function 
  **      coplexMatrixTest()
  ** creates three complex matrices Array<complex<typename T>, 2 >
  ** A, B, C and fills all elements with values
  ** according to the input data and calulateS
  ** the sum C = A + B and product C = A x B and 
  ** writes the result to standard output
  */

template <typename T>
void complexMatrixTest(INPUTDATA input)
{
  int    idum = 10;

  Array<complex<T>,2> A(input.matrRow, input.matrCol), 
                      B(input.matrRow, input.matrCol), 
                      C(input.matrRow, input.matrCol); 

  switch(input.matrType) {
    case 3: randomComplexVal(A, B, input);
            break;
    case 4: readComplexVal(A, B, input);
            break;
  } // end switch

  C = A + B;        // calculate the sum of matrix C

  cout << endl << "A = " << A  << endl;
  cout << endl << "B = " << B  << endl;
  cout << endl << "C = " << C  << endl;

  MatrixMult(A,B,C);   // calculate the matrix product C = A x B

  cout << endl << "C = " << C  << endl;

} // End: template function complexMatrixTest()

  /*
  ** The template function 
  **      randomRealVal()
  ** fills all elements in matrices
  ** A and B with random real numbers
  */

template <typename T>
void randomRealVal(Array<T,2>& A, Array<T,2>& B, INPUTDATA input)
{
  int    idum = 10;

  for(int i = 0; i < input.matrRow ; i++) {
    for(int j = 0; j < input.matrCol ; j++) {
      A(i,j) = static_cast<T>(10.0 * (ran0(idum) - 0.5));
    }
  } 
  for(int i = 0; i < input.matrRow ; i++) {
    for(int j = 0; j < input.matrCol ; j++) {
      B(i,j) = static_cast<T>(10.0 * (ran0(idum) - 0.5));
    }
  }

} // End: template function randomrealVal()

  /*
  ** The template function 
  **      readRealVal()
  ** fills all elements in matrices A 
  ** and B with data from standard input
  */

template <typename T>
void readRealVal(Array<T,2>& A, Array<T,2>& B, INPUTDATA input)
{
  cout << endl << endl 
       << "type in datafor matrix A(" 
       << input.matrRow << "," << input.matrCol
       << " ): " << endl;

  cin >> A;
  cout << endl << endl 
       << "type in date for matrix B(" 
       << input.matrRow << "," << input.matrCol
       << " ): " << endl;

  cin >> B;

} // End: template function readRealVal()

  /*
  ** The template function 
  **      randomComplexVal()
  ** fills all elements in matrices
  ** A and B with random complex numbers
  */

template <typename T>
void randomComplexVal(Array<complex<T>,2>& A, Array<complex<T>,2>& B,
                                                       INPUTDATA input)
{
  int    idum = 10;
  T realComp, imagComp;

  for(int i = 0; i < input.matrRow ; i++) {
    for(int j = 0; j < input.matrCol ; j++) {
      realComp = static_cast<T>(10.0 * (ran0(idum) - 0.5));
      imagComp = static_cast<T>(10.0 * (ran0(idum) - 0.5));
      cout << endl << "real A(" << i <<"," << j <<") = "<< realComp;
      cout << endl << "imag A(" << i <<"," << j <<") = "<< imagComp;
      A(i,j)   = complex<T>(realComp, imagComp); 
    }
  }
  for(int i = 0; i < input.matrRow ; i++) {
    for(int j = 0; j < input.matrCol ; j++) {
      realComp = static_cast<T>(10.0 * (ran0(idum) - 0.5));
      imagComp = static_cast<T>(10.0 * (ran0(idum) - 0.5));
      B(i,j)   = complex<T>(realComp, imagComp); 
    }
  }
} // End: template function randomComplexVal()

  /*
  ** The template function 
  **      readComplexVal()
  ** fills all elements in matrices A 
  ** and B with data from standard input
  */

template <typename T>
void readComplexVal(Array<complex<T>,2>& A, Array<complex<T>,2>& B,
                                                     INPUTDATA input)
{
  cout << endl << endl 
       << "type in date for complex matrix A(" 
       << input.matrRow << "," << input.matrCol
       << " ): " << endl;

  cin >> A;
  cout << endl << endl 
       << "type in date for complex matrix B(" 
       << input.matrRow << "," << input.matrCol
       << " ): " << endl;

  cin >> B;

} // End: template function readVal()

     /*
     ** The function 
     **     ran0()
     ** "Minimal" random number generator of Park and Miller. Returns a uniform
     ** random deviate between 0.0 and 1.0. Set or reset idum to any integer value
     ** (except the unlikely value MASK) to initialize the sequence; idum must not
     ** be altered between calls for successive deviates in a sequence. 
     */

double ran0(int &idum)
{
	const int IA=16807,IM=2147483647,IQ=127773;
	const int IR=2836,MASK=123459876;
	const double AM=1.0/(double)IM;
	int k;
	double ans;

	idum ^= MASK;
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	ans=AM*idum;
	idum ^= MASK;
	return ans;
} // End: function ran0()

   /*
   **  The template function 
   **         MatrixMult()
   ** calculates the matrix product between matrices A and B.
   ** As parameters it takes first matrix A and second matrix B,
   ** calculates and return the matrix product:
   **       C = A x B
   */  
  
template<class T> 
void MatrixMult(Array<T,2>& A, Array<T,2>& B, Array<T,2>& C)
{
  if(  (A.rows() != C.rows()) || (B.columns() != C.columns())
	 ||(A.columns() != B.rows() )) {
    cout << endl << "Error in function Matrix_Mult() :";
    cout << endl << "matrix A(" << A.rows() <<"."<< A.columns() << "), "
	 << "matrix B(" << B.rows() <<"."<< B.columns() << "), "
	 << "matrix C(" << C.rows() <<"."<< C.columns() << "), ";
    exit(1);
  }

  for(int row = 0; row < C.rows(); row++) {
    for(int col = 0; col < C.columns(); col++) {
      T temp = 0;
      for(int k = 0; k < A.columns(); k++) {
	temp += A(row,k) * B(k,col);
      }
      C(row,col) = temp;
    } // end C.columns index
  } // end C.rows index

} // End: template function MatrixMult() 


#endif
