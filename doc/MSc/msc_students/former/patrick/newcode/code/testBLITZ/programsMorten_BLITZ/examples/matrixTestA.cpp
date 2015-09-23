#include <blitz/array.h>


#include <iostream>
#include "matrixTestA.hpp"

using namespace  blitz;

    /*
    ** The program 
    **           matrix_test()
    ** calculates matrix sum and product between matrixes A and B.
    ** The sum 
    **             C = A + B
    ** since the + operator is overloaded in BLITZ
    ** The matrix product 
    **             C = A x B
    ** in implemented through a function MatrixMult()
    */  

int main()
{
  INPUTDATA   data;

  while(true) {
    cout << endl << endl
         << "  Program to test matrix operation in BLITZ";

    inputData(data);       // data from standard input
  
    if(data.matrType < 3) {
      switch(data.dataType) {
        case 1: realMatrixTest<int>(data);
                break;  
        case 2: realMatrixTest<double>(data);
                break;  
        case 3: realMatrixTest<float>(data);
                break;
     }  // end switch()
      }
      else {
	switch(data.dataType) {
          case 1: complexMatrixTest<int>(data);
                  break;  
          case 2: complexMatrixTest<double>(data);
                  break;  
          case 3: complexMatrixTest<float>(data);
                  break;
	} // end switch()
      }
  } // end while() loop


  return 0;

} // End: program main()


void inputData(INPUTDATA& input)
{
  while(true)  {
    cout << endl << endl    
	 << "Type of matrices :"
	 << endl << endl;

    cout << "Type 0 - break program"            << endl
         << "Type 1 - Random real matrices "    << endl
         << "Type 2 - Type in real matrices"    << endl
         << "Type 3 - Random complex matrices " << endl
         << "Type 4 - Type in complex matrices" << endl 
         << endl << "   Choice = ";

    cin >> input.matrType;

    if(input.matrType == 0) exit(0);    // terminate the program

    if((input.matrType < 1) || (input.matrType > 4)) {
      cout << endl << endl <<"Wrong matrix choice - try once more";
      continue;
    }

    cout << "Data  type in the matrices:" << endl
         << " Type 1 - int     matrix"    << endl
         << " Type 2 - double  matrix"    << endl
         << " Type 3 - float   matrix"    << endl
         << endl << "   Choice = ";

    cin >> input.dataType;
 
    if((input.dataType < 1) || (input.dataType > 3)) {
      cout << endl << endl <<"Wrong data choice - try once more";
      continue;
    }
    cout << "Matrix dimensions" << endl
         << endl <<"    rows =";
    cin  >> input.matrRow;
    cout << endl <<"    cols =";
    cin  >> input.matrCol;

    if((input.matrRow <= 0) || (input.matrCol <= 0)) {
      cout << endl << "   Wrong matrix dimension - try once more !";
      continue;
    }
    break;       // all necessary data are receivec
  } // end while() loop

} // End: function inputData()
