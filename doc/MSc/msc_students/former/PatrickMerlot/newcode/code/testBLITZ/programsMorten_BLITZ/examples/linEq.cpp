    /*
    ** The program 
    **           matrix_test()
    ** solves the set of  n linear eguation
    **            A x X = B
    ** The Blitz matrices Array<T,2> A(n,n) and 
    ** Array<T.1> B(n) are specified through input
    ** data.
    ** The result is stored in Array<T,1> X(n).
    ** The program checks the solution
    */

#include "lib.hpp"
#include "linEq.hpp"

        // *** local function ***

void inputData(INPUTDATA& input);

int main()
{
  INPUTDATA   input;

     /*
     ** all data are read from standard input and
     ** returned to main() by the reference pointer
     */

  inputData(input);

  switch(input.dataType) {
    case 1:  solveEq<double>(input);
             break;
    case 2:  solveEq<float>(input);	     
             break;
  } // end switch()

  return 0;

} // End: function main()

    /*
    ** The function 
    **      inputData()
    ** reads all necessary input data and store them in 
    ** struct INPUTDATA input. The data is returned to
    ** function main() using reference pointer
    */

void inputData(INPUTDATA& input)
{
  while(1)  {
    cout << endl << endl    
	 << "Type of matrices :"
	 << endl << endl;

    cout << "Type 0 - break program"            << endl
         << "Type 1 - Random matrices A and B"  << endl
         << "Type 2 - Type in matrices A and B" << endl
         << endl << "   Choice = ";

    cin >> input.matrType;

    if(input.matrType == 0) exit(0);    // terminate the program

    if((input.matrType < 1) || (input.matrType > 2)) {
      cout << endl << endl <<"Wrong matrix choice - try once more";
      continue;
    }
    cout << "Data  type in the matrices?" << endl
         << " Type 1 - double  matrices" << endl
         << " Type 2 - float   matrices" << endl
         << " Type 3 - complex matrices" << endl
         << endl << "   Choice = ";

    cin >> input.dataType;
 
    if((input.dataType < 1) || (input.dataType > 3)) {
      cout << endl << endl <<"Wrong data choice - try once more";
      continue;
    }
    cout << "Matrix dimension n = ";
    cin  >> input.dim;

    if(input.dim <= 0) {
      cout << endl << "   Not allowed  matrix dimension - try once more !";
      continue;
    }
    break;       // all necessary data are receivec
  } // end while() loop

} // End: function inputData()

    /*
    ** The function                           
    **      TID time_step(..)                 
    ** calculates start and stop time and returns the difference.                
    ** Input data:                            
    **    int number  = 1 for start - zero return values     
    **                = 2 for stop  - return time from last start
    **    TID run     - returns the time difference. 
    */

TID time_step(int num)
{
  unsigned long long int
                           num_sec;

  static long
                           zsec = 0, zusec = 0;
  double
                           retval;
  TID
                           ex_time;
  struct timeval
                           tp;
  struct timezone
                           tpz;

  if(num == 1) {              // initialization of time
    gettimeofday(&tp, &tpz);

    zsec  = tp.tv_sec;
    zusec = tp.tv_usec;

    ex_time.sec  = 0;
    ex_time.min  = 0;
    ex_time.hour = 0;
  }
  else if(num == 2) {
    gettimeofday(&tp, &tpz);

    retval = (double)(tp.tv_sec - zsec) + (tp.tv_usec - zusec) * 0.000001;

    num_sec = (unsigned long long int)retval;
    ex_time.sec  = num_sec % 60;
    ex_time.min  = num_sec / 60;
    ex_time.hour = ex_time.min/ 60;
    ex_time.min  = ex_time.min % 60;
  }
  else {
    printf("\n\nError in function time_step(): ");
    printf("\nInput data num = %d is wrong !!\n\n", num);
    exit(1);
  }
  return ex_time;

} // End: function time_step()

