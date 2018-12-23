#ifndef LIN_EQ_IS_INCLUDED
#define LIN_EQ_IS_INCLUDED

    // ******   data declaration  ******* 

typedef   struct  {   // structure definition for execution time   
  unsigned long long int
                          tick,
                           sec,
                           min,
                          hour;
} TID; 

typedef   struct {
  int
           matrType,
           dataType,
                dim;
} INPUTDATA;

         // *** function declaration ***

template <typename T> void solveEq(INPUTDATA data);
template <typename T> void matrixVal(INPUTDATA data, Array<T,2>& A, Array<T,1>& B);
template <typename T> void checkEq(Array<T,2>& A, Array<T,1>& X, Array<T,1>& B);


TID time_step(int num);



      // *** end function declation ***

     /*
     ** The template function 
     **      solveEq()
     ** solves the set of linear equations,
     ** check the results and print it to 
     ** standard output
     */

template <typename T>
void solveEq(INPUTDATA data)
{
  int            row,col;

  TID            ex_time;

  T              val;
  Array<T,2>     A(data.dim, data.dim),  A1(data.dim, data.dim);
  Array<T,1>     B(data.dim), X(data.dim);
  Array<int,1>   Index(data.dim);

      // fill the matrix A( , ) and vector B() with data

  matrixVal(data, A, B);

       /*
       ** In order to check the solution we save the original
       ** matrix A and vector B and use A1 and X
       ** in the calculation since ludcmp() and lubksb() destroy
       ** the original elements
       */

  A1 = A;               // Blitz asignements
  X  = B;

        // necessary parameter if ludcmp() is used to invert a matrix

  T permutation = static_cast<T>(1);

  time_step(1);                          // start clock

  ludcmp(A1, data.dim, Index, permutation);     // LU decomposition of A1()

  ex_time = time_step(2);                // read clock

  cout << endl << "Time used:  " << ex_time.hour << " hour    "
                                 << ex_time.min << " min      "
                                 << ex_time.sec << " sec      "
       << endl;                    

  cout << endl << "The coefficient matrix A1 after  ludcmp(A1 = " << A1  << endl;

  lubksb(A1, data.dim, Index, X);              // sovle the equations

  cout << endl << endl
       << "Solution to " << data.dim <<" linear equations"
       << endl <<endl;

  cout << endl << "The coefficient matrix A = " << A  << endl;  // Blitz object output
  cout << endl << "The right-hand matrix  B = " << B  << endl;
  cout << endl << "The solution           X = " << X  << endl;

  checkEq(A, X, B);

  cout << endl << endl
       << "Check the solution B - (A x X) = " << B
       << endl <<endl;

  A.free();    // release memory   -- Blitz methods     
  A1.free();
  B.free();
  X.free();

} // End: tenplate function solveEq()

  /*
  ** The template function 
  **      matrixVal()
  ** fills matrices A and B with random numbers
  ** and prints the result to standard output
  */

template <typename T>
void matrixVal(INPUTDATA data, Array<T,2>& A, Array<T,1>& B)
{
  int    idum = 10, ok = 1, row, col;

  switch(data.matrType) {
    case 1:  for(row = 0; row < A.rows(); row++) {
               for(col = 0; col < A.columns(); col++) {
		 A(row,col) = static_cast<T>(10.0 * (ran0<T>(idum) - 0.5));
	       }
             } 
             for(row = 0; row < B.rows(); row++) {
	       B(row) = static_cast<T>(10.0 * (ran0<T>(idum) - 0.5));
	     }
             break;
    case 2:  while(ok) {
	       cout << endl << endl
		    << "Type in data for the A matrix:";
	       for(row = 0; row < data.dim; row++) {
		 for(col = 0; col < data.dim; col++) {
                   cout << endl << "A(" << row << "," 
                        << col << ") = ";
		   cin >> A(row, col);
		 }
	       }
	       cout << endl << endl
	            << "Type in data for the B vector:";
	       for(row = 0; row < data.dim; row++) {
		 cout << endl <<"B(" << row <<") = ";
                 cin >> B(row);
	       }
	       cout << endl << "A = " << A  << endl;
	       cout << endl << "B = " << B  << endl;
	    
	       cout << endl << "Data ok (yes = 0, no = 1) ? : "; 
	       cin >> ok;
	     }
	    break;
  }  // end switch()
} // End: template function matrixVal()

   /*
   **  The template function 
   **         CheckEq()
   ** calculates 
   **         (A x X)  - B
   ** If X is the solution to the set of linear equations
   **         A x X  = B
   ** the return values in B are all ZERO 
   */  
  
template<typename T> 
void checkEq(Array<T,2>& A, Array<T,1>& X_vector, Array<T,1>& B)
{
  for(int row = 0; row < A.rows(); row++) {
    T temp = 0;
    for(int k = 0; k < A.columns(); k++) {
      temp += A(row,k) * X_vector(k);
    }
    B(row) -= temp;
  } // end row index

} // end: template function checkEq() 


#endif
