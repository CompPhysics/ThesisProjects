#include <math.h>
#include <iostream> 
#include <stdlib.h>
#include <blitz/array.h>

using namespace blitz;



int main( int argc, char* argv[] ) {
  
  Array<complex<double>,2> y(2,2);
  
  y = 0;

  y(0,0)=(1,2);
  //real(y)(0,0)=1;
  //imag(y)(0,0)=1;
  //imag(y)(0,0)=2;
  //imag(y)=1,1,1,1;
  
  //cout << "imag(y(0,0))=" << imag(y(0,0)) <<"\n";
  //cout << real(y) << endl;
  //cout << imag(y) << endl;
    
  cout << y(0,0) <<endl;
  
 

}//main
