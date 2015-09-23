#include <iostream>
#include <string>
#include <blitz/array.h>   

using namespace std;
using namespace blitz;
 
int main(int argn, char* args[]) {
 
  int num;
  cin >> num;
  double sum = 0;
  double doubleSum = 0;

  Array<double,1> A;

  for (int i=0; i<num; i++) {
    double input;
    cin >> input;
    sum       += input;
    doubleSum += input*input;
  }
  cerr << "Average " << sum/num << " ± " 
       << sqrt((abs(doubleSum/num- sum*sum/num/num))/num) << endl;

  return 0;
}
