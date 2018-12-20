#include <iostream>
#include <fstream>
#include "Complex.h"
#include "State_Vector.h"
#include "All_States.h"
#include "lib.h"

// file        : main.cpp
// description : Definition and solving energy of a quantum dot using Hartree-Fock method

using namespace std;


int main(int argc, char* argv[])
{

  double double_g;
  int int_n;
  while(true)
    {
      printf("test of the gamma function, give a real number: \n");
      cin >> double_g;
      printf("result gamma(%f) = %f \n", double_g,gamma(double_g));
      
      printf("test of the factorial function, give an integer: \n");
      cin >> int_n;
      printf("result (%d)! = %d \n", int_n,factorial(int_n));
    }
  // while(true)
//     {
//       int s_in;
//       printf("Enter an index for a state among all possible state: (so less than %d )\n",totalNbStates);
//       cin >> s_in;
//       printf("state %d: (n=%d,l=%d,ml=%d,ms=%d)\n",s_in, TableOfStates[s_in][0], TableOfStates[s_in][1], TableOfStates[s_in][2], TableOfStates[s_in][3]);

//     }

}
