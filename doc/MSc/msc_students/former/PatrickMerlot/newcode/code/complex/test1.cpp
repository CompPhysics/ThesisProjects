#include <iostream>
#include "Complex.h"

using namespace std;

int main()
{
  Complex a(0,1); // imaginary unit
  Complex b(2), c(3,-1);
  Complex q = b;

  std::cout << "(imaginary unit) a=" << a << ", b(2)=" << b << "q=b= " << q << "c(3,-1)= " << c << "\n";

  // DAXPY test
  q = a*c + b/a;
  std::cout << "q = a*c + b/a;" << q << "\n";
 
  // iterative sum
  Complex d(2,200);
  q += d;
  std::cout << "d(2,200); q += d" << q << "\n";
  
  // Real and Imaginary part
  std::cout << "Re(q)=" << q.Re() << ", Im(q)=" << q.Im() << "\n";

  // substraction of 2 complex
  cout << "b-a= " << b-a  << "\n";

  // division of 2 complex
  Complex e(-1), f(2.0, 0.01);
  cout << "e = " << e << " f= " << f << '\n';
  cout << "e/f" << e/f << '\n';

  // conjugate of a complex
  cout << "conj(f)= " << f.conj() << '\n';
  
  // print out a complex
  std::cout << "enter a Complex number: z=" << "\n";
  Complex z;
  std::cin >> z;
  std::cout << "f + z= " << f+z << "\n";


  return 0;
}
