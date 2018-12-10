#include <iostream>
#include "Vector_Complex.h"
#include "Complex.h"
#include <fstream>

using namespace std;

int main()
{
  // Declaration of a vector
  cout << "Vector_Complex a(5); a= ";
  Vector_Complex a(5);
  a.print(std::cout);
  std::cout << '\n';

  // definition of a
  for (int i=1; i<= a.size(); i++) {
  Complex b(10-i, i*i);
    a(i) = b;
  }
  cout << "fill in vector a with new complex a(i)=b; a= ";
  a.print(std::cout);
  std::cout << '\n';

  // b = a
  Vector_Complex b;
  b = a;
  cout << "b = a; b= ";
  b.print(std::cout);
  std::cout << '\n';

  // e = b(4); read entry from the vector
  Complex e = b(4);
  cout << "e = b(4) = " << e;
  std::cout << '\n'; 

  // Reading the real part of an entry of a complex vector
  double f = b(3).Re();
  cout << "f = b(3).Re() = " << f << '\n'; 

  // Division by a double
  Vector_Complex g=a/7;
  cout << "Vector_Complex g = a/7 = ";
  g.print(std::cout);
  std::cout << '\n';

  // Multiplication by a double and re-assignement
  g=4*g;
  cout << "Vector_Complex g = 4*g = ";
  g.print(std::cout);
  std::cout << '\n';

  // Inner product of 2 Vector_Complex
  Complex h = a.inner(g);
  cout << "Complex h=inner(a,g) = " << h << '\n';

  // Read the size of a Vector_Complex
  int Size = g.size();
  cout << "int Size = g.size() = " << Size << '\n';

  // Change dimension of a vector
  cout << " Enter the new length of vector g (integer please) \n";
  int newLength;
  cin >> newLength;
  g.redim(newLength);
  cout << "g.redim(newLength); g=  ";
  g.print(std::cout);
  std::cout << '\n';

  // read from file a Vector_Complex
  cout << "Test for reading a vector, please enter a vector of complex numbers VEC= ";
  Vector_Complex VEC;
  VEC.scan(cin);
  cout << "\n twice your vector is 2*VEC = ";
  (2*VEC).print(std::cout);
  std::cout << '\n';
  
  // Read from file a Vector_Complex
  cout << "Read what's in the EXAMPLE file \n";
  const char* filename = "example_vector.txt";
  std::ifstream ifile(filename);
  Vector_Complex fileVEC;
  fileVEC.scan(ifile);
  cout << "\n the vector inside the file is:  ";
  fileVEC.print(std::cout);
  std::cout << '\n';

  return 0;
}
