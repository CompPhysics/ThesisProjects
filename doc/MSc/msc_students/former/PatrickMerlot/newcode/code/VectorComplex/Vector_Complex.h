#if !defined(Vector_Complex_h_IS_INCLUDED)
#define Vector_Complex_h_IS_INCLUDED

// file        : Vector_Complex.h
// description : Definition of the class Vector_Complex which contains a vector of complex values

#include <iostream>
#include "Complex.h"

class Vector_Complex;      // tell C++ that this class name exists

class Vector_Complex
{
private:
  Complex* A;              // vectors entries
  int length;              // length of the vector
  void allocate (int n);   // allocate memory, length=n
  void deallocate();       // free memory
public:                    // members visible also outside the class
  Vector_Complex ();                           // Vector_Complex c; (EMPTY CONSTRUCTOR)
  ~Vector_Complex ();                          // (Destructor, clean up dynamic memory)
  Vector_Complex (int n);                      // Vector_Complex v(n); (Vector of a specified length n)
  Vector_Complex (const Vector_Complex& w);    // Vector_Complex v(w); (copy of another vector w)
  int size() const { return length; };         // n = v.size(); (extract the lenght of a vector w)
  bool redim (int n);                          // v.redim(n); (Redimension a vector to length n)


  Vector_Complex& operator= (const Vector_Complex& w); // v = w; (assigment operator)
  Complex operator() (int i) const;                    // Complex e = v(i);(extract an entry)
  Complex& operator() (int i);// v(i)=a; (assign a number to an entry indexing from 1 instead of 0)
  
  Complex inner (const Vector_Complex& w) const;       // a = inner(v,w); (inner product)
  void print (std::ostream& o) const; // v.print(cout);(write a vector to the screen/file)
  void scan (std::istream& is);                        // for input from a file


  Vector_Complex operator+ (const Vector_Complex& w) const; // Sum of 2 vectors
  Vector_Complex operator- (const Vector_Complex& w) const; // Substraction of vector
 
  Vector_Complex operator* (const double& a) const;    // u = v*a; 
  Vector_Complex operator/ (const double& a) const;    // u = v/a;
  Vector_Complex operator* (const Complex& a) const;   // u = v*a;
  Vector_Complex operator/ (const Complex& a) const;   // u = v/a;

};

// Multiplication by a real scalar on the left
Vector_Complex operator* (const double& a, const Vector_Complex& v);  // u = a*v;

// Multiplication by a Complex scalar on the left
Vector_Complex operator* (const Complex& a, const Vector_Complex& v) ;  // u = a*v; 

// Print out a complex to file or to screen
std::ostream& operator<< (std::ostream & o, const Vector_Complex & v);

#endif
