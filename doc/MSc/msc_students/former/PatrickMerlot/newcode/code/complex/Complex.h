#ifndef Complex_h_IS_INCLUDED
#define Complex_h_IS_INCLUDED
#include <iostream>

// file        : Complex.h
// description : definition of class Complex

class Complex;    // tell C++ that this class name exists

class Complex
{
private:
  double re, im; // real and imaginary part
  
public:                // members visible also outside the class
  Complex (double _re = 0.0, double _im = 0.0);// Complex a(4,3);
  Complex (const Complex& c);                  // Complex q(a);
  ~Complex () {}                               // Destructor
  Complex& operator= (const Complex& c);       // a = b;
  double Re () const;                          // double real_part = a.Re();
  double Im () const;                          // double imag_part = a.Im();
  double abs () const;                         // double m = a.abs(); //modulus
  Complex conj () const;                       // Complex m = a.conj(); //conjugate

 
  Complex operator+ (const Complex& a) const;  // Sum of 2 complex
  Complex& operator+= (const Complex& a);      // iterative sum of complex
  Complex& operator-= (const Complex& a);      // iterative substraction of complex
  Complex operator- (const Complex& a) const;  // Substraction of 2 complex
  Complex operator* (const Complex& a)const;   // Multiplication of 2 complex
  Complex operator* (const double& a) const;   // Multiplication by a real
  Complex operator/ (const double& a) const;   // Division by a real
  Complex operator/ (const Complex& a) const;  // Division by a complex
  friend std::ostream& operator<< (std::ostream& os, const Complex& c); // write a complex to screen
  friend std::istream& operator>> (std::istream& is, Complex& c); // read in a complex 
};

#endif
