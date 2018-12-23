#ifndef Complex_CPP
#define Complex_CPP

// file        : Complex.cpp
// description : implentation of class Complex

#include <iostream>
#include <cmath>
#include "Complex.h"

using namespace std;

// Contructors
Complex:: Complex(double _re, double _im):re(_re), im(_im) {}

Complex:: Complex (const Complex& c )           // Complex q(a);
{ *this = c; }  

// a = b; Assignement
Complex& Complex:: operator= (const Complex& b) // a = b;
{
  re = b.re;
  im = b.im;
  return *this;
}

// Real and Imaginary part of a complex
double Complex:: Re () const { return re; }; // a copy of the real part is returned
double Complex:: Im () const { return im; }; // a copy of the imag part is returned

// absolute value of a complex
double Complex:: abs () const { return sqrt(re*re + im*im); } // modulus

// conjugate of a complex
Complex Complex:: conj () const { return Complex(re, -im); } // conjugate 

// addition of a complex
Complex Complex:: operator+ (const Complex& a) const
{
  return Complex(a.re+this->re,a.im+this->im);
}


// substraction of a complex
Complex Complex:: operator- (const Complex& a) const
{
  return Complex((*this).re-a.re,(*this).im-a.im);
}

// iterative sum of complex
Complex& Complex:: operator+= (const Complex& a)
{
  this->re += a.re;
  this->im += a.im;
  return *this;
}

// iterative substraction of complex
Complex& Complex:: operator-= (const Complex& a)
{
  this->re -= a.re;
  this->im -= a.im;
  return *this;
}

// Multiplication by a complex
Complex Complex:: operator* (const Complex& a) const
{
  return Complex(a.re*this->re-a.im*this->im, a.re*this->im+a.im*this->re);
}

// Multiplication by a double
Complex Complex:: operator* (const double& b) const
{
  return Complex(this->re*b, this->im*b);
}

// Division by a double
Complex Complex:: operator/ (const double& b) const
{
  return Complex(this->re/b, this->im/b);
}

// Division by a Complex
Complex Complex:: operator/ (const Complex& b) const
{
  return (*this) * b.conj() / (b.re*b.re + b.im*b.im);
}

// print a Complex to screen/file
ostream& operator<< (ostream& os, const Complex& c)
{
  return os << "(" << c.re << "," << c.im << ")";
}

// read in a Complex from screen/file
istream& operator>> (istream& is, Complex& a)
{
  // input format for a complex: f or (f) or (f,f)
  double re =0, im=0;
  char c=0;
    
  is>>c;
  if(c=='(') {
    is >> re >> c;
    if (c == ',') is >> im >> c;
    if (c != ')') is.clear(ios_base::failbit); // in case of format error: set state
  }
  else {
    is.putback(c);
    is >> re;
  }
  if (is) a= Complex(re, im);
  return is;
}

#endif
