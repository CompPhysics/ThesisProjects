#ifndef Vector_Complex_CPP
#define Vector_Complex_CPP

// file        : Vector_Complex.cpp
// description : implentation of class Vector_Complex

#include <iostream>
#include <cmath>
#include "Vector_Complex.h"
#include "Complex.h"

using namespace std;

// Allocate memory for n Complex in memory
void Vector_Complex:: allocate (int n)
{
  length = n;
  A = new Complex [n];
}

// Free dynamic memory
void Vector_Complex:: deallocate ()
{
  delete [] A;
}

// Declare a vector of length 0
Vector_Complex:: Vector_Complex ()
{
  A = NULL; length  = 0;
}

// Destructor to specify how to destroy the object
Vector_Complex:: ~Vector_Complex ()
{
  deallocate();
}

// Declare a vector of length n
Vector_Complex:: Vector_Complex (int n)
{
  allocate(n);
}

// v.redim(n); redimensioning the length
bool Vector_Complex:: redim (int n)
{
  if (length == n)
    return false; // no need to allocate anything
  else {
    if (A != NULL) {
      // this object has already allcoated memory
      deallocate();
    }
    allocate(n);
    return true; // the length was changed
  }
}

// v = w; (assigment operator, v and w are Vector_Complex objects)
Vector_Complex& Vector_Complex:: operator= (const Vector_Complex& w)
{
  redim (w.size()); // make v as long as w
  int i;
  for(i = 0; i < length; i++)
    A[i] = w.A[i];
  return *this;
  // return of *this allows nested assignement: u = v = u_vec = v_vec ...
}

// Create a new vector as a copy of an existing one
Vector_Complex:: Vector_Complex (const Vector_Complex& w)
{
  allocate(w.size());
  *this = w;
}

// Complex e = v(i);(extract an entry from a vector)
Complex Vector_Complex:: operator() (int i) const
{
  if (i < 1 || i > length)
    std::cerr << "Complex Vector_Complex::operator(), illegal index, i=" << i << '\n';
  return A[i-1]; // base index is 1 (not 0 as in C/C++)
}

// v(i)=a; (assign a number to an entry indexing from 1 instead of 0)
Complex& Vector_Complex:: operator() (int i)
{
 
  if (i < 1 || i > length)
    std::cerr << "Complex& Vector_Complex::operator(), illegal index, i=" << i << '\n';
  return A[i-1]; // base index is 1 (not 0 as in C/C++)
}

// v.print(cout);(write a vector to the screen/file)
void Vector_Complex:: print(ostream& o) const
{
  o << "(";
  for(int i = 1; i < length; i++)
    o << (*this)(i) << ",";
  o << (*this)(length) << ")";
}

ostream& operator<< (ostream & o, const Vector_Complex & v)
{
  v.print(o);
  return o;
}

// a = inner(v,w); (inner product)
Complex Vector_Complex:: inner (const Vector_Complex& w) const
{
  Complex sum(0,0);
  for (int i=0; i < length; i++)
    sum += A[i]*w.A[i];
  return sum;
}


// for reading a vector from a file/screen
void Vector_Complex:: scan (std::istream& is)  
{
  // input format for a Complex: f or (f) or (f,f)
  // input format for a Vector_Complex: (f1, f2, ...)
  //              or ((f1),(f2),...) or ((re1,im1),(re2,im2),...)
  Vector_Complex vec;
  Complex * buf, * buf2;
  int blen = 0x10;// default buffer length (then multiple of it if bigger)
  int sizeVec = 0;// size of the vector in file 
  buf = new Complex[blen];
  Complex comp; // courant complex
  char str = 0; // courant character
  bool OK;// check presence or not of a complex

  is >> str;
  if(str=='(') { //start with a parenthesis, can be an empty vector or not
    OK = (is >> comp);
    if (OK==false) { // is it an empty vector?
      cerr << "The vector is empty" << '\n';
    }
    else {      // there is at leat one complex number?
      while(OK==true) // while there is a complex
	{
	  sizeVec++;
	  if(sizeVec >= blen)
	    {
	      blen *= 2;
	      buf2 = new Complex[blen];
	      for(int j = 0; j < sizeVec; j++) buf2[j]=buf[j];
	      Complex * tmp = buf2;
	      buf2 = buf;
	      buf = tmp;
	      delete [] buf2;
	    }
	  buf[sizeVec-1] = comp;
	  is >> str;
	  if(str == ',') {
	    OK = (is >> comp);
	  }
	  else if(str != ')') {
	    is.clear(ios_base::failbit); // in case of format error: set state
	    OK = false;
	  }
	  else
	    OK = false;
	}
      // no more complex to read
    }
  // now copy buf into a Vector_Complex
  vec.redim(sizeVec);
  for(int j = 0; j < sizeVec; j++) vec(j+1)=buf[j];
  delete [] buf;
  }
  else { // no parenthesis, if the vector contains only one real number or nothing
    if (is) {
      is.putback(str);
      vec.redim(1);
      is >> vec(1);
    }
    else {
      cout << "The file is empty" << '\n'; // empty file
    }
  }


  (*this) = vec;
}


// Sum of 2 vectors
Vector_Complex Vector_Complex:: operator+ (const Vector_Complex& w) const
{
  Vector_Complex tmp(w.size());
  if ( w.size() != (*this).size() )
    std::cerr << "Vector_Complex::operator+(), illegal size, size=" << w.size();
  for (int i=0; i< w.size(); i++)
    tmp.A[i] = A[i] + w.A[i];
  return tmp;
} 

// Substraction of 2 vectors
Vector_Complex Vector_Complex:: operator- (const Vector_Complex& w) const
{
  Vector_Complex tmp(w.size());
  if ( w.size() != (*this).size() )
    std::cerr << "Vector_Complex::operator-(), illegal size, size=" << w.size();
  for (int i=0; i< w.size(); i++)
    tmp.A[i] = A[i] - w.A[i];
  return tmp;
} 

// u = v*a; Multiplication with a double
Vector_Complex Vector_Complex:: operator* (const double& a) const
{
  Vector_Complex tmp( length );
  for (int i=0; i< length; i++)
    tmp.A[i] = A[i]*a;
  return tmp;
} 

// u = v/a; Division with a double
Vector_Complex Vector_Complex:: operator/ (const double& a) const
{
  Vector_Complex tmp( length );
  for (int i=0; i< length; i++)
    tmp.A[i] = A[i]/a;
  return tmp;
} 

// u = v*a; Multiplication with a Complex
Vector_Complex Vector_Complex:: operator* (const Complex& a) const
{
  Vector_Complex tmp( length );
  for (int i=0; i< length; i++)
    tmp.A[i] = A[i]*a;
  return tmp;
} 

// u = v/a; Division with a Complex
Vector_Complex Vector_Complex:: operator/ (const Complex& a) const
{
  Vector_Complex tmp( length );
  for (int i=0; i< length; i++)
    tmp.A[i] = A[i]/a;
  return tmp;
} 


//  u = a*v; Multiplication by a real scalar on the left
Vector_Complex operator* (const double& a, const Vector_Complex& v)
{
  return v*a;
}

// u = a*v; Multiplication by a Complex scalar on the left
Vector_Complex operator* (const Complex& a, const Vector_Complex& v)
{
  Vector_Complex tmp( v.size() );
  for (int i=1; i<= v.size(); i++)
    tmp(i) = v(i)*a;
  return tmp;
}


#endif
