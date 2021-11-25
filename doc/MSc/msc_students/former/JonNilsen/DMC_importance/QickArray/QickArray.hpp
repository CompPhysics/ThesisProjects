#ifndef __QICKARRAY__
#define __QICKARRAY__

#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>

#define _TRUE_ 1
#define _FALSE_ 0
//#define _DEBUG_
//#define _BOUNDS_CHECK_

class QickArray {
private:
  double *data_;              //contains actual data
  int dimensions_[6];        //6 dimensions ought to be enough...
  int n_dims_;                  //contains no. of used dimensions
  int is_initialized_;          //_always_ indicates whether mem is alloced.
  int no_of_elements_;       //stores length of data-array.

  int step1_, step2_, step3_, step4_, step5_;   //stores indexing increments.

  int allocate_memory(void);
  int delete_memory(void);

  std::string sub_python_string( int dim_no, int indices[], int compact ) const;

public:
  QickArray(void);              //default constructor.
  QickArray(int i1, int i2=0, int i3=0, int i4=0, int i5=0, int i6=0); //a more elaborate constructor.
  int redim(int i1, int i2=0, int i3=0, int i4=0, int i5=0, int i6=0); //brand new array.
  ~QickArray(void);
  QickArray(const QickArray& clone_me); //a copy constructor.

  int no_of_elements(void) const;

  int *get_dimensions_ptr(void);      //this one is dangerous!
  int get_dimension_info(int *i1=NULL, int *i2=NULL, int *i3=NULL, int *i4=NULL, int *i5=NULL, int *i6=NULL);
  int is_initialized(void);
  double *get_data_ptr(void);         //also dangerous!

  int fill_data(double value);   //fill data...


  int index(int i1, int i2=0, int i3=0, int i4=0, int i5=0, int i6=0) const;
  int logical_index(int idx, int *i1=NULL, int *i2=NULL, int *i3=NULL, int *i4=NULL, int *i5=NULL, int *i6=NULL);

  const double& get(int i1, int i2=0, int i3=0, int i4=0, int i5=0, int i6=0) const;
  int put(double value, int i1, int i2=0, int i3=0, int i4=0, int i5=0, int i6=0);

  QickArray& operator=(const QickArray &clone_me);  //assignment
  QickArray& operator+=(const QickArray &A);  //memberwise arithmetics.
  QickArray& operator-=(const QickArray &A);  // ...
  QickArray& operator*=(const QickArray &A);  // ...
  QickArray& operator/=(const QickArray &A);  // ...

  QickArray& operator+=(const double x);  //memberwise arithmetics with scalar.
  QickArray& operator-=(const double x);  // ...
  QickArray& operator*=(const double x);  // ...
  QickArray& operator/=(const double x);  // ...


  //  double& operator()(int i1, int i2=0, int i3=0, int i4=0, int i5=0, int i6=0); //indexing
  //  double operator()(int i1, int i2=0, int i3=0, int i4=0, int i5=0, int i6=0) const; //indexing

  inline double operator()(int i1) const; 
  inline double operator()(int i1, int i2) const;
  inline double operator()(int i1, int i2, int i3) const;
  inline double operator()(int i1, int i2, int i3, int i4) const;
  inline double operator()(int i1, int i2, int i3, int i4, int i5) const;
  inline double operator()(int i1, int i2, int i3, int i4, int i5, int i6) const;
  
  inline double& operator()(int i1); 
  inline double& operator()(int i1, int i2);
  inline double& operator()(int i1, int i2, int i3);
  inline double& operator()(int i1, int i2, int i3, int i4);
  inline double& operator()(int i1, int i2, int i3, int i4, int i5);
  inline double& operator()(int i1, int i2, int i3, int i4, int i5, int i6);
  
  //these were dangerous! commented out!
  //QickArray& operator*(const QickArray &B); //matrix multiplication
  //QickArray& operator+(const QickArray &B); //matrix addition
  //QickArray& operator-(const QickArray &B); //matrix subtraction

  //QickArray& operator*(const double x); //scalar multiplication
  //QickArray& operator/(const double x); //scalar division
  //QickArray& operator+(const double x); //scalar addition
  //QickArray& operator-(const double x); //scalar subtraction


  std::string python_string(int compact=_TRUE_) const;

  std::string summary(void);

  //some statistical functions.
  double min(int *i1=NULL, int *i2=NULL, int *i3=NULL, int *i4=NULL, int *i5=NULL, int *i6=NULL);
  double max(int *i1=NULL, int *i2=NULL, int *i3=NULL, int *i4=NULL, int *i5=NULL, int *i6=NULL);
  double average(void);
  double variance(void);
  double stddev(void);

};



inline double QickArray::operator()(int i1) const {
  return data_[i1];
}
inline double QickArray::operator()(int i1, int i2) const {
  return data_[i1 + i2*step1_];
}
inline double QickArray::operator()(int i1, int i2, int i3) const {
  return data_[i1 + i2*step1_ + i3*step2_];
}
inline double QickArray::operator()(int i1, int i2, int i3, int i4) const {
  return data_[i1 + i2*step1_ + i3*step2_ + i4*step3_];
}
inline double QickArray::operator()(int i1, int i2, int i3, int i4, int i5) const {
  return data_[i1 + i2*step1_ + i3*step2_ + i4*step3_ + i5*step4_];
}
inline double QickArray::operator()(int i1, int i2, int i3, int i4, int i5, int i6) const {
  return data_[i1 + i2*step1_ + i3*step2_ + i4*step3_ + i5*step4_ + i6*step5_];
}


inline double& QickArray::operator()(int i1) {
  return data_[i1];
}
inline double& QickArray::operator()(int i1, int i2) {
  return data_[i1 + i2*step1_];
}
inline double& QickArray::operator()(int i1, int i2, int i3) {
  return data_[i1 + i2*step1_ + i3*step2_];
}
inline double& QickArray::operator()(int i1, int i2, int i3, int i4) {
  return data_[i1 + i2*step1_ + i3*step2_ + i4*step3_];
}
inline double& QickArray::operator()(int i1, int i2, int i3, int i4, int i5) {
  return data_[i1 + i2*step1_ + i3*step2_ + i4*step3_ + i5*step4_];
}
inline double& QickArray::operator()(int i1, int i2, int i3, int i4, int i5, int i6) {
  return data_[i1 + i2*step1_ + i3*step2_ + i4*step3_ + i5*step4_ + i6*step5_];
}


//
// this function convert the set of logical indices
// to absolute index to data_[]
//
inline int QickArray::index(int i1, int i2, int i3, int i4, int i5, int i6) const {
#ifdef _DEBUG_
  cout << "entering index(int,...) ";
#endif
  
#ifdef _BOUNDS_CHECK_

  //Make sure user doesn't out-index the array...
  if (dimensions_[0]) i1 %= dimensions_[0];
  if (dimensions_[1]) i2 %= dimensions_[1];
  if (dimensions_[2]) i3 %= dimensions_[2];
  if (dimensions_[3]) i4 %= dimensions_[3];
  if (dimensions_[4]) i5 %= dimensions_[4];
  if (dimensions_[5]) i6 %= dimensions_[5];

#endif


  /*

  int temp; 


  temp = 0;

  switch (n_dims_) {
  case 1:
    temp = i1;
    break;
  case 2:
    temp = (i2*dimensions_[0] + i1);
    break;
  case 3:
    temp = i3*dimensions_[1] + i2;
    temp = temp * dimensions_[0] + i1;
    break;
  case 4:
    temp = (i4*dimensions_[2] + i3);
    temp = temp * dimensions_[1] + i2;
    temp = temp * dimensions_[0] + i1;
    break;
  case 5:
    temp = (i5*dimensions_[3] + i4);
    temp = temp * dimensions_[2] + i3;
    temp = temp * dimensions_[1] + i2;
    temp = temp * dimensions_[0] + i1;
    break;
  case 6:
    temp = (i6*dimensions_[4] + i5);
    temp = temp * dimensions_[3] + i4;
    temp = temp * dimensions_[2] + i3;
    temp = temp * dimensions_[1] + i2;
    temp = temp * dimensions_[0] + i1;
    break;
  }
  
  */

  int temp = 0;

  //if (step1_) temp += step1_ * i2;
  //if (step2_) temp += step2_ * i3;
  //if (step3_) temp += step3_ * i4;
  //if (step4_) temp += step4_ * i5;
  //  if (step5_) temp += step5_ * i6;

  temp = i1 + i2*step1_ + i3*step2_ + i4*step3_ + i5*step4_ + i6*step5_;

#ifdef _DEBUG_
  cout << "exiting index(int,...). (" << temp << ")" << endl;
#endif

  return temp;

}


#endif
