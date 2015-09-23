#ifndef SolidHarmonicsFuncs_IS_INCLUDED
#define SolidHarmonicsFuncs_IS_INCLUDED

#include "../Domain/Domain.h"
#include <iostream>
#include <iomanip>

#ifndef sqr_IS_INCLUDED
#define sqr_IS_INCLUDED
inline double sqr(double x) {return x*x;};
#endif


// ****************************************************************
// *                     SOLIDHARMONICSFUNCS                      *
// ****************************************************************
template <class Param>
class SolidHarmonicsFuncs {
 protected:
  double constant;
 public:
  SolidHarmonicsFuncs() {}
  virtual double SH(Param& coordinate, int currentNumber) {return 0;}
};

// ***************************** SH00 ****************************
template <class Param>
class SH00 : public SolidHarmonicsFuncs<Param> {
 public:
  SH00() {}
  virtual double SH(Param& coordinate, int currentNumber) {return 1;}
};

// **************************** SH1m1 ****************************
template <class Param>
class SH1m1 : public SolidHarmonicsFuncs<Param> {
 public:
  SH1m1() {}
  virtual double SH(Param& coordinate, int currentNumber) 
    {return coordinate(0, currentNumber);}
};

// ***************************** SH10 ****************************
template <class Param>
class SH10 : public SolidHarmonicsFuncs<Param> {
 public:
  SH10() {}
  virtual double SH(Param& coordinate, int currentNumber) 
    {return coordinate(2, currentNumber);}
};

// **************************** SH1p1 ****************************
template <class Param>
class SH1p1 : public SolidHarmonicsFuncs<Param> {
 public:
  SH1p1() {}
  virtual double SH(Param& coordinate, int currentNumber) 
    {return coordinate(1, currentNumber);}
};

// **************************** SH2m2 ****************************
template <class Param>
class SH2m2 : public SolidHarmonicsFuncs<Param> {
 public:
  SH2m2() {constant = 1/2*sqrt(3.);}
  virtual double SH(Param& coordinate, int currentNumber) 
    {return constant*(sqr(coordinate(0, currentNumber)) - sqr(coordinate(1, currentNumber)));}
};

// **************************** SH2m1 ****************************
template <class Param>
class SH2m1 : public SolidHarmonicsFuncs<Param> {
 public:
  SH2m1() {constant = sqrt(3.);}
  virtual double SH(Param& coordinate, int currentNumber) 
    {return constant*coordinate(0, currentNumber)*coordinate(2, currentNumber);}
};

// ***************************** SH20 ****************************
template <class Param>
class SH20 : public SolidHarmonicsFuncs<Param> {
 public:
  SH20() {}
  virtual double SH(Param& coordinate, int currentNumber) 
    {return (3*sqr(coordinate(2, currentNumber)) - sqr(coordinate.r(currentNumber)))/2;}
};

// **************************** SH2p1 ****************************
template <class Param>
class SH2p1 : public SolidHarmonicsFuncs<Param> {
 public:
  SH2p1() {constant = sqrt(3.);}
  virtual double SH(Param& coordinate, int currentNumber) 
    {return constant*coordinate(1, currentNumber)*coordinate(2, currentNumber);}
};

// **************************** SH2p2 ****************************
template <class Param>
class SH2p2 : public SolidHarmonicsFuncs<Param> {
 public:
  SH2p2() {constant = sqrt(3.);}
  virtual double SH(Param& coordinate, int currentNumber) 
    {return  constant*coordinate(0, currentNumber)*coordinate(1, currentNumber);}
};

#endif
