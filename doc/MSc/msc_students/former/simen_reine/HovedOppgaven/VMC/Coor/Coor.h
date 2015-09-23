#ifndef Coor_IS_INCLUDED
#define Coor_IS_INCLUDED

#include <iostream>
#include <cmath>
using namespace std;


// ****************************************************************
// *                            COOR                              *
// ****************************************************************
class Coor {

 protected:
  double* x;    // Cartesian coordinates
  double* _x;   // Pointer to current coordinate
  int     len;  // Number of cartesian dimensions


 public:
  Coor(double* x, int _len);
  Coor() {}
  virtual void           attach(double* __x, int _len);
  inline virtual void    resetPtr()          {_x = x;}
  inline virtual double& operator()()        {return *_x;}
  inline virtual double& operator()(int num) {return x[num];}
  inline virtual void    operator++(int)     {_x++;}
  virtual int            getLen()            {return len;}
};


// ****************************************************************
// *                           COORR                              *
// ****************************************************************
class CoorR : public Coor {

 protected:
  double _r;
  int    rIsCalculated;

 public:
  CoorR(double* x, int _len);
  CoorR();  
  virtual void    attach(double* x, int _len);
  inline virtual double& r()      {return _r;}
  inline virtual void    calculateR();
};


// ****************************************************************
// *                          COORSPIN                            *
// ****************************************************************
class CoorSpin : public CoorR {

 protected:
  int* _spin;

 public:
  CoorSpin(double* x, int _len, int* __spin);
  CoorSpin() {}
  virtual void attach(double* x, int _len, int* __spin);
  inline virtual int& spin()     {return *_spin;}
  inline virtual int  flipSpin() {return (*_spin *= -1);}

};


// ****************************************************************
// *                        COORSPINDIFF                          *
// ****************************************************************
class CoorSpinDiff : public CoorSpin {

 protected:
  // rDiff[0]             = sqrt( (x+h)^2 + y^2     + ... )
  // rDiff[1]             = sqrt( x^2     + (y+h)^2 + ... )
  // ...
  // rDiff[numDimensions] = sqrt( (x-h)^2 + y^2     + ... )
  // ...
  double* rDiff;
  double h, twoh;
  int twoLen;

 public:
  CoorSpinDiff(double* x, int _len, int* __spin, double _h);
  CoorSpinDiff() {}
  virtual void attach(double* x, int _len, int* __spin, double _h);
  virtual void calculateDiffs();
  // rdiff(0)             = sqrt( (x+h)^2 + y^2     + ... )
  // rdiff(1)             = sqrt( x^2     + (y+h)^2 + ... )
  // ...
  // rdiff(numDimensions) = sqrt( (x-h)^2 + y^2     + ... )
  // ...
  inline virtual double& r()           {return _r;}
  inline virtual double& r(int number) {return rDiff[number];}
  inline virtual double& operator()()  {return *_x;}
  inline virtual double& operator()(int num);
  inline virtual double  operator()(int num, int currentNr);
  inline virtual double* getX()        {return x;}
  inline virtual void    copy(CoorSpinDiff* copyCoor);
  inline virtual void    calculateR();
  virtual double rdiff(int currentNr)  {return rDiff[currentNr];}
};

#endif
