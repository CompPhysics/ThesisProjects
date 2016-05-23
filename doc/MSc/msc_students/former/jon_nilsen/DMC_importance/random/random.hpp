#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <cmath>
#include <iostream>

class Random{
private:
  long idum;

public:
  
  Random(){idum=-1;}
  Random(long idum_){idum=-(long)(fabs(idum_));}
  
  // normal random variate generator
  double gran(double s=.1, double m=.5);

  // uniform dist random generator
  double ran1();
  
  // slower and better uniform dist random generator
  double ran2();
  
};

#endif
