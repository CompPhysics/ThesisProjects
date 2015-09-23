#ifndef Random2_IS_INCLUDED
#define Random2_IS_INCLUDED

#include <iostream>
#include <cstdio>

// ****************************************************************
// *                          RANDOM2                             *
// ****************************************************************
class Random2 {
 protected:
  long dummy;

 public:
  Random2(long seed) : dummy(seed) {}
  virtual ~Random2() {}
  virtual double getNum(void) {return 0;}
  virtual void   demo(int num, FILE* stream);
};

// ****************************************************************
// *                            RAN0                              *
// ****************************************************************
class Ran0 : public Random2 {
 public:
  Ran0(long seed) : Random2(seed) {}
  double getNum(void);
};

// ****************************************************************
// *                            RAN1                              *
// ****************************************************************
class Ran1 : public Random2 {
 public:
  Ran1(long seed) : Random2(seed) {}
  double getNum(void);
};

#endif
