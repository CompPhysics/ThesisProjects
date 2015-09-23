#ifndef Functor_IS_INCLUDED
#define Functor_IS_INCLUDED

#include "../Coor/Coor.h"
#include <string>
#include <iostream>
#include <iostream>
using namespace std;


template <class Param, class Return>
class Functor {
  typedef Return (Function)(Param&);
  Function* function;
  double* memory;

 public:
  Functor() {}
  Functor(Function& _function) : function(&_function) {}
  void attach(Function& _function) {function = &_function;}
  // memory=new double[100];}
  inline Return operator()(Param& coor) {return (*function)(coor);}
  void dumpFunctionPointer() {cerr << function;}
};

#endif
