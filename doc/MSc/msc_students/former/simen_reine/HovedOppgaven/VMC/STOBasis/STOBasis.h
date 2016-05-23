#ifndef STOBasis_IS_INCLUDED
#define STOBasis_IS_INCLUDED

#include "../STOBasisFuncs/STOBasisFuncs.h"

// ****************************************************************
// *                          STOBASIS                            *
// ****************************************************************
template <class Param>
class STOBasis {

 protected:
  STOBasisFuncs<Param>* STOfuncs;
  double constant, expConst;

 public:
  STOBasis() {}
  void init(int num, int* orbN, int orbL, double* Xsi) 
    {
      STOfuncs    = new STOBasisFuncs<Param>[num];
      for (int i=0; i<num; i++) 
	STOfuncs[i].initConsts(orbN[i], orbL, Xsi[i]);
    }

  virtual double operator ()(Param& coordinate, int num, int currentNumber)
    { return STOfuncs[num](coordinate, currentNumber); }

};

#endif
