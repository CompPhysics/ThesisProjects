#ifndef SolidHarmonics_IS_INCLUDED
#define SolidHarmonics_IS_INCLUDED

#include "../Domain/Domain.h"
#include <iostream>
#include <iomanip>
#include "../SolidHarmonicsFuncs/SolidHarmonicsFuncs.h"

#ifndef sqr_IS_INCLUDED
#define sqr_IS_INCLUDED
inline double sqr(double x) {return x*x;};
#endif

// ****************************************************************
// *                       SOLIDHARMONICS                         *
// ****************************************************************
template <class Param>
class SolidHarmonics {

 protected:
  SolidHarmonicsFuncs<Param>** SHfuncs;

 public:
  SolidHarmonics() {}
  void init(int num, int orbL, int orbM) 
    {
      SHfuncs    = new SolidHarmonicsFuncs<Param>*[num];
      for (int i=0; i<num; i++) {
	if ( (orbL == 0)&(orbM == 0) )
	  SHfuncs[i] = new SH00<Param>;
	else if ( (orbL == 1)&(orbM == -1) )
	  SHfuncs[i] = new SH1m1<Param>;
	else if ( (orbL == 1)&(orbM == 0) )
	  SHfuncs[i] = new SH10<Param>;
	else if ( (orbL == 1)&(orbM == 1) )
	  SHfuncs[i] = new SH1p1<Param>;
	else if ( (orbL == 2)&(orbM == -2) )
	  SHfuncs[i] = new SH2m2<Param>;
	else if ( (orbL == 2)&(orbM == -1) )
	  SHfuncs[i] = new SH2m1<Param>;
	else if ( (orbL == 2)&(orbM == 0) )
	  SHfuncs[i] = new SH20<Param>;
	else if ( (orbL == 2)&(orbM == 1) )
	  SHfuncs[i] = new SH2p1<Param>;
	else if ( (orbL == 2)&(orbM == 2) )
	  SHfuncs[i] = new SH2p2<Param>;
      }
    }
  virtual double operator ()(Param& coordinate, int num, int currentNumber)
    { return SHfuncs[num]->SH(coordinate, currentNumber); }

};

#endif
