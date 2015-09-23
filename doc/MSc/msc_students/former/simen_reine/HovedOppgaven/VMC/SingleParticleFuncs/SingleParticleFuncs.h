#ifndef SingleParticleFuncs_IS_INCLUDED
#define SingleParticleFuncs_IS_INCLUDED

#include <cmath>
#include <iostream>
#include <fstream>
#include "../Ref/Ref.h"
#include "../Coor/Coor.h"
#include "../Domain/Domain.h"
#include "../SolidHarmonics/SolidHarmonics.h"
#include "../STOBasis/STOBasis.h"

// ****************************************************************
// *                      SINGELPARTICLEFUNC                      *
// ****************************************************************
template <class Param>
class SingleParticleFunc {

 protected:
  double*  param;
  virtual  inline double phi(Param& coordinate) {return 0;}
  int      centerNr;
  int      currentNr;

 public:
  SingleParticleFunc() {}
  virtual ~SingleParticleFunc() {}
  virtual void   attach(int _numParams, double* _param) 
    { param=_param;}
  virtual void   setNumberDimensions(int numDimensions)
    { centerNr = 2*numDimensions;}
  virtual double operator()(Param& coordinate) 
    { currentNr = centerNr;
    return phi(coordinate);}
  virtual double operator()(Param& coordinate, int nr) 
    { currentNr = nr;
    return phi(coordinate);}
  virtual void readHFparams(ifstream& ifile) {;}
};


//************************************************************************
//*                  3-dimensional hydrogen orbitals                     *
//************************************************************************
//*********** Hydr1s *********
template <class Param>
class Hydr1s : public SingleParticleFunc<Param> {
 protected:
  virtual inline double phi(Param& coordinate) {
    return exp(- param[0] * coordinate.r(currentNr));
  }
 public:
  Hydr1s() {}
  virtual ~Hydr1s() {}
};
//*********** Hydr2s *********
template <class Param>
class Hydr2s : public SingleParticleFunc<Param> {
 protected:
  virtual inline double phi(Param& coordinate) {
    return ((2 - param[0]*coordinate.r(currentNr)) * exp(-0.5*coordinate.r(currentNr)));
  }
 public:
  Hydr2s() {}
  virtual ~Hydr2s() {}
};
//*********** Hydr2px *********
template <class Param>
class Hydr2px : public SingleParticleFunc<Param> {
 protected:
  virtual inline double phi(Param& coordinate) {
    return param[0]*coordinate(0, currentNr) * exp(-0.5*coordinate.r(currentNr)*param[0]);
  }
 public:
  Hydr2px() {}
  virtual ~Hydr2px() {}
};
//*********** Hydr2py *********
template <class Param>
class Hydr2py : public SingleParticleFunc<Param> {
 protected:
  virtual inline double phi(Param& coordinate) {
    return param[0]*coordinate(1, currentNr) * exp(-0.5*coordinate.r(currentNr)*param[0]);
  }
 public:
  Hydr2py() {}
  virtual ~Hydr2py() {}
};
//*********** Hydr2pz *********
template <class Param>
class Hydr2pz : public SingleParticleFunc<Param> {
 protected:
  virtual inline double phi(Param& coordinate) {
    return param[0]*coordinate(2, currentNr) * exp(-0.5*coordinate.r(currentNr)*param[0]);
  }
 public:
  Hydr2pz() {}
  virtual ~Hydr2pz() {}
};


//************************************************************************
//*                3-dimensional Hartree-Fock orbitals                   *
//************************************************************************
//
template <class Param>
class HF : public SingleParticleFunc<Param> {
 protected:
  double*                     coefficients;
  int                         numCoeff;
  int                         numParams;
  Ref<SolidHarmonics<Param> > solidHarmonics;
  Ref<STOBasis<Param> >       STOfuncs;

  virtual inline double z1(Param& coordinate) {return 0;}
  virtual inline double z2(Param& coordinate) {return 0;}
  virtual inline double z3(Param& coordinate) {return 0;}
  virtual inline double z4(Param& coordinate) {return 0;}
  virtual inline double z5(Param& coordinate) {return 0;}
  virtual inline double z6(Param& coordinate) {return 0;}
  virtual inline double z7(Param& coordinate) {return 0;}
  virtual inline double z8(Param& coordinate) {return 0;}
  virtual inline double z9(Param& coordinate) {return 0;}
  virtual inline double z10(Param& coordinate) {return 0;}

  virtual inline double phi2(Param& coordinate) {
    return param[0]*z1(coordinate) + param[1]*z2(coordinate) +
      param[2]*z3(coordinate) + param[3]*z4(coordinate) +
      param[4]*z5(coordinate) + param[5]*z6(coordinate);
  }
  virtual inline double phi(Param& coordinate) {
    double result = 0;
    for (int i=0; i<numCoeff; i++) 
      result += coefficients[i] * solidHarmonics()(coordinate,i, currentNr) 
	* STOfuncs()(coordinate, i, currentNr);
    return result;
  }
 public:
  HF() {}

  virtual void readHFparams(ifstream& ifile) 
    {
      ifile >> numCoeff;

      int*    OrbN = new int[numCoeff];
      int     OrbL;
      int     OrbM;
      double* Xsi  = new double[numCoeff];
      coefficients = new double[numCoeff];

      for (int i=0; i<numCoeff; i++)
	ifile >> OrbN[i];
      ifile >> OrbL;
      ifile >> OrbM;
      for (int i=0; i<numCoeff; i++)
	ifile >> Xsi[i];
      for (int i=0; i<numCoeff; i++)
	ifile >> coefficients[i];

      solidHarmonics = new SolidHarmonics<Param>;
      STOfuncs       = new STOBasis<Param>;
      solidHarmonics().init(numCoeff, OrbL, OrbM);
      STOfuncs().init(numCoeff, OrbN, OrbL, Xsi);
      delete OrbN;
      delete Xsi;
    }
};

#endif
