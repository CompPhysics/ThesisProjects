#ifndef Func_IS_INCLUDED
#define Func_IS_INCLUDED

#include <string>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "../Coor/Coor.h"
#include "../Functor/Functor.h"
#include "../Ref/Ref.h"
#include "../Domain/Domain.h"
#include "../SingleParticleFuncs/SingleParticleFuncs.h"

using namespace std;


// ****************************************************************
// *                            FUNC                              *
// ****************************************************************
template <class Param>
class Func {

 protected:
  Ref<Param>     coordinate;
  SingleParticleFunc<Param>** function;
  double         h, h_half, h_double, h_inv, h_inv2, h_invdouble;
  double*        values;
  double*        result;
  double*        diffResult;
  double*        ddiffResult;
  int            numDimensions;

 public:
  Func(Param& _coordinate);
  Func();
  virtual ~Func();
  virtual void   init();
  virtual void   init(Param& _coordinate);
  virtual void   attachNumDimensions(int _numDimensions) 
    { numDimensions = _numDimensions; }
  virtual void   setCoordinate(Param& _coordinate);
  virtual void   calcValueCenter();
  virtual void   calcValueCenter(Param& _coordinate);
  virtual void   calcValueSides();
  virtual void   calcValueSides(Param& _coordinate);
  virtual double valuePt();
  virtual double diff();
  virtual double ddiff();
  virtual void   attachResult(double* __result);
  virtual void   attachDiffResult(double* __diffResult);
  virtual void   attachDdiffResult(double* __ddiffResult);
  virtual Param& getCoordinate() {return coordinate();}
  virtual void   summary();

 protected:
  virtual void   calcH();

};


// ****************************************************************
// *                           FUNCSET                            *
// ****************************************************************
template <class Param>
class FuncSet : public Func<Param> {

 protected:
  int         len;
  double*     valuesM;
  double*     valuesP;

 public:
  FuncSet(Param& _coordinate, int _len);
  FuncSet();
  virtual ~FuncSet();
  virtual void   init();
  virtual void   init(Param& _coordinate, int _len);
  virtual void   calcValueCenter();
  virtual void   calcValueCenter(Param& _coordinate);
  virtual void   calcValueSides();
  virtual void   calcValueSides(Param& _coordinate);
  virtual double valuePt();
  virtual double diff();
  virtual double ddiff();
  virtual void   summary();

};


// ****************************************************************
// *                       FUNCSETMULTIVAR                        *
// ****************************************************************
template <class Param>
class FuncSetMultivar : public FuncSet<Param> {
 protected:
  int      numVar;
  int      len_double;
  double*  _diffResult;

 public:
  FuncSetMultivar(Param& _coordinate, int _len);
  FuncSetMultivar();
  virtual ~FuncSetMultivar();
  virtual void   init();
  virtual void   init(Param& _coordinate, int _len);
  virtual void   calcValueSides();
  virtual void   calcValueSides(Param& _coordinate);
  virtual double diff();
  virtual double diff(int v);
  virtual double ddiff();
  virtual void   attachDiffResult(double* __diffResult);
  virtual void   summary();
};


// ****************************************************************
// *                           FUNCSET                            *
// ****************************************************************
/*
template <class Param>
class FuncDiff : public FuncSetMultivar<Param> {
 protected:
  SingleParticleFunc<Param>** dfunction;
  SingleParticleFunc<Param>** ddfunction;

 public:
  FuncDiff(Param& _coordinate, int _len);
  FuncDiff();
  virtual ~FuncDiff();
  virtual void   init();
  virtual void   init(Param& _coordinate, int _len);
  virtual void   calcValueCenter();
  virtual void   calcValueCenter(Param& _coordinate);
  virtual void   calcValueSides();
  virtual void   calcValueSides(Param& _coordinate);
  virtual double valuePt();
  virtual double valuePt(Param& _coordinate);
  virtual double diff();
  virtual double diff(Param& _coordinate);
  virtual double ddiff();
  virtual double ddiff(Param& _coordinate);

};
*/

#include "Func.cpp"

#endif
