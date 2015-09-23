#ifndef SlaterDet_IS_INCLUDED
#define SlaterDet_IS_INCLUDED

#include <iostream>
#include "../SlaterMatrix/SlaterMatrix.h"
#include "../Coor/Coor.h"
#include "../Domain/Domain.h"
#include "../Func/Func.h"
#include "../SingleParticleFuncs/SingleParticleFuncs.h"
#include "../Random/Random.h"
#include "../Ref/Ref.h"


// ****************************************************************
// *                         SLATERDET                            *
// ****************************************************************
template <class FuncUp, class FuncDown>
class SlaterDet {
   
 protected:
   
  SlaterMatrix_Det *smUp, *smDown;
   
  int          numParticles;
  int          numDimensions;
  int          dimUp, dimDown;
  double       *newValuesColumnUp, *newValuesColumnDown;
  int          *columnIndexOfParticle;
  
  Domain       *domain;
  CoorSpinDiff *coors;
  Ref<CoorSpinDiff> trialCoor;
  FuncUp       *funcsUp;
  FuncDown     *funcsDown;
  Ref<Random2> random;
  
  Ref<double>  ratioUp, ratioDown;
  double       ratio;
  double       dDiffRatio;
  double       *diffRatios, *_diffRatios;
  double       det;
  Ref<double>  detUp, detDown;
  double       trialDet;
  Ref<double>  trialDetUp, trialDetDown;
  
  int          *spinFlip;
  int          currentParticle;
  int          otherParticle;
  int          spin;
  int          currentColumn;
  int          otherColumn;
  int          moveAcceptance, flipAcceptance, moveFlipAcceptance;
  
#ifdef _DEBUG_
  int          moved, flipped, valueSidesCalculated;
#endif
  
 public:
  
  SlaterDet();
  ~SlaterDet();

  void      init(Domain* _domain, int alphaVar);
  void      initNewVmcRun();
  void      setToNextParticle();
  void      setCurrentParticle(int _currentParticle);
  
  void      suggestMove();
  void      suggestFlip();
  void      suggestMoveFlip();
  
  void      acceptMove();
  void      acceptFlip();
  void      acceptMoveFlip();
  
  void      rejectMove();
  void      rejectFlip();
  void      rejectMoveFlip();
  
  void      calcDiffRatios();
  void      calcDDiffRatio();
  
  void      setSpinFlip(int _sf)    {spinFlip=_sf;}
  int&      getSpinFlip()           {return spinFlip;}
  int&      getCurrentParticle()    {return currentParticle;}
  int&      getOtherParticle()      {return otherParticle;}
  double&   getDet()                {return det;}
  double&   getTrialDet()           {return trialDet;}
  double&   getRatio()              {return ratio;}
  double*   getDiffRatiosPtr()      {return diffRatios;}
  double&   getDDiffRatio()         {return dDiffRatio;}
  
  FuncUp*   getFuncsUpPtr()         {return funcsUp;}
  FuncDown* getFuncsDownPtr()       {return funcsDown;}
  int       getMoveAcceptance()     {return moveAcceptance;}
  int       getFlipAcceptance()     {return flipAcceptance;}
  int       getMoveFlipAcceptance() {return moveFlipAcceptance;}
  void      resetAcceptances();
  
  void      summary();
  void      initDiff();
  
 protected:
  
  
#ifdef _DEBUG_
  int       checkIfValueSidesCalculated();
#endif
  
  
};


#include "SlaterDet.cpp"

#endif
