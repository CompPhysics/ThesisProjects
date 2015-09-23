#ifndef Walker_IS_INCLUDED
#define Walker_IS_INCLUDED

#include <iostream>
#include <cmath>
#include "../Random/Random.h"
#include "../Domain/Domain.h"
#include "../SlaterDet/SlaterDet.h"
#include "../Ref/Ref.h"
#include "../Correlation/Correlation.h"
#include "../Variations/Variations.h"

// ****************************************************************
// *                           WALKER                             *
// ****************************************************************
template <class SlaterDeterminant>
class Walker {
 protected:

  int             accepted;         // The number of accepted steps.
  double          ratio;            // The product of slaterRatio 
                                    // and correlationRatio.
  RefFund<double> slaterRatio;
  RefFund<double> correlationRatio; 

  Ref<Random2>                        random;
  Ref<Domain>                         domain;
  Ref<SlaterDeterminant>              slaterDeterminant;
  Ref<Correl>                         correlation;
  Ref<Variations<SlaterDeterminant> > variations;

 public:
  Walker(Domain& _domain,  
	 SlaterDeterminant& _slaterDeterminant, 
	 Random2& _random, Correl& _correlation, 
	 Variations<SlaterDeterminant>& _variations);

  int                doThermalizationStep();
  int                doRandomStep();
  void               resetAcceptance();
  int                getAcceptance();
  SlaterDeterminant* getSlaterDeterminantPtr() {return &(slaterDeterminant());}
};

#include "Walker.cpp"

#endif
