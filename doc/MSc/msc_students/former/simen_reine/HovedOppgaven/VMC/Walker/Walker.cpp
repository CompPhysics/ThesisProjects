#ifndef WalkerCPP_IS_INCLUDED
#define WalkerCPP_IS_INCLUDED

#include "Walker.h"


template <class SlaterDeterminant>
Walker<SlaterDeterminant>::Walker(Domain& _domain, SlaterDeterminant& _slaterDeterminant, Random2& _random, Correl& _correlation, Variations<SlaterDeterminant>& _variations) {
  random            = _random;
  domain            = _domain;
  slaterDeterminant = _slaterDeterminant;
  correlation       = _correlation;
  slaterRatio       = &slaterDeterminant().getRatio();
  correlationRatio  = correlation().getRatioPtr();
  variations        = _variations;
}


template <class SlaterDeterminant>
int Walker<SlaterDeterminant>::doThermalizationStep() {
  domain().suggestMove();
  //domain().proposeFlip();
  variations().suggestMove();
  ratio = slaterRatio()*correlationRatio();

  // THE METROPOLIS TEST!
  if ( (ratio*ratio) > random().getNum()) {
    domain().acceptThermalizedMove();
    variations().acceptThermalization();
    accepted = 1;
  }
  else {
    domain().rejectThermalizedMove();
    variations().rejectThermalization();
    accepted = 0;
  }
  variations().setToNextParticle();
  domain().setToNextParticle();
  return accepted;
}

template <class SlaterDeterminant>
int Walker<SlaterDeterminant>::doRandomStep() {
  domain().suggestMove();
  //domain().proposeFlip();
  variations().suggestMove();
  ratio = slaterRatio()*correlationRatio();
  // THE METROPOLIS TEST!
  if ( (ratio*ratio) > random().getNum()) {
    domain().acceptMove();
    variations().acceptMove();
    accepted = 1;
  }
  else {
    //domain().rejectTrialPosition();
    domain().rejectMove();
    variations().rejectMove();
    accepted = 0;
  }  
  variations().setToNextParticle();
  domain().setToNextParticle();
  return accepted;
}


template <class SlaterDeterminant>
void Walker<SlaterDeterminant>::resetAcceptance() {
  slaterDeterminant().resetAcceptances();
}


template <class SlaterDeterminant>
int Walker<SlaterDeterminant>::getAcceptance() {
  return slaterDeterminant().getMoveAcceptance();
}

#endif
