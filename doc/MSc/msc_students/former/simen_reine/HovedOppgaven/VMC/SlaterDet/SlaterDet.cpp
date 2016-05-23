#ifndef SlaterDetCPP_IS_INCLUDED
#define SlaterDetCPP_IS_INCLUDED

#include "SlaterDet.h"

// ****************************************************************
// *                         SLATERDET                            *
// ****************************************************************
//
// ************************* SlaterDet ****************************
template <class FuncUp, class FuncDown>
SlaterDet<FuncUp, FuncDown>::SlaterDet() {}


// ************************* ~SlaterDet ***************************
template <class FuncUp, class FuncDown>
SlaterDet<FuncUp, FuncDown>::~SlaterDet() 
{
  delete[] columnIndexOfParticle;
  delete[] newValuesColumnUp;
  delete[] newValuesColumnDown;
  delete[] diffRatios;
  delete[] funcsUp;
  delete[] funcsDown;
  delete[] smUp;
  delete[] smDown;
}


// **************************** init ******************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::init(Domain* _domain, int alphaVar) {
 
  domain          = _domain; 
  numParticles    = domain->getNumParticles();
  numDimensions   = domain->getNumDimensions();
  coors           = domain->getCoors();
  trialCoor       = domain->getTrialCoor();
  spinFlip        = domain->getSpinFlip();
 
  //Calculate the number of particles with spin up and down
  //and generate initial particle index to column index transform array
  columnIndexOfParticle = new int[numParticles];
  dimUp   = 0;
  dimDown = 0;
  for (int i = 0; i < numParticles; i++)
    if (coors[i].spin() > 0)
      columnIndexOfParticle[i] = dimUp++;
    else  // if (coors[i].spin() == -1)
      columnIndexOfParticle[i] = dimDown++;

  //Allocate newValueColumn arrays
  newValuesColumnUp   = new double[dimUp];
  newValuesColumnDown = new double[dimDown];
  
  //Allocate array for first derivative ratios
  diffRatios = new double[numParticles*numDimensions + 1];
  
  //Generate function objects
  funcsUp   = new FuncUp[dimUp];
  funcsDown = new FuncDown[dimDown];
  for (int i = 0; i < dimUp; i++) {
    funcsUp[i].init(trialCoor(), dimUp, domain, alphaVar);
    funcsUp[i].attachResult(newValuesColumnUp);
    funcsUp[i].attachDdiffResult(newValuesColumnUp);
  }
  for (int i = 0; i < dimDown; i++) {
    funcsDown[i].init(trialCoor(), dimDown, domain, alphaVar);
    funcsDown[i].attachResult(newValuesColumnDown);
    funcsDown[i].attachDdiffResult(newValuesColumnDown);
  }
  
  //Generate SlaterMatrix objects
  smUp   = new SlaterMatrix_Det[1];
  smDown = new SlaterMatrix_Det[1];
  smUp->redim(dimUp);
  smUp->setPtrToNewColumn(newValuesColumnUp);
  smDown->redim(dimDown);
  smDown->setPtrToNewColumn(newValuesColumnDown);
  ratioUp      = smUp->getRatio();
  ratioDown    = smDown->getRatio();
  detUp        = smUp->getDet();
  detDown      = smDown->getDet();
  trialDetUp   = smUp->getTrialDet();
  trialDetDown = smDown->getTrialDet();
  
  //Initialize SlaterMatrix objects
  for (int i = 0; i < dimUp; i++) {
    funcsUp[i].calcValueCenter(coors[i]);
    funcsUp[i].valuePt();
    smUp->calcRatio();
    smUp->importNewColumn();
    smUp->setToNextColumn();
  }
  for (int i = 0; i < dimDown; i++) {
    funcsDown[i].calcValueCenter(coors[i+dimUp]);
    funcsDown[i].valuePt();
    smDown->calcRatio();
    smDown->importNewColumn();
    smDown->setToNextColumn();
  }
  
#ifdef _DEBUG_
  moved = flipped = valueSidesCalculated = 0;
#endif
  
  setCurrentParticle(0);
  resetAcceptances();
  initDiff();
}


// *********************** initNewVmcRun **************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::initNewVmcRun() {

  //Calculate the number of particles with spin up and down
  //and generate initial particle index to column index transform array
  dimUp   = 0;
  dimDown = 0;
  for (int i = 0; i < numParticles; i++)
    if (coors[i].spin() > 0)
      columnIndexOfParticle[i] = dimUp++;
    else   // if (coors[i].spin() == -1)
      columnIndexOfParticle[i] = dimDown++;
  
  //Initialize SlaterMatrix objects
  smUp->init();
  for (int i = 0; i < dimUp; i++) {
    funcsUp[i].calcValueCenter(coors[i]);
    funcsUp[i].valuePt();
    smUp->calcRatio();
    smUp->importNewColumn();
    smUp->setToNextColumn();
  }

  smDown->init();
  for (int i = 0; i < dimDown; i++) {
    funcsDown[i].calcValueCenter(coors[i+dimUp]);
    funcsDown[i].valuePt();
    smDown->calcRatio();
    smDown->importNewColumn();
    smDown->setToNextColumn();
  }
  setCurrentParticle(0);
  resetAcceptances();
  initDiff();
}
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::suggestMove() {
#ifdef _DEBUG_
  if (moved == 1 || flipped == 1) {
    cerr << "\nCannot move. Change already suggested.\nAccept or reject first.\n";
    return;
  }
#endif
  
  currentColumn = columnIndexOfParticle[currentParticle];
  
  if (spin > 0) {
    funcsUp[currentColumn].calcValueCenter();
    funcsUp[currentColumn].valuePt();
    smUp->setCurrentColumn(currentColumn);
    smUp->calcRatio();
    ratio    = ratioUp();
    trialDet = trialDetUp();
  }
  else {
    funcsDown[currentColumn].calcValueCenter();
    funcsDown[currentColumn].valuePt();
    smDown->setCurrentColumn(currentColumn);
    smDown->calcRatio();
    ratio    = ratioDown();
    trialDet = trialDetDown();
  }
  
#ifdef _DEBUG_
  moved = 1;
  valueSidesCalculated = 0;
#endif
}


// ************************ suggestFlip ***************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::suggestFlip() {
#ifdef _DEBUG_
  if (moved == 1 || flipped == 1) {
    cerr << "\nCannot flip. Change already suggested.\nAccept or reject first.\n";
    return;
  }
#endif
  
  //Find a particle with a spin opposite of the current particle
  while ( coors[otherParticle = (int)(random().getNum() * numParticles)].spin() == spin );
  
  currentColumn = columnIndexOfParticle[currentParticle];
  otherColumn   = columnIndexOfParticle[otherParticle];
  
  if (spin > 0) {
    funcsUp[currentColumn].calcValueCenter(coors[otherParticle]);
    funcsDown[otherColumn].calcValueCenter(coors[currentParticle]);
    funcsUp[currentColumn].valuePt();
    funcsDown[otherColumn].valuePt();
    smUp->setCurrentColumn(currentColumn);
    smDown->setCurrentColumn(otherColumn);
  }
  else {
    funcsDown[currentColumn].calcValueCenter(coors[otherParticle]);
    funcsUp[otherColumn].calcValueCenter(coors[currentParticle]);
    funcsDown[currentColumn].valuePt();
    funcsUp[otherColumn].valuePt();
    smDown->setCurrentColumn(currentColumn);
    smUp->setCurrentColumn(otherColumn);
  }
  
  smUp->calcRatio();
  smDown->calcRatio();
  
  ratio    = ratioUp()    * ratioDown();
  trialDet = trialDetUp() * trialDetDown();
  
#ifdef _DEBUG_
  flipped = 1;
  valueSidesCalculated = 0;
#endif
}


// ********************** suggestMoveFlip *************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::suggestMoveFlip() {
#ifdef _DEBUG_
  if (moved == 1 || flipped == 1) {
    cerr << "\nCannot move-flip. Change already suggested.\nAccept or reject first.\n";
    return;
  }
#endif
  
  //Find a particle with a spin opposite of the current particle
  //while ( coors[otherParticle = (int)(random().getNum() * numParticles)].spin() == spin );

  // Chose otherParticle randomly, and if the spins differ flip spin
  // otherParticle = (int)(random().getNum() * numParticles);
  // spinFlip = (coors[otherParticle].spin() == spin);
  if (!*spinFlip) suggestMove();
  else {
    otherParticle = domain->getOtherParticle();
    currentColumn = columnIndexOfParticle[currentParticle];
    otherColumn   = columnIndexOfParticle[otherParticle];
  
    if (spin > 0) {
      funcsUp[currentColumn].calcValueCenter(coors[otherParticle]);
      funcsDown[otherColumn].calcValueCenter();
      funcsUp[currentColumn].valuePt();
      funcsDown[otherColumn].valuePt();
      smUp->setCurrentColumn(currentColumn);
      smDown->setCurrentColumn(otherColumn);
    }
    else {
      funcsDown[currentColumn].calcValueCenter(coors[otherParticle]);
      funcsUp[otherColumn].calcValueCenter();
      funcsDown[currentColumn].valuePt();
      funcsUp[otherColumn].valuePt();
      smDown->setCurrentColumn(currentColumn);
      smUp->setCurrentColumn(otherColumn);
    }
  
    smUp->calcRatio();
    smDown->calcRatio();
  
    ratio    = ratioUp()    * ratioDown();
    trialDet = trialDetUp() * trialDetDown();
  
#ifdef _DEBUG_
  moved = flipped = 1;
  valueSidesCalculated = 0;
#endif
  }

}


// ************************ acceptMove ****************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::acceptMove() {
#ifdef _DEBUG_
  if (flipped == 1 || moved == 0) {
    cerr << "\nCannot accept move. A flip or a move-flip (or none) was suggested.\n";
    return;
  }
#endif

  if (spin > 0) {
    smUp->importNewColumn();
    det = detUp();
    //Calculate and store differentiated values of single orbital
    //functions. These values will be needed to calculate the first and
    //second derivatives later on.
    funcsUp[currentColumn].calcValueSides();
  }
  else {
    smDown->importNewColumn();
    det = detDown();
    //Calculate and store differentiated values of single orbital
    //functions. These values will be needed to calculate the first
    //and second derivatives later on.
    funcsDown[currentColumn].calcValueSides();
  }
  
  moveAcceptance++;
  
#ifdef _DEBUG_
  moved = valueSidesCalculated = 0;
#endif
}


// ************************ acceptFlip ****************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::acceptFlip() {
#ifdef _DEBUG_
  if (moved == 1 || flipped == 0) {
    cerr << "\nCannot accept flip. A move or a move-flip (or none) was suggested.\n";
    return;
  }
#endif

  if (spin > 0) {
    funcsUp[currentColumn].calcValueSides(coors[otherParticle]);
    funcsDown[otherColumn].calcValueSides(coors[currentParticle]);
  }
  else {
    funcsDown[currentColumn].calcValueSides(coors[otherParticle]);
    funcsUp[otherColumn].calcValueSides(coors[currentParticle]);
  }
  
  //Let the slater matrices import their new columns
  smUp->importNewColumn();
  smDown->importNewColumn();
  det = detUp() * detDown();
  
  //Flip the spin coordinates of the particles involved
  coors[currentParticle].flipSpin();
  coors[otherParticle].flipSpin();
  
  //Interchange indices in particle to column transform array
  columnIndexOfParticle[currentParticle] = otherColumn;
  columnIndexOfParticle[otherParticle]   = currentColumn;
  
  flipAcceptance++;
  
#ifdef _DEBUG_
  flipped = valueSidesCalculated = 0;
#endif
}


// ********************** acceptMoveFlip **************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::acceptMoveFlip() {
  if (!*spinFlip) acceptMove();
  else {
#ifdef _DEBUG_
    if (moved == 0 || flipped == 0) {
      cerr << "\nCannot accept move-flip. A move or a flip (or none) was suggested.\n";
      return;
    }
#endif
    if (spin > 0) {
      funcsUp[currentColumn].calcValueSides(coors[otherParticle]);
      funcsDown[otherColumn].calcValueSides();
    }
    else {
      funcsDown[currentColumn].calcValueSides(coors[otherParticle]);
      funcsUp[otherColumn].calcValueSides();
    }
  
    //Let the slater matrices import their new columns
    smUp->importNewColumn();
    smDown->importNewColumn();
    det = detUp() * detDown();
    
    //Flip the spin coordinates of the particles involved
    coors[currentParticle].flipSpin();
    coors[otherParticle].flipSpin();
    
    //Interchange indices in particle to column transform array
    columnIndexOfParticle[currentParticle] = otherColumn;
    columnIndexOfParticle[otherParticle]   = currentColumn;
  }
  
  moveFlipAcceptance++;
  
#ifdef _DEBUG_
  moved = flipped = valueSidesCalculated = 0;
#endif

}


// ************************ rejectMove ****************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::rejectMove() {
#ifdef _DEBUG_
  if (flipped == 1 || moved == 0) {
    cerr << "\nCannot reject move. A flip or a move-flip (or none) was suggested.\n";
    return;
  }
  moved = 0;
#endif
  
  if (spin > 0) {
    //Recalculate and store values of single orbital functions. These
    //values will be needed to calculate the first and second
    //derivatives later on.
    funcsUp[currentColumn].calcValueCenter(coors[currentParticle]);
  }
  else {
    //Calculate and store values of single orbital functions. These
    //values will be needed to calculate the first and second
    //derivatives later on.
    funcsDown[currentColumn].calcValueCenter(coors[currentParticle]);
  }
  
}


// ************************ rejectFlip ****************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::rejectFlip() {
#ifdef _DEBUG_
  if (moved == 1 || flipped == 0) {
    cerr << "\nCannot reject flip. A move or a move-flip (or none) was suggested.\n";
    return;
  }
  flipped = 0;
#endif
  
  if (spin > 0) {
    funcsUp[currentColumn].calcValueCenter(coors[currentParticle]);
    funcsDown[otherColumn].calcValueCenter(coors[otherParticle]);
  }
  else {
    funcsDown[currentColumn].calcValueCenter(coors[currentParticle]);
    funcsUp[otherColumn].calcValueCenter(coors[otherParticle]);
  }
  
}


// ********************** rejectMoveFlip **************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::rejectMoveFlip() {

  if (!*spinFlip) rejectMove();
  else {
#ifdef _DEBUG_
    if (moved == 0 || flipped == 0) {
      cerr << "\nCannot reject move-flip. A move or a flip (or none) was suggested.\n";
      return;
    }
    moved = flipped = 0;
#endif
  
    if (spin > 0) {
      funcsUp[currentColumn].calcValueCenter(coors[currentParticle]);
      funcsDown[otherColumn].calcValueCenter(coors[otherParticle]);
    }
    else {
      funcsDown[currentColumn].calcValueCenter(coors[currentParticle]);
      funcsUp[otherColumn].calcValueCenter(coors[otherParticle]);
    }
  }
}


// ********************* setToNextParticle ************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::setToNextParticle() {
#ifdef _DEBUG_
  if (moved == 1 || flipped == 1) {
    cerr << "\nA move, flip or move-flip was suggested.\nAccept or reject before setting to next particle.\n";
    return;
  }
#endif
  
  if (++currentParticle == numParticles)
    currentParticle = 0;
  
  spin = coors[currentParticle].spin();
}


// ********************* setCurrentParticle ***********************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::setCurrentParticle(int _currentParticle) {
#ifdef _DEBUG_
  if (moved == 1 || flipped == 1) {
    cerr << "\nA move, flip or move-flip was suggested.\nAccept or reject before setting particle index.\n";
    return;
  }
#endif
  
  currentParticle = _currentParticle;
  spin = coors[currentParticle].spin();
}


// ************************** initDiff ****************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::initDiff() {
  for (int i = 0; i < numParticles; i++)
    if (coors[i].spin() > 0) {
      funcsUp[columnIndexOfParticle[i]].calcValueCenter(coors[i]);
      funcsUp[columnIndexOfParticle[i]].calcValueSides(coors[i]);
    }
    else {
      funcsDown[columnIndexOfParticle[i]].calcValueCenter(coors[i]);
      funcsDown[columnIndexOfParticle[i]].calcValueSides(coors[i]);
   }

#ifdef _DEBUG_
  valueSidesCalculated = 1;
#endif
}


// *********************** calcDiffRatios *************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::calcDiffRatios() {
#ifdef _DEBUG_
  if (checkIfValueSidesCalculated() == 0) {
    cerr << "\nCannot calculate diffRatios.\n";
    exit(1);
  }
#endif
  
  _diffRatios = diffRatios;
  for (int i = 0; i < numParticles; i++)
    if (coors[i].spin() > 0) {
      smUp->setCurrentColumn(columnIndexOfParticle[i]);
      for (int j = 0; j < numDimensions; j++) {
        funcsUp[columnIndexOfParticle[i]].diff(j);
        smUp->calcRatio();
        *_diffRatios++ = ratioUp();
      }
    }
    else {
      smDown->setCurrentColumn(columnIndexOfParticle[i]);
      for (int j = 0; j < numDimensions; j++) {
        funcsDown[columnIndexOfParticle[i]].diff(j);
        smDown->calcRatio();
        (*_diffRatios++) = ratioDown();
      }
    }
}


// *********************** calcDDiffRatio *************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::calcDDiffRatio() {
#ifdef _DEBUG_
  if (checkIfValueSidesCalculated() == 0) {
    cerr << "\nCannot calculate dDiffRatio.\n";
    exit(1);
  }
#endif

  dDiffRatio = 0;
  for (int i = 0; i < dimUp; i++) {
    smUp->setCurrentColumn(i);
    funcsUp[i].ddiff();
    smUp->calcRatio();
    dDiffRatio += ratioUp();
  }
  for (int i = 0; i < dimDown; i++) {
    smDown->setCurrentColumn(i);
    funcsDown[i].ddiff();
    smDown->calcRatio();
    dDiffRatio += ratioDown();
  }
}


// ************************** summary *****************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::summary() {
  cerr << "\n\n--------------------------------------------------";
  cerr << "\n SUMMARY OF SlaterDet INSTANCE\n\n";
  
  smUp->summary();
  smDown->summary();
}


// ********************** resetAcceptances ************************
template <class FuncUp, class FuncDown>
void SlaterDet<FuncUp, FuncDown>::resetAcceptances() {
  moveAcceptance = flipAcceptance = moveFlipAcceptance = 0;
}


#ifdef _DEBUG_
template <class FuncUp, class FuncDown>
int SlaterDet<FuncUp, FuncDown>::checkIfValueSidesCalculated() {
  /*
  if (valueSidesCalculated == 0) {
    cerr << "\nvalueSides not yet calculated.";
    return 0;
  }
  */
  return 1;
}
#endif

#endif
