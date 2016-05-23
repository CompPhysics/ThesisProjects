#include "Jastrow.h"

// ****************************************************************
// *                           JASTROW                            *
// ****************************************************************
//
// *************************** attach *****************************
void Jastrow::attach(Domain* _domain, fFunction* _f) {
  domain           = _domain;
  spinFactors      = domain->getSpinFactors();
  distance         = domain->getDistance();
  f                = _f;
  numParticles     = domain->getNumParticles();
  Nm1              = numParticles - 1;

  jastrowMatrix.resize(numParticles, numParticles);
  trialColumn.resize(numParticles);
  trialRow.resize(numParticles);

  trialDistanceColumn.resize(numParticles);
  trialDistanceRow.resize(numParticles);

  spinFlip         = domain->getSpinFlip();
  newFactor        = spinFactors->getNewFactor();
  oldFactor        = spinFactors->getOldFactor();
  otherDifference  = spinFactors->getOtherDifference();

  coors            = domain->getCoors();
  trialCoor        = domain->getTrialCoor();
}


// ************************* initialize ***************************
void Jastrow::initialize() {

  Array<double, 2> interElectronicDistances 
    = distance->getInterElectronicDistances();
  for(int i=0; i<numParticles; i++)
    for (int j=0; j<numParticles; j++)
      jastrowMatrix(i,j) = (*f)( interElectronicDistances(i,j), 
				 coors[i].r(), coors[j].r());

  spinFactors->calculateMatrix();
  jastrowian = 0;
  int k=0;
  for(int i=0; i<numParticles-1; i++) {
    for (int j=i+1; j<numParticles; j++)
      jastrowian += jastrowMatrix(i,j) * spinFactors->getMatrix()[k++];
  }
  setCurrentParticle(0);
  difference = 0;
}


// ********************** setToNextParticle ***********************
void Jastrow::setToNextParticle() {
  currentParticle++;
  if (currentParticle==numParticles) currentParticle=0;
  difference = 0;
}


// ********************* setCurrentParticle ***********************
void Jastrow::setCurrentParticle(int _currentParticle) {
  currentParticle = _currentParticle;
  difference = 0;
}


// ************************ suggestMove ***************************
void Jastrow::suggestMove() {
  if (*spinFlip) suggestMoveFlip();
  else suggestMoveNoFlip();
}


// ********************** suggestMoveNoFlip ***********************
void Jastrow::suggestMoveNoFlip() {
  trialDistanceColumn = distance->getTrialColumn();
  trialDistanceRow    = distance->getTrialRow();
  for ( int i=0; i < currentParticle; i++ ) {
    trialColumn(i) = trialRow(i) =
      (*f)( trialDistanceColumn(i), coors[i].r(), trialCoor->r());
  }
  for ( int i=currentParticle+1; i < numParticles; i++ ) {
    trialColumn(i) = trialRow(i) =
      (*f)( trialDistanceColumn(i), trialCoor->r(), coors[i].r());
  }
  trialColumn(currentParticle) = trialRow(currentParticle) = 0;
  calculateNoFlipDifference();
}


// ****************** calculateNoFlipDifference *******************
void Jastrow::calculateNoFlipDifference() {
  difference = 0;
  for (int i = 0; i<currentParticle; i++)
    difference += 
      newFactor[i] * ( trialColumn(i) - jastrowMatrix(i, currentParticle) );
  for (int i = currentParticle+1; i < numParticles; i++ )
    difference += 
      newFactor[i-1] * ( trialColumn(i) - jastrowMatrix(i, currentParticle) );
}


// *********************** suggestMoveFlip ************************
void Jastrow::suggestMoveFlip() {
  trialDistanceColumn = distance->getTrialColumn();
  trialDistanceRow    = distance->getTrialRow();
  otherParticle = domain->getOtherParticle();
  for ( int i=0; i < currentParticle; i++ ) {
    trialColumn(i) = trialRow(i) =
      (*f)(trialDistanceColumn(i),  coors[i].r(), trialCoor->r());
  }
  for ( int i=currentParticle+1; i < numParticles; i++ ) {
    trialColumn(i) = trialRow(i) =
      (*f)( trialDistanceColumn(i), trialCoor->r(), coors[i].r());
  }
  trialColumn(currentParticle) = trialRow(currentParticle) = 0;
  calculateFlipDifference();
}


// ******************* calculateFlipDifference ********************
void Jastrow::calculateFlipDifference() {
  difference = 0;
  for (int i = 0; i<currentParticle; i++)
    difference += 
      newFactor[i] * trialColumn(i) 
      - oldFactor[i]*jastrowMatrix(i, currentParticle)
      + otherDifference[i] * jastrowMatrix(i, otherParticle);
  for (int i = currentParticle+1; i < numParticles; i++ )
    difference += 
      newFactor[i-1] * trialColumn(i) 
      - oldFactor[i] * jastrowMatrix(i, currentParticle)
      + otherDifference[i] * jastrowMatrix(i, otherParticle);
}


// ************************* acceptMove ***************************
void Jastrow::acceptMove() {
  Range N(0,Nm1);
  jastrowMatrix( currentParticle, N )  = trialColumn( N );
  jastrowMatrix( N , currentParticle ) = trialRow( N );

  jastrowian+=difference;
}


// ************************** operator ****************************
double& Jastrow::operator()() {
  return jastrowian;
}

// ************************ getDifference *************************
double Jastrow::getDifference() {
  return difference;
}



// ************   CLASS JastrowPropose   ************

// ****************************************************************
// *                         JASTROWDIFF                          *
// ****************************************************************
//
// *************************** attach *****************************
void JastrowDiff::attach(Domain* _domain, fFunction* _f, 
			 int _differentiate) {
  domain           = _domain;
  int differentiate    = _differentiate;
  distanceDiff     = &(domain->getDistanceDiff()[differentiate]);
  spinFactors      = domain->getSpinFactors();
  f                = _f;
  numParticles     = domain->getNumParticles();
  numDimensions    = domain->getNumDimensions();
  Nm1              = numParticles - 1;

  jastrowMatrixPlus.resize(numParticles, numParticles);
  jastrowMatrixMinus.resize(numParticles, numParticles);
  trialColumnPlus.resize(numParticles);
  trialRowPlus.resize(numParticles);
  trialColumnMinus.resize(numParticles);
  trialRowMinus.resize(numParticles);

  trialDistanceColumn.resize(numParticles);
  trialDistanceRow.resize(numParticles);

  coors            = domain->getCoors();     // HER
  trialCoor        = domain->getTrialCoor(); // HER
  diffPlus         = differentiate;
  diffMinus        = differentiate + numDimensions;
}


// ************************* initialize ***************************
void JastrowDiff::initialize() {
  Array<double, 2> interElectronicDistances(numParticles,numParticles);
  interElectronicDistances  = distanceDiff->getInterElectronicDistances();
  for(int i=0; i<numParticles; i++)
    for (int j=0; j<numParticles; j++)
      jastrowMatrixPlus(i,j) 
	= (*f)( interElectronicDistances(i,j), 
		coors[i].r(diffPlus), coors[j].r(diffPlus));
  for(int i=0; i<numParticles; i++)
    for (int j=0; j<numParticles; j++) {
      jastrowMatrixMinus(i,j) 
	= (*f)( interElectronicDistances(j,i), 
		coors[i].r(diffMinus), coors[j].r(diffMinus));
    }
  setCurrentParticle(0);
}


// ********************** setToNextParticle ***********************
void JastrowDiff::setToNextParticle() {
  currentParticle++;
  if (currentParticle==numParticles) currentParticle=0;
}


// ********************* setCurrentParticle ***********************
void JastrowDiff::setCurrentParticle(int _currentParticle) {
  currentParticle = _currentParticle;
}


// ************************ suggestMove ***************************
void JastrowDiff::suggestMove() {
  trialDistanceColumn = distanceDiff->getTrialColumn();
  trialDistanceRow    = distanceDiff->getTrialRow();
  for ( int i=0; i < currentParticle; i++ ) {
    trialColumnPlus(i) = 
      (*f)( trialDistanceColumn(i), coors[i].r(diffPlus)
	    , trialCoor->r());
    trialRowPlus(i) =
      (*f)( trialDistanceRow(i), coors[i].r()
	    , trialCoor->r(diffPlus));

    trialColumnMinus(i) = 
      (*f)( trialDistanceRow(i), coors[i].r(diffMinus)
	    , trialCoor->r());
    trialRowMinus(i) =
      (*f)( trialDistanceColumn(i), coors[i].r()
	    , trialCoor->r(diffMinus));
  }
  for ( int i=currentParticle+1; i < numParticles; i++ ) {
    trialColumnPlus(i) = 
      (*f)( trialDistanceColumn(i), trialCoor->r()
	    , coors[i].r(diffPlus) );
    trialRowPlus(i) =
      (*f)( trialDistanceRow(i), trialCoor->r(diffPlus), 
	    coors[i].r() );

    trialColumnMinus(i) = 
      (*f)( trialDistanceRow(i), trialCoor->r(), 
	    coors[i].r(diffMinus) );
    trialRowMinus(i) =
      (*f)( trialDistanceColumn(i) , trialCoor->r(diffMinus), 
	    coors[i].r());
  }
  trialColumnPlus(currentParticle)  = trialRowPlus(currentParticle)  = 0;
  trialColumnMinus(currentParticle) = trialRowMinus(currentParticle) = 0;

}


// ************************* acceptMove ***************************
void JastrowDiff::acceptMove() {
  Range N(0,Nm1);
  jastrowMatrixPlus( currentParticle, N )   = trialColumnPlus( N );
  jastrowMatrixPlus( N , currentParticle )  = trialRowPlus( N );

  jastrowMatrixMinus( currentParticle, N )  = trialColumnMinus( N );
  jastrowMatrixMinus( N , currentParticle ) = trialRowMinus( N );
}
