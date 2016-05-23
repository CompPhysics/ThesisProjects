#ifndef VariationCPP_IS_INCLUDED
#define VariationCPP_IS_INCLUDED

#include "Variations.h"



// ****************************************************************
// *                          VARIATIONS                          *
// ****************************************************************
//
// ************************** Vatiations **************************
template <class SlaterDeterminant>
Variations<SlaterDeterminant>::Variations(Domain& _domain) : domain(_domain) { 
  numParticles         = domain().getNumParticles();
  numSlaterVariations  = domain().getNumSlaterVariations();
  numCorrelVariations  = domain().getNumCorrelVariations();
  allowSpinFlip        = domain().getAllowSpinFlip();
  outputFile           = domain().getOutputFile();
  varianceOptimization = domain().getVarianceOptimization();

  createSlaterDet();
  createCorrel();
  createLocalSurface();



  // !!! Temporary !!!
  distanceFileIndex = 0;
  distanceFile.open("distanceFile");
}


// *********************** createSlaterDet ************************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::createSlaterDet() { 
  slater = new SlaterDeterminant[numSlaterVariations];
  for (int i=0; i< numSlaterVariations; i++ ) 
    slater[i].init(&domain(), i);
  centralSlater = domain().getCentralSlater();
  slaterDeterminant = &(slater[centralSlater]);
}


// ************************* createCorrel *************************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::createCorrel() {
  correl = new Correl[numCorrelVariations];
  for (int i=0; i<numCorrelVariations; i++) {
    correl[i].attach(domain());
    correl[i].attachFFunction(&(domain().getF()[i]));
    correl[i].createJastrowAndJastrowDiff();
  }
  centralCorrel = domain().getCentralCorrel();
  correlation = &(correl[centralCorrel]);
}


// ********************** setToNextParticle ***********************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::setToNextParticle() { 
  for (int i=0; i<numCorrelVariations; i++)
    correl[i].setToNextParticle();
  for (int i=0; i< numSlaterVariations; i++ )
    slater[i].setToNextParticle();
}


// ********************** setCurrentParticle **********************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::setCurrentParticle(int _currentParticle) { 
  for (int i=0; i<numCorrelVariations; i++)
    correl[i].setCurrentParticle(_currentParticle);
  for (int i=0; i< numSlaterVariations; i++ )
    slater[i].setCurrentParticle(_currentParticle);
}


// ******************* initializeThermalization *******************
// Called prior to starting the thermalization procedure.
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::initializeThermalization() { 
  correlation().initializeThermalization();
  setCurrentParticle(0);
}


// ************************ initializeVMC *************************
// Called inbetween VMC and thermalization
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::initializeVMC() { 
  for (int i=0; i<numCorrelVariations; i++) {
    correl[i].initializeThermalization();
    correl[i].initializeVMC();
  }
  for (int i=0; i< numSlaterVariations; i++)
    slater[i].initDiff();

  for (int i=0; i<numVariations; i++ )
    localWaveFunctions[i].init();
  setCurrentParticle(0);
}


// ************************ initNewVmcRun *************************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::initNewVmcRun() { 
  for ( int i=0; i<numSlaterVariations; i++ ) {
    slater[i].initNewVmcRun();
  }
}


// ********************** createLocalSurface **********************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::createLocalSurface() {
  numVariations       = numSlaterVariations*numCorrelVariations;
  double* _centralSlaterDet   = &slaterDeterminant().getDet();
  double* _centralCorrelation = correlation().getCorrelationPtr();
  localWaveFunctions = new LocalWaveFunction[numVariations];
  for ( int i=0; i<numVariations; i++ ) {
    localWaveFunctions[i].initialize( &domain() );
    localWaveFunctions[i].attachCentralSlaterDet(_centralSlaterDet);
    localWaveFunctions[i].attachCentralCorrelation(_centralCorrelation);
  }

  int k=0;
  for ( int i=0; i<numSlaterVariations; i++ ) {
    for ( int j=0; j<numCorrelVariations; j++ ) {
      localWaveFunctions[k]
	.attachSlater(&slater[i].getDet(),
		      slater[i].getDiffRatiosPtr(),
		      &slater[i].getDDiffRatio(),
		      domain().getAlphaParam(i));
      localWaveFunctions[k]
	.attachCorrelation(correl[j].getCorrelationPtr(),
			   correl[j].getGradRatio(),
			   correl[j].getLaplaceRatio(),
			   domain().getBetaParam(j));
      k++;
    }
  }
}


// ************************* suggestMove **************************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::suggestMove() {
  correlation().suggestMove();
  if (allowSpinFlip) slaterDeterminant().suggestMoveFlip();
  else slaterDeterminant().suggestMove();
}


// ********************* acceptThermalization *********************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::acceptThermalization() {
  correlation().acceptThermalizedMove();
  if (allowSpinFlip) {
    for (int i=0; i<numSlaterVariations; i++ )
      if (i!=centralSlater)
	slater[i].suggestMoveFlip();
    for (int i=0; i<numSlaterVariations; i++ )
      slater[i].acceptMoveFlip();
  }
  else {
    for (int i=0; i< numSlaterVariations; i++ )
      if (i!=centralSlater)
	slater[i].suggestMove();
    for (int i=0; i< numSlaterVariations; i++ )
      slater[i].acceptMove();
  }
}


// ********************* rejectThermalization *********************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::rejectThermalization() {
  correlation().rejectThermalizedMove();
  if (allowSpinFlip) {
    for (int i=0; i<numSlaterVariations; i++ )
      if (i!=centralSlater)
	slater[i].suggestMoveFlip();
    for (int i=0; i< numSlaterVariations; i++ )
      slater[i].rejectMoveFlip();
  }
  else {
    for (int i=0; i< numSlaterVariations; i++ )
      if (i!=centralSlater)
	slater[i].suggestMove();
    for (int i=0; i< numSlaterVariations; i++ )
      slater[i].rejectMove();
  }
}


// ************************** acceptMove **************************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::acceptMove() {
  for (int i=0; i<numCorrelVariations; i++)
    if (i!=centralCorrel)
      correl[i].suggestMove();
  for (int i=0; i<numCorrelVariations; i++) {
    correl[i].acceptMove();
  }
  if (allowSpinFlip) {
    for (int i=0; i<numSlaterVariations; i++ )
      if (i!=centralSlater)
	slater[i].suggestMoveFlip();
    for (int i=0; i<numSlaterVariations; i++ )
      slater[i].acceptMoveFlip();
  }
  else {
    for (int i=0; i< numSlaterVariations; i++ )
      if (i!=centralSlater)
	slater[i].suggestMove();
    for (int i=0; i< numSlaterVariations; i++ )
      slater[i].acceptMove();
  }
}


// ************************** rejectMove **************************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::rejectMove() {

  for (int i=0; i<numCorrelVariations; i++) {
    correl[i].rejectMove();
  }
  if (allowSpinFlip) {
    for (int i=0; i<numSlaterVariations; i++ )
      if (i!=centralSlater)
	slater[i].suggestMoveFlip();
    for (int i=0; i< numSlaterVariations; i++ )
      slater[i].rejectMoveFlip();
  }
  else {
    for (int i=0; i< numSlaterVariations; i++ )
      if (i!=centralSlater)
	slater[i].suggestMove();
    for (int i=0; i< numSlaterVariations; i++ )
      slater[i].rejectMove();
  }
}


// **************************** sample ****************************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::sample() {
  for (int i=0; i<numCorrelVariations; i++) {
    correl[i].calculateCorrelation();
    correl[i].calculateGradAndLaplacianRatios();
  }

  for (int i=0; i< numSlaterVariations; i++ ) {
    slater[i].calcDDiffRatio();
    slater[i].calcDiffRatios();
  }
  double potentialEnergy = domain().getNucleusElectronPotential() 
    + domain().getInterElectronicPot();



  

  for (int i=0; i<numVariations; i++) {
    localWaveFunctions[i].sample(potentialEnergy);
  }

  // !!! Temporary !!!
  /*
  distanceFile << distanceFileIndex++ << " ";
  for (int i=0; i<numParticles; i++) {
    distanceFile << (domain().getCoors()[i]).r() << " ";
  }
  distanceFile << " " << localWaveFunctions[0].getLocalEnergy() << endl;
  */
}


// *************************** summary ****************************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::summary() {
  for (int i=0; i<numVariations; i++)
    localWaveFunctions[i].summary();
}


// ************************ summaryLowest *************************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::summaryLowest() {
  if (varianceOptimization)  *outputFile << "Lowest Variance for: ";
  else  *outputFile << "Lowest Energy for: ";
  localWaveFunctions[findLowestWaveFunction()].summary();
}


// ******************** findLowestWaveFunction ********************
// Returns the integer value of the lovalWaveFunction with the lowest 
// energy or variance
template <class SlaterDeterminant>
int Variations<SlaterDeterminant>::findLowestWaveFunction() {
  if (varianceOptimization) {
    int lowest=-1;
    double minVariance = 999;
    for (int i=0; i<numVariations; i++) {
      double variance = localWaveFunctions[i].getVariance();
      if ( variance < minVariance) {
	minVariance = variance;
	lowest = i;
      }
    }
    if (lowest ==-1) {
      *outputFile << "Error finding variance minimum in class Variation" 
		  << endl << "Minimum variance >= 999!" << endl;
      return -1;
    }
    else return lowest;
  }
  else {
    int lowest=-1;
    double minEnergy = 0;
    for (int i=0; i<numVariations; i++) {
      double energy = localWaveFunctions[i].getAverageEnergy();
      if ( energy < minEnergy) {
	minEnergy = energy;
	lowest = i;
      }
    }
    if (lowest ==-1) {
      *outputFile << "Error finding energy minimum in class Variation" 
		  << endl << "Minimum energy >= 0!" << endl;
      return -1;
    }
    else return lowest;
  }
}


// *********************** getLowestEnergy ************************
template <class SlaterDeterminant>
double Variations<SlaterDeterminant>::getLowestEnergy() {
  return localWaveFunctions[findLowestWaveFunction()].getAverageEnergy();
}


// ********************** getLowestVariance ***********************
template <class SlaterDeterminant>
double Variations<SlaterDeterminant>::getLowestVariance() {
  return localWaveFunctions[findLowestWaveFunction()].getVariance();
}


// ********************* getLowestAlphaParams *********************
template <class SlaterDeterminant>
double* Variations<SlaterDeterminant>::getLowestAlphaParams() {
  return localWaveFunctions[findLowestWaveFunction()].getAlphaParams();
}


// ********************* getLowestBetaParams **********************
template <class SlaterDeterminant>
double* Variations<SlaterDeterminant>::getLowestBetaParams() {
  return localWaveFunctions[findLowestWaveFunction()].getBetaParams();
}


// ********************** setReferenceEnergy **********************
template <class SlaterDeterminant>
void Variations<SlaterDeterminant>::setReferenceEnergy(double E) {
  for (int i=0; i<numVariations; i++)
    localWaveFunctions[i].setReferenceEnergy(E);
}

#endif
