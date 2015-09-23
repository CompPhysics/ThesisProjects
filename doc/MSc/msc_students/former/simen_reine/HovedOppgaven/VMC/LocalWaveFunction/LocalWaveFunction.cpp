#include "LocalWaveFunction.h"

// ****************************************************************
// *                      LOCALWAVEFUNCTION                       *
// ****************************************************************
//
// Keeps track of one local variation
//
// **************************** init ******************************
void LocalWaveFunction::init() {
  numSamples                = 0;
  accumulatedEnergy         = 0;
  accumulatedEnergySquared  = 0;
  accumulatedNorm           = 0;
  accumulatedVariance       = 0;
}


// ************************* initialize ***************************
void LocalWaveFunction::initialize(Domain* domain) {
  numParticles              = domain->getNumParticles();
  numDimensions             = domain->getNumDimensions();
  prodPartDim               = numParticles*numDimensions;
  numAlpha                  = domain->getNumAlpha();
  numBeta                   = domain->getNumBeta();
  outputFile                = domain->getOutputFile();

  setWeightToUnity          = domain->getSetWeightToUnity();
  varianceOptimization      = domain->getVarianceOptimization();
  referenceEnergy           = domain->getReferenceEnergy();
  deltaE                    = domain->getDeltaE();

  init();
  deleteBlocks();
  maxNumBlocks              = 15;
  allocateBlocks();

  addBlock(2);
  addBlock(5);
  addBlock(10);
  addBlock(20);
  addBlock(25);

  addBlock(50);
  addBlock(100);
  addBlock(200);
  addBlock(500);
  addBlock(1000);

  addBlock(2000);
  addBlock(5000);
  addBlock(10000);
  addBlock(20000);
  addBlock(50000);
}


// *********************** allocateBlocks *************************
void LocalWaveFunction::allocateBlocks() {
  blockSize                     = new int[maxNumBlocks];
  currentBlockNum               = new int[maxNumBlocks];
  numSamplesBlock               = new int[maxNumBlocks];
  accumulatedEnergyBlock        = new double[maxNumBlocks];

  accumulatedEnergySquaredBlock = new double[maxNumBlocks];
  accumulatedNormBlock          = new double[maxNumBlocks];
  tempEnergyBlock               = new double[maxNumBlocks];
  tempNormBlock                 = new double[maxNumBlocks];
}


// ******************** attachCentralSlaterDet ********************
// Central == the origin of the local variation
// centralSlaterDet = pointer to the value of the central SD
void LocalWaveFunction::attachCentralSlaterDet(double* _centralSlater) {
  centralSlater   = _centralSlater;
}


// ****************** attachCentralCorrelation ********************
// Central == the origin of the local variation
// centralCorrelation = pointer to the value of the central 
// correlation
void LocalWaveFunction::attachCentralCorrelation(double* _centralCorrelation) {
  centralCorrelation = _centralCorrelation;
}


// ************************ attachSlater **************************
// These are the values of the (varied) SD
void LocalWaveFunction::attachSlater(double* _localSlater, 
				     double* _localDiffSlater,
				     double* _localDDiffSlater, 
				     double* _alphaParams) {
  localSlater      = _localSlater;
  localDiffSlater  = _localDiffSlater;
  localDDiffSlater = _localDDiffSlater;
  alphaParams      = _alphaParams;
}


// ********************* attachCorrelation ************************
void LocalWaveFunction::attachCorrelation(double* _localCorrelation, 
					  double* _localDiffCorrelation,
					  double* _localDDiffCorrelation, 
					  double* _betaParams) {
  localCorrelation      = _localCorrelation;
  localDiffCorrelation  = _localDiffCorrelation;
  localDDiffCorrelation = _localDDiffCorrelation;
  betaParams            = _betaParams;
}


// ************************ calculateNorm *************************
void LocalWaveFunction::calculateNorm() {
  if (setWeightToUnity) norm=1;
  else {
    norm =  (*localSlater)*(*localCorrelation) 
      /(*centralSlater)/(*centralCorrelation);
    norm *= norm;
  }
}


// ******************** calculateKineticEnergy ********************
void LocalWaveFunction::calculateKineticEnergy() {
  kineticEnergy = -0.5*( (*localDDiffSlater) + (*localDDiffCorrelation) )
    - dotProduct(localDiffSlater, localDiffCorrelation, prodPartDim);
}


// *************************** sample *****************************
void  LocalWaveFunction::sample(double potentialEnergy ) {
  calculateNorm();
  calculateKineticEnergy();
  numSamples                += 1;
  double localEnergy      = kineticEnergy + potentialEnergy;


  // !!! Temporary !!!
  localEnergyReturn       = localEnergy;

  double localEnergyNorm  = localEnergy*norm;
  double localEnergy2norm = localEnergyNorm*localEnergy;

  accumulatedEnergy         += localEnergyNorm;
  accumulatedEnergySquared  += localEnergy2norm;
  accumulatedNorm           += norm;
  accumulatedVariance       
    += pow(localEnergy-referenceEnergy-deltaE, 2)*norm;
  for (int i=0; i<numBlocks; i++) {
    tempEnergyBlock[i]        += localEnergy;
    tempNormBlock[i]          += norm;
    if (++currentBlockNum[i]==blockSize[i]) {
      currentBlockNum[i]=0;
      numSamplesBlock[i]+=1;
      double averageNorm   =tempNormBlock[i]/blockSize[i];
      accumulatedNormBlock[i] +=averageNorm;
      double averageEnergy = tempEnergyBlock[i]/blockSize[i];
      accumulatedEnergyBlock[i] +=averageEnergy*averageNorm;
      accumulatedEnergySquaredBlock[i]
	+=averageEnergy*averageEnergy*averageNorm;
      tempEnergyBlock[i]=tempNormBlock[i]=0;     
    }
  }

}


// ********************** getAverageEnergy ************************
double LocalWaveFunction::getAverageEnergy() {
  return accumulatedEnergy/accumulatedNorm;
}


// ******************** getStandardDeviation **********************
double LocalWaveFunction::getStandardDeviation() {
  return sqrt ( (accumulatedEnergySquared/accumulatedNorm 
		- pow( accumulatedEnergy/accumulatedNorm, 2) ) / numSamples);
}


// ************************* getVariance **************************
double LocalWaveFunction::getVariance() {
  return sqrt ( (accumulatedVariance/accumulatedNorm)
		/ numSamples);
}


// ******************** getAverageEnergyBlock *********************
double LocalWaveFunction::getAverageEnergyBlock(int block) {
  return accumulatedEnergyBlock[block]/accumulatedNormBlock[block];
}


// ****************** getStandardDeviationBlock *******************
double LocalWaveFunction::getStandardDeviationBlock(int block) {
  return sqrt ( (accumulatedEnergySquaredBlock[block]
		 /accumulatedNormBlock[block] 
		- pow( accumulatedEnergyBlock[block]
		       /accumulatedNormBlock[block], 2) ) 
		/ numSamplesBlock[block]);
}


// ************************** addBlock ****************************
void LocalWaveFunction::addBlock(int _blockSize) {
  if (numBlocks<maxNumBlocks) {
    blockSize[numBlocks]                     = _blockSize;
    currentBlockNum[numBlocks]               = 0;
    accumulatedEnergyBlock[numBlocks]        = 0;
    accumulatedEnergySquaredBlock[numBlocks] = 0;
    accumulatedNormBlock[numBlocks]          = 0;
    tempEnergyBlock[numBlocks]               = 0;
    tempNormBlock[numBlocks]                 = 0;
    numSamplesBlock[numBlocks]               = 0;
    numBlocks++;
  }
  else
    *outputFile << "Could not create an additional block in LocalWaveFunction.\n" 
		<< "Maximun number of blocks set to " << maxNumBlocks 
		<< " (in code).\n";
}


// ************************ deleteBlocks **************************
void LocalWaveFunction::deleteBlocks() {
  numBlocks=0;
}


// *************************** summary ****************************
void LocalWaveFunction::summary() {
  *outputFile << "*** Alpha = [ ";
  for (int i=0; i<numAlpha; i++)
    *outputFile << setw(5) << alphaParams[i] << " ";
  *outputFile << "]  Beta = [ ";
  for (int i=0; i<numBeta; i++)
    *outputFile << setw(5) << betaParams[i] << " ";
  *outputFile << "]" <<  " Energy = " << setw(12) << getAverageEnergy()
	      << " ±" << setw(12) << getStandardDeviation();
  if (varianceOptimization) *outputFile << " Variance = " << setw(12) 
					<< getVariance();
  *outputFile << " ***" <<endl;
  *outputFile << "Blocksize ";
  for (int i=0; i<numBlocks; i++) 
    *outputFile << setw(9) << blockSize[i] << " ";
  *outputFile << "\n";
  *outputFile << "S.d.        [ ";
  for (int i=0; i<numBlocks; i++)
    *outputFile << setw(9) << getStandardDeviationBlock(i) << " ";
  *outputFile << "]\n";
}
