#ifndef LocalWaveFunction_IS_INCLUDED
#define LocalWaveFunction_IS_INCLUDED

#include "../Domain/Domain.h"
#include <iostream>
#include <iomanip>

// ****************************************************************
// *                      LOCALWAVEFUNCTION                       *
// ****************************************************************
class LocalWaveFunction {

 protected:

  ofstream *outputFile;

  int     numParticles;     // Number of particles.
  int     numDimensions;    // Number of dimension.
  int     prodPartDim;      // numParticles * numDimensions.

  double* centralSlater;    // Pointer to the value of the central 
                            // Slater-determinant.
  double* localSlater;      // Pointer to the value of the local
                            // Slater-determinant.
  double* localDiffSlater;  // Pointer to the gradient of the Slater-
                            // determinent to the Slater-deterimant.
  double* localDDiffSlater; // Pointer to the Laplacian of the Slater-
                            // determinent to the Slater-deterimant.
  int     numAlpha;         // Number of Slater-determinant parameters.
  double* alphaParams;      // Pointer to the local Slater-determinant 
                            // parameters.

  double* centralCorrelation;    // Pointer to the value of the  
                                 // central correlation G.
  double* localCorrelation;      // Pointer to the value of the local
                                 // correlation G.
  double* localDiffCorrelation;  // Pointer to the gradient of the 
                                 // correlation to the correlation.
  double* localDDiffCorrelation; // Pointer to the Laplacian of the 
                                 // correlation to the correlation.
  int     numBeta;               // Number of correlation parameters.
  double* betaParams;            // Pointer to the local  correlation
                                 // parameters.

  double  norm;                  // localSlater*localCorrelation/ 
                                 // (centralSlater*centralCorrelation).
  double  kineticEnergy;         // The kinetic energy.

  int     numSamples;            // Keep track of the number of samples.
  double  accumulatedEnergy;     // Sum of local energies E_L
  double  accumulatedEnergySquared; // Sum of E_L*E_L
  double  accumulatedNorm;       // Sum of the norms
  double  accumulatedVariance;   // Sum of (E_L - E_ref - deltaE)^2

  // Blocks are included to give an estimate of the true standard 
  // deviation of the energy. This approch is due to auto-correltion 
  // effects. The results should be tested by for example producing 
  // several independent VMC runs, and finding their standard 
  // deviation.
  int     numBlocks;             // The actual number of blocks.
  int     maxNumBlocks;          // Number of blocks allocated.
  int*    blockSize;             // Number of samples per block.
  int*    currentBlockNum;       // Index used for performing block
                                 // sampling.
  int*    numSamplesBlock;       // The number of samples for in each
                                 // of the different blocks.
  double* accumulatedEnergyBlock;// The sum of the energy samples for
                                 // the individual blocks.
  double* accumulatedEnergySquaredBlock; // The sum of the square of the
                                         // block-energies.
  double* accumulatedNormBlock;  // Sum of the norms for each block.
  double* tempEnergyBlock;       // Accumulates energies for a sample.
  double* tempNormBlock;         // Accumulates norms to a sample.

  int     varianceOptimization;  // Boolean determining whether we 
                                 // optimize with respect to energy(0) 
                                 // or variance(1).
  int     setWeightToUnity;      // Boolean true(1) weights (the norms)
                                 // are set to unity.
  double  referenceEnergy;       // Reference energy, E_ref (used for
                                 // variance optimization).
  double  deltaE;                // Small energry displacement (used
                                 // for variance optimization).

  void    calculateNorm();
  void    calculateKineticEnergy();
  void    allocateBlocks();
  double  getAverageEnergyBlock(int block);
  double  getStandardDeviationBlock(int block);
  void    addBlock(int _blockSize);
  void    deleteBlocks();

 public:
  LocalWaveFunction() {}
  void    init();
  void    initialize(Domain* domain);
  void    attachCentralSlaterDet(double* _centralSlaterDet);
  void    attachCentralCorrelation(double* _centralCorrelation);
  void    attachSlater(double* _localSlater, 
		     double* _localDiffSlater,
		     double* _localDDiffSlater, 
		       double* _alphaParams);
  void    attachCorrelation(double* _localCorrelation, 
			    double* _localDiffCorrelation,
			    double* _localDDiffCorrelation, 
			    double* _betaParams);
  void    sample(double potentialEnergy);
  void    summary();

  double  getAverageEnergy();
  double  getStandardDeviation();
  void    setReferenceEnergy(double E) {referenceEnergy=E;}
  double  getVariance();

  double* getAlphaParams() {return alphaParams;}
  double* getBetaParams()  {return betaParams;}

  // !!! Temporary !!!
  double localEnergyReturn;
  double getLocalEnergy()  {return localEnergyReturn;}
};

#endif
