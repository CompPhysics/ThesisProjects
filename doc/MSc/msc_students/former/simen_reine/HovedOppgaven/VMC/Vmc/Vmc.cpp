#ifndef VmcCPP_IS_INCLUDED
#define VmcCPP_IS_INCLUDED

#include "Vmc.h"



// ****************************************************************
// *                             VMC                              *
// ****************************************************************
//
// ***************************** Vmc ******************************
template <class SlaterDeterminant>
Vmc<SlaterDeterminant>::Vmc(Domain& _domain) : domain(_domain) { 

  wf                = new Variations<SlaterDeterminant>( domain() );
  slaterDeterminant = &(wf().getSlaterDet());
  correlation       = &(wf().getCorrelation());
  randomMetro       = domain().getRandomMetro();
  walker            = new Walker<SlaterDeterminant>( domain(),  
						     slaterDeterminant(), 
						     randomMetro(), 
						     correlation(), 
						     wf() );

  outputFile           = domain().getOutputFile();
  numThermalization    = domain().getNumThermalization();
  numCycles            = domain().getNumCycles();
  rank                 = domain().getRank();
  centerRank           = (domain().getSize() - 1.0) / 2.0;

  thermalizationType            = domain().getThermalizationType();
  vmcType                       = domain().getVmcType();
  uniDirectionalMovement        = domain().getUniDirectionalMovement();
  numberVmcRuns                 = domain().getNumberVmcRuns();
  numberOfUniDirectionalMoves   = domain().getNumberOfUniDirectionalMoves();
  *outputFile << endl;
}


// *********************** thermalization *************************
template <class SlaterDeterminant>
void Vmc<SlaterDeterminant>::thermalization() { 
  *outputFile << " Thermalization: " << numThermalization;
  wf().initializeThermalization();
  for (int j = 0; j < numThermalization; j++) {
    walker().doThermalizationStep();
  }
}


// ***************** adaptiveStepThermalization *******************
template <class SlaterDeterminant>
void Vmc<SlaterDeterminant>::adaptiveStepThermalization() { 
  *outputFile << " Thermalization: " << numThermalization;
  wf().initializeThermalization();
  int frac = domain().getNumberSeekAcceptance();
  double soughtAcceptance = domain().getSoughtAcceptance();
  if (domain().getVarySeekWithRank()) 
    soughtAcceptance -= 
      (centerRank - rank) * domain().getVarySeekWithRankStep();
  int numThermFrac = numThermalization/frac;
  if (numThermFrac>0) {
    for (int i=0; i<frac; i++) {
      int acceptance = 0;
      for (int j = 0; j < numThermFrac; j++) {
	acceptance += walker().doThermalizationStep();
      }
      domain().changeStepLen(numThermFrac, acceptance, soughtAcceptance);
    }
  }
  *outputFile << " StepLenght changed to: " <<  domain().getStepLen();
}


// ************************** vmcSome *****************************
template <class SlaterDeterminant>
void Vmc<SlaterDeterminant>::vmcSome() { 
  int numParticles = domain().getNumParticles();
  domain().initVMC();
  wf().initializeThermalization();
  wf().initializeVMC();
  walker().resetAcceptance();
  *outputFile << " VMC run: " << numCycles;
  int acceptance = 0;
  for (int i = 0; i < numCycles; i++) {
    for (int j = 0; j < numParticles; j++) {
      acceptance += walker().doRandomStep();
    }
    wf().sample();
  }
  *outputFile << " Acceptance: " 
	      << (((double)acceptance) / (numCycles) / numParticles)
	      << endl << endl;
  wf().summary();
  *outputFile << endl;
  wf().summaryLowest();
}


// ******************* vmcOneParticleAtATime **********************
template <class SlaterDeterminant>
void Vmc<SlaterDeterminant>::vmcOneParticleAtATime() { 
  domain().initVMC();
  wf().initializeThermalization();
  wf().initializeVMC();
  walker().resetAcceptance();
  *outputFile << " VMC run: " << numCycles;
  int acceptance = 0;
  for (int i = 0; i < numCycles; i++) {
    acceptance += walker().doRandomStep();
    wf().sample();
  }

  *outputFile << " Acceptance: " << 
    (((double)acceptance) / (numCycles))   << endl << endl;
  wf().summary();
  *outputFile << endl;
  wf().summaryLowest();
}


// ************************ initNextRun ***************************
template <class SlaterDeterminant>
void Vmc<SlaterDeterminant>::initNextRun() { 
  domain().setReferenceEnergy(wf().getLowestEnergy());

  domain().setCentralAlphaParam(wf().getLowestAlphaParams());
  domain().calculateAlphaParamArray();

  domain().setCentralBetaParam(wf().getLowestBetaParams());
  domain().calculateBetaParamArray();

  domain().init();
  wf().initNewVmcRun();
}


// ***************************** run ******************************
template <class SlaterDeterminant>
void Vmc<SlaterDeterminant>::run() { 

  domain().initSummary();

  for (int i=0; i<numberOfUniDirectionalMoves; i++) {

    int numRuns=0;
    int uniDirection=1;
    *outputFile << "::::::::::::::::::::::::::::::::" 
		<< "::::::::::::::::::::::::::::::::" << endl;

    while ( (uniDirection) & (numRuns<numberVmcRuns) ) {
      domain().Summary();
      *outputFile << endl;

      // Thermalization
      if (thermalizationType == "thermalization" ) thermalization();
      else if (thermalizationType == "adaptiveStepThermalization" ) 
	adaptiveStepThermalization();
      else *outputFile 
	<< "Error in Vmc::run. Could not perform thermalization." 
	<< endl << "  Input thermalizationType not known." << endl;

      // VMC
      if (vmcType == "vmcSome" ) vmcSome();
      else if (vmcType == "vmcOneParticleAtATime" ) vmcOneParticleAtATime();
      else *outputFile << "Error in Vmc::run. Could not perform vmc." 
		       << endl << "  Input vmcType not known." << endl;

      *outputFile << "One VMC run finished." << endl;

      if (uniDirectionalMovement)
	uniDirection = domain().isMovementUniDirectional
	  (wf().getLowestAlphaParams(), wf().getLowestBetaParams() );
      if (uniDirection) initNextRun();
      numRuns++;
    }  

    if ( (uniDirectionalMovement) & (numRuns==numberVmcRuns) ) 
      *outputFile << "Error in Vmc:run! Could not find minima in numRuns = " 
		  << numberVmcRuns << " runs." << endl;

    if (uniDirectionalMovement) domain().initUniDirectionalIndicator();

    *outputFile << endl << endl;

    domain().increaseNumCycles();
    numCycles = domain().getNumCycles();
    domain().increaseNumThermalization();
    numThermalization = domain().getNumThermalization();
    domain().reduceLocalArea();
    initNextRun();

  }

}

#endif
