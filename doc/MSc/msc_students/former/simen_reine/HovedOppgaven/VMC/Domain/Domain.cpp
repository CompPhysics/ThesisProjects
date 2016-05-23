#include "Domain.h"

  // Constructor: Input file name, prosess number, and number of MPI
  // processes
Domain::Domain(char* initFileName, int _rank, int _size) {
  ifstream initIfile(initFileName); // Opening initIFile
  rank = _rank;
  size = _size;

  // ******************* Input and Output Files ********************
  standardInput         = new char[100];
  electronConfiguration = new char[100];
  randomConfig          = new char[100];
  uniDirectionalConfig  = new char[100];
  slaterParam           = new char[100];
  fixedParamsUp         = new char[100];
  fixedParamsDown       = new char[100];
  correlParam           = new char[100];
  output                = new char[100];
  // Reading from initIFile
  if (!(initIfile >> standardInput)) 
    cerr << "Error reading file: " << initFileName << endl;
  if (!(initIfile >> electronConfiguration))
    cerr << "Error reading file: " << initFileName << endl;
  if (!(initIfile >> randomConfig))
    cerr << "Error reading file: " << initFileName << endl;
  if (!(initIfile >> uniDirectionalConfig))
    cerr << "Error reading file: " << initFileName << endl;
  if (!(initIfile >> slaterParam))
    cerr << "Error reading file: " << initFileName << endl;
  if (!(initIfile >> fixedParamsUp))
    cerr << "Error reading file: " << fixedParamsUp << endl;
  if (!(initIfile >> fixedParamsDown))
    cerr << "Error reading file: " << fixedParamsDown << endl;
  if (!(initIfile >> correlParam))
    cerr << "Error reading file: " << initFileName << endl;
  if (!(initIfile >> output))
    cerr << "Error reading file: " << initFileName << endl;
  // Opening the files
  ifstream standardInputIfile(standardInput);
  ifstream electronConfigurationIfile(electronConfiguration);
  ifstream randomConfigIfile(randomConfig);
  ifstream uniDirectionalConfigIfile(uniDirectionalConfig);
  sprintf(output, "%s.run%d", output, rank);
  outputFile.open(output);


  // *********************** Standard Input ************************
  if (!(standardInputIfile >> numDimensions))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: numDimensions" << endl;
  if (!(standardInputIfile >> numThermalization))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: numThermalization" << endl;
  if (!(standardInputIfile >> numCycles))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: numCycles" << endl;
  if (!(standardInputIfile >> stepLen))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: stepLen" << endl;
  if (!(standardInputIfile >> allowSpinFlip))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: allowSpinFlip" << endl;
  if (!(standardInputIfile >> quarterCusp))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: quarterCusp" << endl;
  thermalizationType = new char[100];
  if (!(standardInputIfile >> thermalizationType))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: thermalizationType" << endl;
  if (!(standardInputIfile >> soughtAcceptance))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: soughtAcceptance" << endl;
  if (!(standardInputIfile >> numberSeekAcceptance))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: numberSeekAcceptance" << endl;
  if (!(standardInputIfile >> varySeekWithRank))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: varySeekWithRank" << endl;
  if (!(standardInputIfile >> varySeekWithRankStep))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: varySeekWithRankStep" << endl;
  vmcType = new char[100];
  if (!(standardInputIfile >> vmcType))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: vmcType" << endl;
  if (!(standardInputIfile >> numberVmcRuns))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: numberVmcRuns" << endl;
  if (!(standardInputIfile >> varianceOptimization))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: varianceOptimization" << endl;
  if (!(standardInputIfile >> referenceEnergy))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: referenceEnergy" << endl;
  if (!(standardInputIfile >> deltaE))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: deltaE" << endl;
  if (!(standardInputIfile >> setWeightToUnity))
    cerr << "Error reading file: " << standardInput 
	 << ". Couldn't read: setWeightToUnity" << endl;
  twoD                 = numDimensions*2;


  // ******************* Electron Configuration ********************
  if (!(electronConfigurationIfile >> numParticles))
    cerr << "Error reading file: " << electronConfiguration 
	 << ". Couldn't read: numParticles" << endl;
  if (!(electronConfigurationIfile >> numParticlesSpinUp ))
    cerr << "Error reading file: " << electronConfiguration 
	 << ". Couldn't read: numParticlesSpinUp" << endl;
  if (!(electronConfigurationIfile >> up1s >> up2s 
	>> up2px >> up2py >> up2pz >> up3s >> up3px 
	>> up3py >> up3pz >> up4s)) 
    cerr << "Error reading file: " << electronConfiguration 
	 << ". Couldn't read: up1s >> up2s >> up2px >> up2py >> up2pz" << endl;
  if (!(electronConfigurationIfile >> numParticlesSpinDown ))
    cerr << "Error reading file: " << electronConfiguration 
	 << ". Couldn't read: numParticlesSpinDown" << endl;
  if (!(electronConfigurationIfile >> down1s >> down2s 
	>> down2px >> down2py >> down2pz >> down3s 
	>> down3px >> down3py >> down3pz >> down4s))
    cerr << "Error reading file: " << electronConfiguration 
	 << ". Couldn't read: down1s >> down2s >> down2px >> down2py >> down2pz" << endl;


  // *********************** Random Config ************************
  randomGenerator = new char[100];
  if (!(randomConfigIfile >> randomGenerator))
    cerr << "Error reading file: " << randomConfig 
	 << ". Couldn't read: randomGenerator" <<  endl;
  int seed;
  if (!(randomConfigIfile >> seed))
    cerr << "Error reading file: " << randomConfig 
	 << ". Couldn't read: seed" <<  endl;
  if ((!seed < 0) & (seed > -1000000)) 
    cerr << "Error! The seed should be in the range (-1 000 000, -1)" << endl;
  else {
    initRanMove  = seed - rank*3000000;
    initRanFlip  = seed - rank*3000000 - 1000000;
    initRanMetro = seed - rank*3000000 - 2000000;
  }
  if ( randomGenerator=="Ran0") {
    randomMove           = new Ran0(initRanMove);
    randomFlip           = new Ran0(initRanFlip);
    randomMetro          = new Ran0(initRanMetro);
  }
  else if ( randomGenerator=="Ran1") {
    randomMove           = new Ran1(initRanMove);
    randomFlip           = new Ran1(initRanFlip);
    randomMetro          = new Ran1(initRanMetro);
  }
  else cerr << "Error input randomGenerator, only Ran0 and Ran1 implemented" 
	    << endl;


  // **************** Uni-Directional Configuration ***************
  if (!(uniDirectionalConfigIfile >> uniDirectionalMovement))
    cerr << "Error reading file: " << uniDirectionalConfig 
	 << ". Couldn't read: uniDirectionalMovement(1=True,0=False)" << endl;
  if (!(uniDirectionalConfigIfile >> numberOfUniDirectionalMoves))
    cerr << "Error reading file: " << uniDirectionalConfig 
	 << ". Couldn't read: numberOfUniDirectionalMoves" << endl;
  if (!(uniDirectionalConfigIfile >> reduceLocalAreaByFraction))
    cerr << "Error reading file: " << uniDirectionalConfig 
	 << ". Couldn't read: reduceLocalAreaByFraction" << endl;
  if (!(uniDirectionalConfigIfile >> increaseNumCyclesByFactor))
    cerr << "Error reading file: " << uniDirectionalConfig 
	 << ". Couldn't read: increaseNumCyclesByFactor" << endl;
  if (!(uniDirectionalConfigIfile >> increaseNumThermalizationByFactor))
    cerr << "Error reading file: " << uniDirectionalConfig 
	 << ". Couldn't read: increaseNumThermalizationByFactor" << endl;

  // *************************** Other *****************************
  paramIndex = new int[10];


  // *************** Slater-Determinant Parameters ***************
  createSlaterParams();


  // ******************* Correlation Parameters *******************
  createCorrelParams();


  // *************************** Other *****************************
  numVariations = numSlaterVariations * numCorrelVariations;
  h = 0.0001;   // Numerical differentiate parameter


  // ********************* Coordinates and Spin ********************
  coorArray      = new double[numDimensions*numParticles];
  trialCoorArray = new double[numDimensions];
  spinArray      = new int[numParticles + 1];
  coors          = new CoorSpinDiff[numParticles];

  trialCoor      = new CoorSpinDiff[1];
  setCoorsToOrigon();
  for (int i = 0; i <= numParticles; i++) spinArray[i]=1;
  for (int i = 0; i < numParticles; i++)
    coors[i].attach((coorArray + i*numDimensions), numDimensions, 
		    &spinArray[i], h);
  trialCoor[0].attach(trialCoorArray, numDimensions, 
		      &spinArray[numParticles], h);


  // **************** Uni-Directional Configuration ***************
  uniDirectionalIndicator = new int[2*(numAlpha+numBeta)];
  initUniDirectionalIndicator();


  // ************************ Spin-Factors *************************
  createSpinFactors();
  createDistanceAndDistanceDiff();
}



// *************************************************************
// *                       Standard Input                      *
// *************************************************************
//
// ************************ changeStepLen **********************
void Domain::changeStepLen(int numSteps, int acceptedSteps, 
			   double soughtAcceptance) {
  double acceptance = ((double) acceptedSteps)/numSteps;
  if (acceptance == 1) cerr << "Error in changing stepLen, either stepLen=0 or numThermalization small" << endl;
  else stepLen *= (1-soughtAcceptance)/(1-acceptance);
}



// *************************************************************
// *                Uni-Directional Configuration              *
// *************************************************************
//
// ********************** increaseNumCycles ********************
void Domain::increaseNumCycles() {
  numCycles*=increaseNumCyclesByFactor;

}

// *****************  increaseNumThermalization ****************
void Domain::increaseNumThermalization() {
  numThermalization*=increaseNumThermalizationByFactor;

}

// *********************** reduceLocalArea *********************
void Domain::reduceLocalArea() {
  for (int i=0; i<numAlpha; i++)
    alphaStep[i]*=reduceLocalAreaByFraction;
  for (int i=0; i<numBeta; i++)
    betaStep[i]*=reduceLocalAreaByFraction;
}

/// ***************** initUniDirectionalIndicator ***************
void Domain::initUniDirectionalIndicator() {
  int num = numAlpha + numBeta;
  for (int i=0; i<2*num; i++)
    uniDirectionalIndicator[i] = 1;
}

// ******************* isMovementUniDirectional ****************
int Domain::isMovementUniDirectional(double* alphaParams, double* betaParams) {
  // Update the uniDirectionalIndicator
  int num = numAlpha + numBeta;
  double small = 1e-15;
  for (int i=0; i<numAlpha; i++) {
    if ( alphaParams[i] <= centralAlphaParams[i] + small ) 
      uniDirectionalIndicator[i] = 0;
    if ( alphaParams[i] >= centralAlphaParams[i] - small) 
      uniDirectionalIndicator[i+num] = 0;
  }
  for (int i=0; i<numBeta; i++) {
    if ( betaParams[i] <= centralBetaParams[i] + small) 
      uniDirectionalIndicator[i+numAlpha] = 0;
    if ( betaParams[i] >= centralBetaParams[i] - small) 
      uniDirectionalIndicator[i+numAlpha+num] = 0;
  }

  // Here the test is performed
  int uni = 0;
  for (int i=0; i<2*num; i++)
    if ( uniDirectionalIndicator[i] != 0 ) uni = 1;

  return uni;
}



// *************************************************************
// *               Slater-Determinant Parameters               *
// *************************************************************
//
// ********************* createSlaterParams ********************
void Domain::createSlaterParams() {
  // Opening slaterParamIfile
  ifstream slaterParamIfile(slaterParam);
  //    cerr << "Error! Could not open file: " << slaterParam << endl;
  orbitalType = new char[100];
  if (!(slaterParamIfile >> orbitalType))
    cerr << "Error reading file: " << slaterParam 
	 << ". Couldn't read: orbitalType" << endl;
  if (!(slaterParamIfile >> numAlpha))
    cerr << "Error reading file: " << slaterParam 
	 << ". Couldn't read: numAlpha" << endl;
  if (numAlpha<1) cerr << "numAlpha must be at least 1" << endl;
  centralAlphaParams  = new double[numAlpha];
  numAlphaVar         = new int[numAlpha];
  alphaStep           = new double[numAlpha];
  for (int i=0; i<numAlpha; i++ ) {
    if (!(slaterParamIfile >> centralAlphaParams[i] >> 
	  numAlphaVar[i] >> alphaStep[i])) {
      cerr << "Error reading file: " << slaterParam 
	   << ". Couldn't read: centralAlphaParams[i] >>" 
	   << " numAlphaVar[i] >> alphaStep[i] (i=" << i << ")" << endl;
      cerr << "Possible missmach between numAlpha and" 
	   << " number of additional lines?" << endl;
    }
  }
  numSlaterVariations = 1;
  for (int i=0; i<numAlpha; i++)
    numSlaterVariations *= numAlphaVar[i];
  if ( (numSlaterVariations<1) || 
       ( !(odd(numSlaterVariations)) ) )
    cerr << "alpha#numVar must be odd and >0" << cerr << endl;
  alphaParam = new double[numAlpha*numSlaterVariations];
  calculateAlphaParamArray();
  centerSlater = centerOfOddIntegerMinusOne(numSlaterVariations);
}

// ****************** calculateAlphaParamArray *****************
void Domain::calculateAlphaParamArray() {
  for (int i=0; i<numAlpha; i++)
    centralAlphaParams[i] 
      -=  centerOfOddIntegerMinusOne(numAlphaVar[i]) * alphaStep[i];
  for (int i=0;i<numAlpha;i++) paramIndex[i]=0;
  bool increaseParam;
  int index;
  int j = 0;
  for (int i=0; i<numSlaterVariations; i++) {
    for (int k=0; k<numAlpha; k++) {
      alphaParam[j++] = centralAlphaParams[k];
    }
    index = 0;
    increaseParam = false;
    while (!increaseParam) {
      paramIndex[index] += 1;
      if (paramIndex[index] == numAlphaVar[index]) {
	centralAlphaParams[index] -= (numAlphaVar[index]-1)*alphaStep[index];
	paramIndex[index] = 0;
	index +=1;
      }
      else {
	centralAlphaParams[index] += alphaStep[index];
	increaseParam = true;
      }
    }
  }
  for (int i=0; i<numAlpha; i++)
    centralAlphaParams[i] += 
      centerOfOddIntegerMinusOne(numAlphaVar[i])*alphaStep[i];
}

// ******************** setCentralAlphaParam *******************
void Domain::setCentralAlphaParam(double* param) {
  copyArray(param, centralAlphaParams, numAlpha);
}



// *************************************************************
// *                   Correlation Parameters                  *
// *************************************************************
//
// ********************* createCorrelParams ********************
void Domain::createCorrelParams() {
  // Opening correlParamIfile
  ifstream correlParamIfile(correlParam);
  //    cerr << "Error! Could not open file: " << correlParam << endl;
  fBetaType = new char[100];
  if (!(correlParamIfile >> fBetaType))
    cerr << "Error reading file: " << correlParam 
	 << ". Couldn't read: fBetaType" << endl;
  if (!(correlParamIfile >> numBeta))
    cerr << "Error reading file: " << correlParam 
	 << ". Couldn't read: numBeta" << endl;
  centralBetaParams = new double[numBeta];
  numBetaVar    = new int[numBeta];
  betaStep      = new double[numBeta];
  for (int i=0; i<numBeta; i++ ) {
    if (!(correlParamIfile >> centralBetaParams[i] >> 
	  numBetaVar[i] >> betaStep[i])) {
      cerr << "Error reading file: " << correlParam 
	   << ". Couldn't read: centralBetaParams[i] >>" 
	   << " numBetaVar[i] >> betaStep[i] (i=" << i << ")" << endl;
      cerr << "Possible missmach between numBeta and" 
	   << " number of additional lines?" << endl;
    }
  }

  numCorrelVariations = 1;
  for (int i=0; i<numBeta; i++)
    numCorrelVariations *= numBetaVar[i];
  if ( (numCorrelVariations<1) || 
       ( !(odd(numCorrelVariations)) ) )
    cerr << "beta#numVar must be odd and >0" << cerr << endl;

  betaParam = new double[numBeta*numCorrelVariations];
  calculateBetaParamArray();
  centerCorrel = centerOfOddIntegerMinusOne(numCorrelVariations);

  if      ( fBetaType == "fNone" )        createFNone();
  else if ( fBetaType == "fBeta" )        createFBeta(1);
  else if ( fBetaType == "fBetaLinear" )  createFBeta(0);
  else if ( fBetaType == "fBeta2" )       createFBeta2(1);
  else if ( fBetaType == "fBeta3" )       createFBeta3(1);
  else if ( fBetaType == "fBeta2r2" )     createFBeta2r2(1);
  else if ( fBetaType == "fExtended" )    createFExtended(1);
  else if ( fBetaType == "fBeta2Linear" ) createFBeta2(0);
  else if ( fBetaType == "fBetaMany" )    createFBetaMany(1);
  else    cerr << "fBetaType not identified!" << endl;
}

// ******************* calculateBetaParamArray *****************
void Domain::calculateBetaParamArray() {
  for (int i=0; i<numBeta; i++)
    centralBetaParams[i] 
      -= centerOfOddIntegerMinusOne(numBetaVar[i]) * betaStep[i];
  for (int i=0;i<numBeta;i++) paramIndex[i]=0;
  bool increaseParam;
  int index;
  int j = 0;
  for (int i=0; i<numCorrelVariations; i++) {
    for (int k=0; k<numBeta; k++) {
      betaParam[j++] = centralBetaParams[k];
    }
    index = 0;
    increaseParam = false;
    while (!increaseParam) {
      paramIndex[index] += 1;
      if (paramIndex[index] == numBetaVar[index]) {
	centralBetaParams[index] -= (numBetaVar[index]-1)*betaStep[index];
	paramIndex[index] = 0;
	index +=1;
      }
      else {
	centralBetaParams[index] += betaStep[index];
	increaseParam = true;
      }
    }
  }
  for (int i=0; i<numBeta; i++)
    centralBetaParams[i] += 
      centerOfOddIntegerMinusOne(numBetaVar[i]) * betaStep[i];
}

// ******************** setCentralBetaParam ********************
void Domain::setCentralBetaParam(double* param) {
  copyArray(param, centralBetaParams, numBeta);
}

// ************************ attachFParams **********************
void Domain::attachFParams(int _expBool) {
  for (int i=0; i<numCorrelVariations; i++)
    f[i].attach(numBeta, &(betaParam[i*numBeta]),_expBool);
}

// ************************* createFNone ***********************
void Domain::createFNone() {
  if (numBeta != 1) cerr << "Need numBeta=1 to create fNone!" << endl;
  numCorrelVariations = 1;
  f = new fNone[1];
  attachFParams(0);
}

// ************************* createFBeta ***********************
void Domain::createFBeta(int _expBool) {
  if (numBeta != 1) 
    cerr << "Could not create fBeta, numBeta must be 1!" << endl;
  else {
    f = new fBeta[numCorrelVariations];
    attachFParams(_expBool);
  }
}

// *********************** createFBetaMany *********************
void Domain::createFBetaMany(int _expBool) {
  if (numBeta != 9) 
    cerr << "Could not create fBetaMany, numBeta must be 9!" << endl;
  else {
    f = new fBetaMany[numCorrelVariations];
    attachFParams(_expBool);
  }
}

// ************************ createFBeta2 ***********************
void Domain::createFBeta2(int _expBool) {
  if (numBeta != 2) 
    cerr << "Could not create fBeta2, numBeta must be 2!" << endl;
  else {
    f = new fBeta2[numCorrelVariations];
    attachFParams(_expBool);
  }
}

// ************************ createFBeta3 ***********************
void Domain::createFBeta3(int _expBool) {
  if (numBeta != 3) 
    cerr << "Could not create fBeta3, numBeta must be 3!" << endl;
  else {
    f = new fBeta3[numCorrelVariations];
    attachFParams(_expBool);
  }
}

// *********************** createFBeta2r2 **********************
void Domain::createFBeta2r2(int _expBool) {
  if (numBeta != 2) 
    cerr << "Could not create fBeta2r2, numBeta must be 2!" << endl;
  else {
    f = new fBeta2r2[numCorrelVariations];
    attachFParams(_expBool);
  }
}

// ********************** createFExtended **********************
void Domain::createFExtended(int _expBool) {
  if (numBeta != 8) 
    cerr << "Could not create fExtended, numBeta must be 8!" << endl;
  else {
    f = new fExtended[numCorrelVariations];
    attachFParams(_expBool);
  }
}



// *************************************************************
// *                 Inter-electronic distances                *
// *************************************************************
//
// **************** createDistanceAndDistanceDiff **************
void Domain::createDistanceAndDistanceDiff() {
  distance        = new Distance;
  distance->attach(coors, trialCoor, 
		   numParticles);
  distanceDiff    = new DistanceDiff[numDimensions];
  for (int i=0; i<numDimensions; i++)
      distanceDiff[i]
	.attach(coors, trialCoor, numParticles, h,i); 
  
}



// *************************************************************
// *                    Coordinates and Spin                   *
// *************************************************************
//
// ********************* setCurrentParticle ********************
void Domain::setCurrentParticle(int __currentParticle) {
  currentParticle = __currentParticle;
  _coorArray      = (coorArray + currentParticle*numDimensions);
  spinFactors->setCurrentParticle(currentParticle);
  distance->setCurrentParticle(currentParticle);
  for (int i=0; i< numDimensions; i++)
    distanceDiff[i].setCurrentParticle(currentParticle);
}

// ********************** setOtherParticle *********************
void Domain::setOtherParticle(int __otherParticle) {
  otherParticle  = __otherParticle;
  spinFactors->setOtherParticle(otherParticle);
}

// ********************** setToNextParticle ********************
void Domain::setToNextParticle() {
  if (++currentParticle == numParticles) {
    currentParticle = 0;
    resetPtr();
  }
  else
    _coorArray += numDimensions;
  spinFactors->setToNextParticle();
  distance->setToNextParticle();
  for (int i=0; i< numDimensions; i++)
    distanceDiff[i].setToNextParticle();
}

// **************************** init ***************************
void Domain::init() {
  setCurrentParticle(0);
  for (int i = 0; i <= numParticles; i++)
    spinArray[i]=1;

  for (int i = 0; i < numParticles; i++) {
    setPositionToOrigin();
    computeTrialPosition();
    acceptTrialPosition();
    if (currentParticle >= numParticlesSpinUp)
      coors[currentParticle].flipSpin();
    coors[currentParticle].calculateR();
    setToNextParticle();
  }
  spinFactors->init();
  spinFactors->setOtherParticle(0);
  distance->initialize();
}

// ************************** initVMC **************************
void Domain::initVMC() {
  setCurrentParticle(0);
  distance->initialize();
  for (int i = 0; i < numDimensions; i++)
    distanceDiff[i].initialize();
}

// ************************ suggestMove ************************
void Domain::suggestMove() {
  computeTrialPosition();
  distance->suggestMove();
  proposeFlip();
}

// ************************ acceptMove *************************
void Domain::acceptMove() {
  for (int i = 0; i < numDimensions; i++)
    distanceDiff[i].suggestMove();
  acceptThermalizedMove();
  for (int i = 0; i < numDimensions; i++)
    distanceDiff[i].acceptMove();
}

// ************************ rejectMove *************************
void Domain::rejectMove() {
  rejectThermalizedMove();
 for (int i = 0; i < numDimensions; i++)
    distanceDiff[i].rejectMove();
}

// ******************* acceptThermalizedMove *******************
void Domain::acceptThermalizedMove() {
  acceptTrialPosition();
  distance->acceptMove();
}


// ******************* rejectThermalizedMove *******************
void Domain::rejectThermalizedMove() {
  rejectTrialPosition();
  distance->rejectMove();
}

// ************************ proposeFlip ************************
void Domain::proposeFlip() {
  if (allowSpinFlip) {
    setOtherParticle ( (int) (randomFlip().getNum()*numParticles) );
    spinFlip = (*(spinArray+currentParticle) != *(spinArray+otherParticle));
  }
  if (quarterCusp) {
    if (spinFlip) spinFactors->calculateFlipFactors();
    else spinFactors->calculateNoFlipFactors();
  }
}

// ************************** resetPtr *************************
void Domain::resetPtr() {
  _coorArray = coorArray;
}

// ******************** computeTrialPosition *******************
void Domain::computeTrialPosition() {
  for (int i = 0; i < numDimensions; i++) {
    trialCoorArray[i] = _coorArray[i] + coorStep();
  }
  trialCoor[0].spin() = coors[currentParticle].spin();
  trialCoor[0].calculateR();
  trialCoor[0].calculateDiffs();
}

// ******************** acceptTrialPosition ********************
void Domain::acceptTrialPosition() {
  for (int i = 0; i < numDimensions; i++)
    _coorArray[i] = trialCoorArray[i];
  coors[currentParticle].copy(trialCoor);
}

// ******************** rejectTrialPosition ********************
void Domain::rejectTrialPosition() {
}

// ******************** setPositionToOrigin ********************
void Domain::setPositionToOrigin() {
   for (int i = 0; i < numDimensions; i++)
    _coorArray[i] = 0;
   coors[currentParticle].calculateR();
}

// ************************** getCoors *************************
CoorSpinDiff* Domain::getCoors() {
  return coors;
}

// ************************ getTrialCoors **********************
CoorSpinDiff* Domain::getTrialCoor() {
  return trialCoor;
}

// *************************** coorStep ************************
double Domain::coorStep() {
  return (stepLen * (randomMove().getNum() - 0.5)/ numParticles);
}


// ********************** setCoorsToOrigon *********************
void Domain::setCoorsToOrigon() {
  for (int i=0; i<numDimensions*numParticles; i++)
    coorArray[i]=0;

  for (int i=0; i<numDimensions; i++)
    trialCoorArray[i]=0;
}


// *************************************************************
// *                      Spin-Factors                         *
// *************************************************************
//
// ******************** createSpinFactors **********************
void Domain::createSpinFactors() {
  spinFactors = new SpinFactors;
  spinFlip    = 0;
  spinFactors->allocate(numParticles, spinArray);
  spinFactors->setQuarterCusp(quarterCusp);
  spinFactors->init();
}


// *************************************************************
// *                     Potential Energy                      *
// *************************************************************
//
// ***************** getNucleusElectronPotential ***************
double Domain::getNucleusElectronPotential() {
  double nucleusElectronPotential = 0;
  for (int i=0; i<numParticles; i++)
    nucleusElectronPotential -= 1/coors[i].r();
  return ( (double) numParticles )*nucleusElectronPotential;
}


// *************************************************************
// *                          Summary                          *
// *************************************************************
//
// ************************* initSummary ***********************
void Domain::initSummary() {
  outputFile << "!!!!!!!!!!!!   INITIAL CONFIGURATION   !!!!!!!!!!!!!!!!!!" 
	     << endl
	     << "numDimensions          " << numDimensions << endl
	     << "numParticles           " << numParticles << endl
	     << "numParticlesSpinUp     " << numParticlesSpinUp << endl
	     << "spinUpConfig           [ "
	     << up1s << " " << up2s << " " << up2px << " " << up2py 
	     << " " << up2pz << " ]" << endl
	     << "numParticlesSpinDown   " << numParticlesSpinDown << endl
	     << "spinDownConfig         [ "
	     << down1s << " " << down2s << " " << down2px << " " << down2py 
	     << " " << down2pz << " ]" << endl
	     << "randomGenerator        " << randomGenerator << endl
	     << "initRanMove/Flip/Metro [ " << initRanMove << " " 
	     << initRanFlip << " " << initRanMetro << " ]\n"
	     << "orbitalType            " << orbitalType << endl
	     << "fBetaType              " << fBetaType << endl;
  if (quarterCusp) 
    outputFile << "Quarter-cusp for like-spin particles turned on!\n";
  else 
    outputFile << "Quarter-cusp for like-spin particles turned off!\n";
  if (allowSpinFlip) 
    outputFile << "Spin-flip turned on!\n";
  else 
    outputFile << "Spin-flip turned off!\n";
  outputFile << "Thermalization set to type '" << thermalizationType 
	     << "', and VMC type set\n" << "to '" << vmcType << "'.\n";
  if (uniDirectionalMovement) 
    outputFile << "The minima in parameter space "
	       << "will be sought (uniDirectionalMovement=1).\n"
	       << "There will be " << numberOfUniDirectionalMoves 
	       << " unidirectional moves. " 
	       << "One unidirectional move will consist\n" 
	       << "of maximum " << numberVmcRuns << " VMC runs.\n"
	       << "In-between each unidirectional move" 
	       << " the number of cycles will be increased by \n" 
	       << "a factor " << increaseNumCyclesByFactor 
	       << ", and the number" 
	       << " of thermalization by a factor " 
	       << increaseNumThermalizationByFactor << ". The "
	       << "variational \nparameters step length will be reduced" 
	       << " by fraction " << reduceLocalAreaByFraction << ".\n";
  else outputFile << numberVmcRuns << " VMC runs (with differing" 
		  << " seeds) will be conducted.\n";
  outputFile << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" 
	     << endl;
}

// **************************** Summary ****************************
void Domain::Summary() {
  outputFile << "*********************************************************" 
	     << endl
	     << "Summary of domain:" << endl
	     << "     numThermalization  " << numThermalization << endl
	     << "     numCycles          " << numCycles << endl
	     << "     stepLen            " << stepLen << endl
	     << "     centralAlphaParams [ ";

  for (int i=0; i<numAlpha; i++) 
    outputFile << centralAlphaParams[i]  << " ";

  outputFile << "]" << endl
	     << "     numAlphaVar        [ ";

  for (int i=0; i<numAlpha; i++) outputFile << numAlphaVar[i]  << " ";
  outputFile << "]" << endl
	     << "     alphaStep          [ ";
  for (int i=0; i<numAlpha; i++) outputFile << alphaStep[i]  << " ";
  outputFile << "]" << endl
	     << "     centralBetaParams  [ ";
  for (int i=0; i<numBeta; i++) outputFile << centralBetaParams[i]  << " ";
  outputFile << "]" << endl
	     << "     numBetaVar         [ ";
  for (int i=0; i<numBeta; i++) outputFile << numBetaVar[i]  << " ";
  outputFile << "]" << endl
	     << "     betaStep           [ ";
  for (int i=0; i<numBeta; i++) outputFile << betaStep[i]  << " ";
  outputFile << "]" << endl
	     << "*********************************************************" 
	     << endl;
}
