#ifndef Domain_IS_INCLUDED
#define Domain_IS_INCLUDED

#include "../Paramizer/Paramizer.h"
#include "../Coor/Coor.h"
#include "../Random/Random.h"
#include "../Distance/Distance.h"
#include "../SpinFactors/SpinFactors.h"
#include "../Ref/Ref.h"
#include "../Random/Random.h"
#include "../fFunction/fFunction.h"

#include <fstream>
#include <string>
#include <iostream>
#include <cstdio>
using namespace std;

#ifndef copyArray_IS_INCLUDED
#define copyArray_IS_INCLUDED
// *************************** copyArray ***********************
inline void copyArray(double* fromArray, double* toArray, int dim) {
  double* _fromArray = fromArray - 1;
  double* _toArray = toArray - 1;
  for (int i=0; i<dim; i++)
    (*++_toArray) = (*++_fromArray);
}
#endif

#ifndef dotProduct_IS_INCLUDED
#define dotProduct_IS_INCLUDED
// ************************** dotProduct ***********************
inline double dotProduct(double* array1, double* array2, int dim) {
  double* a1 = array1 - 1;
  double* a2 = array2 - 1;
  double product = 0;
  for (int i=0; i<dim; i++)
    product += (*++a1)*(*++a2);
  return product;
}
#endif

#ifndef odd_IS_INCLUDED
#define odd_IS_INCLUDED
// ***************************** odd ***************************
inline int odd(int number) {
  return (number - number/2);
}
// ********************** centerOfOddInteger *******************
inline int centerOfOddInteger(int oddInteger) {
  if (!odd(oddInteger)) cerr << "Error finding center of integer "
			     << "in Domain; integer is even (=" 
			     << oddInteger << ")!" << endl;
  else return ((oddInteger+1)/2);
}
// ****************** centerOfOddIntegerMinusOne ***************
inline int centerOfOddIntegerMinusOne(int oddInteger) {
  return (centerOfOddInteger(oddInteger)-1);
}
#endif





// *************************************************************
// *                                                           *
// *                           DOMAIN                          *
// *                                                           *
// *************************************************************
class Domain {
 // ************************************************************
 // *                 PROTECTED ALGORITHMS                     *
 // ************************************************************
 protected:
  // ************************** MPI ****************************
  int rank; // Prosess number
  int size; // Number of prosesses


  // ****************** Input and Output Files *****************
  char*    standardInput;
  char*    electronConfiguration;
  char*    randomConfig;
  char*    uniDirectionalConfig;
  char*    slaterParam;
  char*    fixedParamsUp;
  char*    fixedParamsDown;
  char*    correlParam;
  char*    output;
  ofstream outputFile;


  // ********************* Standard Input **********************
  int    numDimensions;
  int    numThermalization;
  int    numCycles;
  double stepLen;              // The step length of the Metropolis algorithm.
  int    allowSpinFlip;        // Boolean turns on(1) spin-flipping routines.
  int    quarterCusp;          // Boolean turns on(1) cusp=1/4 off (0).
  string thermalizationType;   // Normal or adaptive.
  double soughtAcceptance;     // What acceptance ratio we wish to have.
  int    numberSeekAcceptance; // Number of updates of the step length
                               // to get sought acceptance.
  int    varySeekWithRank;     // MPI option to vary the acceptance 
                               // ratio with the rank (process number).
  double varySeekWithRankStep; // How much we vary the acceptance.
  string vmcType;              // One at a time or Some.
  int    numberVmcRuns;        // Number of VMC runs.
  int    varianceOptimization; // Boolean 1=True, 0=False.
  double referenceEnergy;      // Reference energy for the variance-
                               // optimization. Updated between different VMC 
                               // moves.
  double deltaE;               // Quantity to allow both energy and variance 
                               // optimization.
  int    setWeightToUnity;     // Boolean 1=True, 0=False

  // Accosiated parameters
  int    twoD;


  // ****************** Electron Configuration ******************
  int numParticles;
  int numParticlesSpinUp;
  int numParticlesSpinDown;
  // Booleans indicating which electrons are allowed
  int up1s, up2s, up2px, up2py, up2pz;
  int up3s, up3px, up3py, up3pz, up4s;
  int down1s, down2s, down2px, down2py, down2pz;
  int down3s, down3px, down3py, down3pz, down4s;


  // ********************** Random Config ***********************
  string randomGenerator; // Either Ran0 or Ran1
  // int seed; The seed is given by the input file

  // The different rendom generator seeds
  int             initRanFlip;  // = seed - rank*3mill
  int             initRanMove;  // = seed - rank*3mill - 1 mill
  int             initRanMetro; // = seed - rank*3mill - 2 mill
  // References to the random generator
  Ref<Random2>    randomFlip;
  Ref<Random2>    randomMove;
  Ref<Random2>    randomMetro;


  // *************** Uni-Directional Configuration **************
  int    uniDirectionalMovement;       // Do we seek for minima in parameter
                                       // space (1=True, 0=False)?
  int    continueAfterUniDirectional;  // Whem minima is found, do we 
                                       // wish to make additional
                                       // sophisticated searches
                                       // (1=True, 0=False)?
  int    numberOfUniDirectionalMoves;  // Total number of sophisticated
                                       // minima searches.
  int    reduceLocalAreaBetweenMoves;  // Reduce the area in parameter-
                                       // space (1=True, 0=False)?
  double reduceLocalAreaByFraction;    // Reduce the area in parameter-
                                       // space by fraction.
  int    increaseNumCyclesBetweenMoves;// Increase the number of
                                       // cycles (1=True, 0=False)?
  int    increaseNumCyclesByFactor;    // Increase by factor.
  int    increaseNumThermalizationBetweenMoves; // Increase thermalization 
                                                // (1=True, 0=False)?
  int    increaseNumThermalizationByFactor;     // Increase by factor.

  // Indicator specifying whether or not minima is found
  int*   uniDirectionalIndicator;


  // *************************** Other *****************************
  int    numVariations;         // Total number of local variations 
                                // numSlaterVariations*numCorrelVariations
  double h;                     // Numberical difference (for calculating 
                                // derivatives)
  int*   paramIndex;            // Index used for initialization and updating 
                                // of alphaParam and betaParam in the 
                                // neighborhood of some central values for 
                                // the parameters.


  // *************** Slater-Determinant Parameters ***************
  string  orbitalType; // For example hydrogen, HF, etc.
  int     numAlpha;    // Number of variational (Slater) parameters (min 1).
  double* centralAlphaParams; // The central parameters (with no local 
                              // variation this is the parameter we use).
  int*    numAlphaVar; // Number of local variations for each iniviual
                       // variational parameter.
  double* alphaStep;   // The step (parameter difference) between  
                       // each local variation.
  double* alphaParam;  // Array containing all the different local 
                       // configurations in (Slater) parameter space.

  int     numSlaterVariations;  // Total number of local Slater variations.
  int     centerSlater;         // Index indicating which of the Slater 
                                // variations is the central one.
  void    createSlaterParams(); // Routine to initialize the Slater 
                                // determinant parameters.


  // ******************* Correlation Parameters *******************
  string  fBetaType;         // For example fNone, fBeta, fBetaLinear, etc.
  int     numBeta;           // Number of variational (correlation) 
                             // parameters (min 1).
  double* centralBetaParams; // The central parameters (with no local 
                             // variation this is the parameter we use).
  int*    numBetaVar;        // Number of local variations for each iniviual
                             // variational parameter.
  double* betaStep;          // The step (parameter difference) between  
                             // each local variation.
  
  double*    betaParam;            // Array containing all the different local 
                                   // configurations in (correlation) parameter 
                                   // space.
  int        numCorrelVariations;  // Total number of local correlation 
                                   // variations.
  int        centerCorrel;         // Index indicating which of the correlation
                                   // variations is the central one.
  fFunction* f;                    // Function set in accordance to the input
                                   // of fBetaType.
  void       createCorrelParams(); // Routine to initialize the correlation
                                   // parameters.
  void attachFParams(int _expBool);// Attach paramaters to given fFunction
  // Routines to create different fFunctions
  void createFBeta(int _expBool);
  void createFBeta2(int _expBool);
  void createFBeta3(int _expBool);
  void createFBetaMany(int _expBool);

  void createFBeta2r2(int _expBool);
  void createFExtended(int _expBool);
  void createFNone();


  // ********************* Coordinates and Spin ********************
  int           currentParticle;// The current particle (being moved)
  int           otherParticle;  // The proposed particle to interchange spin 
                                // with (if spin flip is allowed)
  double*       coorArray;      // Cartesian coordinates of all coors.
  double*       _coorArray;      
  double*       trialCoorArray; // Cartesian coordinates of trialCoor.
  int           spinFlip;       // Boolean, is spin being proposed flipped.

  int*          spinArray;      // Array containing the different spins.
  CoorSpinDiff* coors;          // The coordinate objects.
  CoorSpinDiff* trialCoor;      // The trial coordinate.
  void          computeTrialPosition(); // Proposes a move of the current
                                        // particle.
  double        coorStep();     // Random walk of one coordinate.
  void          proposeFlip();  // Proposes an interchange of electron spin
                                // between current particle and (randomly)
                                // generated other particle.
  void          acceptTrialPosition(); // Accept the proposed move.
  void          rejectTrialPosition(); // Rejects the proposed move.
  void          resetPtr();            // Sets coorArray = _coorArray.
  void          setPositionToOrigin(); // Sets all coordinates equal zero.


  // ****************** Inter-electronic distances ******************
  Distance*     distance;     // Object that keeps track of the distances.
  DistanceDiff* distanceDiff; // Keeps track of the differentiated distances.
  void          createDistanceAndDistanceDiff(); // Creates the above objects.


  // ************************ Spin-Factors *************************
  // Object to impose the different cusp conditions for like-spin and opposite-
  // spin electrons, 1/4 and 1/2 respectively. Is turned on and off by the 
  // value of quarterCusp (on(1) cusp=1/4 off (0)).
  SpinFactors* spinFactors;
  void         createSpinFactors();




 // ****************************************************************
 // *                     PUBLIC ALGORITHMS                        *
 // ****************************************************************
 public:
  // Constructor: Input file name, prosess number, and number of MPI
  // processes
  Domain(char* _initFileName, int _rank, int _size);


  // ***************************** MPI *****************************
  int           getRank()                  {return rank;}
  int           getSize()                  {return size;}


  // ******************* Input and Output Files ********************
  ofstream*     getOutputFile()            {return &outputFile;}
  char*         getFixedParamsUp()         {return fixedParamsUp;}
  char*         getFixedParamsDown()       {return fixedParamsDown;}


  // *********************** Standard Input ************************
  int           getNumDimensions()         {return numDimensions;}
  int           getNumThermalization()     {return numThermalization;}
  int           getNumCycles()             {return numCycles;}
  double        getStepLen()               {return stepLen;}

  int           getAllowSpinFlip()         {return allowSpinFlip;}
  string        getThermalizationType()    {return thermalizationType;}
  double        getSoughtAcceptance()      {return soughtAcceptance;}
  int           getNumberSeekAcceptance()  {return numberSeekAcceptance;}

  int           getVarySeekWithRank()      {return varySeekWithRank;}
  double        getVarySeekWithRankStep()  {return varySeekWithRankStep;}
  string        getVmcType()               {return vmcType;}
  int           getVarianceOptimization()  {return varianceOptimization;}

  double        getReferenceEnergy()       {return referenceEnergy;}
  double        getDeltaE()                {return deltaE;}
  int           getSetWeightToUnity()      {return setWeightToUnity;}
  void          changeStepLen(int numSteps, int acceptedSteps, 
			      double soughtAcceptance);

  void          setReferenceEnergy(double E) {referenceEnergy=E;}


  // ******************* Electron Configuration ********************
  int getNumParticles()   {return numParticles;}
  int getUp1s()           {return up1s;}
  int getUp2s()           {return up2s;}
  int getUp2px()          {return up2px;}

  int getUp2py()          {return up2py;}
  int getUp2pz()          {return up2pz;}
  int getDown1s()         {return down1s;}
  int getDown2s()         {return down2s;}

  int getDown2px()        {return down2px;}
  int getDown2py()        {return down2py;}
  int getDown2pz()        {return down2pz;}


  // *********************** Random Config ************************
  Ref<Random2> getRandomMetro() {return randomMetro;}
  Ref<Random2> getRandomMove()  {return randomMove;}


  // **************** Uni-Directional Configuration ***************
  void increaseNumCycles(); 
  void increaseNumThermalization(); 
  void reduceLocalArea();
  int  getUniDirectionalMovement(){return uniDirectionalMovement;}

  int  getNumberVmcRuns()         {return numberVmcRuns;}
  int  getNumberOfUniDirectionalMoves() {return numberOfUniDirectionalMoves;}
  void initUniDirectionalIndicator();
  int  isMovementUniDirectional(double* alphaParams, double* betaParams);


  // *************************** Other *****************************
  int           getNumVariations()         {return numVariations;}
  double        getH()                     {return h;}


  // *************** Slater-Determinant Parameters ****************
  string  getOrbitalType()         {return orbitalType;}
  int     getNumAlpha()            {return numAlpha;} 
  int     getNumSlaterVariations() {return numSlaterVariations;}
  double* getAlphaParam(int i)     {return &(alphaParam[i*numAlpha]);}

  int     getCentralSlater()       {return centerSlater;}
  void    setCentralAlphaParam(double* param);
  void    calculateAlphaParamArray();


  // ******************* Correlation Parameters *******************
  int        getNumBeta()             {return numBeta;} 
  double*    getBetaParam(int i)      {return &(betaParam[i*numBeta]);}
  void       setCentralBetaParam(double* param);
  int        getCentralCorrel()       {return centerCorrel;}

  fFunction* getF()                   {return f;}
  void       calculateBetaParamArray();
  int        getNumCorrelVariations() {return numCorrelVariations;}


  // ******************** Coordinates and Spin *********************
  double*       getCoorArray()             {return coorArray;}
  CoorSpinDiff* getCoors();
  CoorSpinDiff* getTrialCoor();
  int*          getSpinArray()             {return spinArray;}

  int*          getSpinFlip()              {return &spinFlip;}
  int           getOtherParticle()         {return otherParticle;}
  void          setCoorsToOrigon();
  void          setCurrentParticle(int __currentParticle);

  void          setOtherParticle(int __otherParticle);
  void          setToNextParticle();
  void          init();
  void          initVMC();

  void          suggestMove();
  void          acceptMove();
  void          rejectMove();
  void          acceptThermalizedMove();

  void          rejectThermalizedMove();


  // ****************** Inter-electronic distances *****************
  Distance*     getDistance()              {return distance;}
  DistanceDiff* getDistanceDiff()          {return distanceDiff;}


  // ************************ Spin-Factors *************************
  SpinFactors*  getSpinFactors()           {return spinFactors;}


  // ********************** Potential Energy ***********************
  double getNucleusElectronPotential();
  double getInterElectronicPot() {return distance->getPotential();}


  // ************************** Summary ****************************
  void    initSummary();
  void    Summary();
};

#endif
