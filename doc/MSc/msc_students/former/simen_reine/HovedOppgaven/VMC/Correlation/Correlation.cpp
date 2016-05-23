#include "Correlation.h"

//*****************************************************************
//*                           CORREL                              *
//*****************************************************************
//
//*************************** Correl *******************************
Correl::Correl() {}


//*************************** attach ******************************
void Correl::attach(Domain& _domain) {
  domain          = &_domain;
  numParticles    = domain->getNumParticles();
  numDimensions   = domain->getNumDimensions();
  Nm1             = numParticles-1;
  f               = domain->getF();

  spinFactors     = domain->getSpinFactors();
  spinFactorMatrix= spinFactors->getMatrix();

  gradJastrow                  = new double[numParticles*numDimensions];

  distance        = domain->getDistance();
  distanceDiff    = domain->getDistanceDiff();
  h               = domain->getH();
  hDoubleInverse  = 1./2./h;
  hSquaredInverse = 1./h/h;
}


//*********************** attachFFunction *************************
void Correl::attachFFunction(fFunction* _f) {
  f = _f;
  expBool=f->getExpBool();
}


//****************** createJastrowAndJastrowDiff ******************
// Create jastrow and jastrowDiff 
// NB: This must be preformed AFTER create/attachDistanceAndDistanceDiff 
//     AND AFTER createFBeta/attachFFunction
void Correl::createJastrowAndJastrowDiff() {
  // Create jastrow and jastrowDiff
  jastrow           = new Jastrow[1];
  jastrow->attach(domain, f);
  jastrowDiff       = new JastrowDiff[ numDimensions ]; 
  for (int i=0; i<numDimensions; i++)
      jastrowDiff[i].attach(domain, f, i);
  jastrowMatrix.resize(numParticles, numParticles);
  jastrowMatrixPlus  = new Array<double, 2>[3];
  jastrowMatrixMinus = new Array<double, 2>[3];
  for (int i=0; i<numDimensions; i++) {
    jastrowMatrixPlus[i].resize(numParticles, numParticles);
    jastrowMatrixMinus[i].resize(numParticles, numParticles);
  }
}



//***********************setToNextParticle ************************
void Correl::setToNextParticle() {
  currentParticle++;
  if (currentParticle==numParticles) currentParticle=0;
  jastrow->setToNextParticle();
  for (int i=0; i<numDimensions; i++)
    jastrowDiff[i].setToNextParticle();
}


//********************* setCurrentParticle ************************
void  Correl::setCurrentParticle(int _currentParticle){
  currentParticle=_currentParticle;
  jastrow->setCurrentParticle(currentParticle);
  for (int i=0; i<numDimensions; i++)
    jastrowDiff[i].setCurrentParticle(currentParticle);
}


//************************* suggestMove ***************************
void Correl::suggestMove() {
  jastrow->suggestMove();
  calculateRatio();
}


//******************* initializeThermalization ********************
void Correl::initializeThermalization() {
  jastrow->initialize();
  currentParticle=0;
}


//******************** acceptThermalizedMove **********************
void Correl::acceptThermalizedMove() {
  jastrow->acceptMove();
}


//******************** rejectThermalizedMove **********************
void Correl::rejectThermalizedMove() {
  jastrow->rejectMove();
}


//************************ initializeVMC **************************
void Correl::initializeVMC() {
  for (int i=0; i<numDimensions; i++)
    jastrowDiff[i].initialize();
  spinFactors->calculateMatrix();
  setCurrentParticle(0);
}


//************************** acceptMove ***************************
void Correl::acceptMove() {
  for (int i=0; i<numDimensions; i++)
    jastrowDiff[i].suggestMove();
  acceptThermalizedMove();
  for (int i=0; i<numDimensions; i++)
    jastrowDiff[i].acceptMove();
}


//************************** rejectMove ***************************
void Correl::rejectMove() {
  rejectThermalizedMove();
  for (int i=0; i<numDimensions; i++)
    jastrowDiff[i].rejectMove();
}


//*************** calculateGradAndLaplacianRatios *****************
void Correl::calculateGradAndLaplacianRatios() {

  laplacian = 0.0;
  int gjIndex = 0;
  double J = 1;
  if (!expBool) J = (*jastrow)();
  for (int i=0; i<numParticles; i++) {
    for (int k=0; k<numDimensions; k++) {
      jastrowMatrixPlus[k]  = jastrowDiff[k].getJastrowMatrixPlus();
      jastrowMatrixMinus[k] = jastrowDiff[k].getJastrowMatrixMinus();

      double gradPlus = 0, gradMinus = 0;
      int t=Nm1;
      int l=i-1;
      for (int j=0; j<i; j++) {
	gradPlus  += spinFactorMatrix[l]*jastrowMatrixPlus[k](j,i);
	gradMinus += spinFactorMatrix[l]*jastrowMatrixMinus[k](j,i);
	l += --t;
      }
      for (int j=i+1; j<numParticles; j++) {
	l++;
	gradPlus  += spinFactorMatrix[l]*jastrowMatrixPlus[k](j,i);
	gradMinus += spinFactorMatrix[l]*jastrowMatrixMinus[k](j,i);
      }
      gradJastrow[gjIndex++] = (gradMinus-gradPlus)*hDoubleInverse/J;
      laplacian += gradPlus + gradMinus;
    }
  }

  jastrowMatrix = jastrow->getJastrowMatrix();

  double grad = 0;
  for (int i=0; i<numParticles; i++) {
    int k=Nm1;
    int l=i-1;
    for (int j=0; j<i; j++) {
      grad  += spinFactorMatrix[l]*jastrowMatrix(j,i);
      l += --k;
    }
    for (int j=i+1; j<numParticles; j++) {
      l++;
      grad  += spinFactorMatrix[l]*jastrowMatrix(j,i);
    }
  }

  laplacian -= numDimensions*2.*grad;
  laplacian *= hSquaredInverse/J;

  if (expBool)
    for (int i=0; i<numParticles*numDimensions; i++) 
        laplacian += sqr( gradJastrow[i] );

}


//************************ calculateRatio *************************
// Only applicabel prior to accepting move
void Correl::calculateRatio() {
  if (expBool) { ratio = exp( (*jastrow).getDifference() );}
  else { ratio =  ((*jastrow)()-(*jastrow).getDifference())/(*jastrow)(); }
}


//********************* calculateCorrelation **********************
// To be used to get direct access to the value of the correlation
void Correl::calculateCorrelation() {
  if (expBool) {correlation = exp( (*jastrow)() );}
  else {correlation = (*jastrow)(); }  
}


//*************************** operator ****************************
// Only applicable after accepting move
double Correl::operator()() {
  if (expBool) {return correlation = exp( (*jastrow)() );}
  else {return correlation =  (*jastrow)(); }
}




