#ifndef SpinFactors_IS_INCLUDED
#define SpinFactors_IS_INCLUDED

#include <iostream>
using namespace std;

#ifndef dotProduct_IS_INCLUDED
#define dotProduct_IS_INCLUDED

inline double dotProduct(double* array1, double* array2, int dim) {
//double dotProduct(double* array1, double* array2, int dim) {
  double* a1 = array1 - 1;
  double* a2 = array2 - 1;
  double product = 0;
  for (int i=0; i<dim; i++)
    product += (*++a1)*(*++a2);
  return product;
}
#endif

class SpinFactors {

protected:
  int     numParticles;
  int     Nmatrix;

  int     currentParticle;
  int     otherParticle;

  int*    spinArray;
  int*    _spinArray;
  int*    currentSpin;
  int*    otherSpin;

  double* newFactor;
  double* oldFactor;
  double* otherDifference;
  double* matrix;
  double* _newFactor;
  double* _oldFactor;
  double* _otherDifference;
  double* _matrix;

  int     quarterCusp;
  int     spinFlip;

  void    dumpNewFactor();
  void    dumpOldFactor();
  void    dumpOtherDifference();
  void    setNewFactor(double value);
public:
  SpinFactors() {}
  void    allocate(int _numParticles, int* _spinArray);
  void    init();
  // Call this routine if no spins are flipped 
  // (OR the two spins are equal) 
  void    calculateNoFlipFactors();
  // Call this routine ONLY if the two spins differ
  void    calculateFlipFactors();
  void    setToNextParticle();
  void    setCurrentParticle(int _currentParticle) 
    {currentParticle = _currentParticle;}
  void    setOtherParticle(int _otherParticle) 
    {otherParticle= _otherParticle;}
  void    setSpinFlip(int _spinFlip)       {spinFlip = _spinFlip;}
  void    setQuarterCusp(int _quarterCusp) {quarterCusp = _quarterCusp;}
  double* getNewFactor()       {return newFactor;}
  double* getOldFactor()       {return oldFactor;}
  double* getOtherDifference() {return otherDifference;}
  double* getMatrix()          {return matrix;}
  int*    getSpinFlip()        {return &spinFlip;}
  int*    getOtherParticle()   {return &otherParticle;}
  void    calculateMatrix();
  void    summaryNoFlip();
  void    summaryFlip();
  void    summayMatrix();
};

#endif
