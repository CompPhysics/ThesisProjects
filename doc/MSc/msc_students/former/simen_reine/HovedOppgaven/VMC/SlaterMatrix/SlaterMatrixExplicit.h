#ifndef SlaterMatrix_IS_INCLUDED
#define SlaterMatrix_IS_INCLUDED

#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

// ****************************************************************
// *                        SLATERMATRIX                          *
// ****************************************************************
class SlaterMatrix {

 protected:
  //Primary variables
  double* inverseMatrix;
  double* inverseMatrixCurrentRow;
  double* newColumn;
  double  ratio;
  int     dim;
  int     currentColumn;
  int     isAllocated;
#ifdef _DEBUG_
  int     ratioIsCalculated;
#endif

  //Iterating and secondary variables
  double* _inverseMatrix;
  double* _inverseMatrixCurrentRow;
  double* _newColumn;
  double  result;
  double  invRatio;

 public:
  SlaterMatrix();
  SlaterMatrix(int __dim);
  virtual        ~SlaterMatrix();
  void            init();
  void            redim(int __dim);
  void    setPtrToNewColumn(double* __newColumn);
  virtual void    setCurrentColumn(int __currentColumn);
  virtual void    setToNextColumn();
  virtual void    resetPtr();
  virtual void    calcRatio();
  virtual void    importNewColumn();
  void            importNewMatrixColumnwise(double* __matrix);
  virtual void    setMatrixToUnity();

  virtual void    summary();

  double&         getRatio()        {return ratio;}
  int             getDim()          {return dim;}
  double*         getInvMatrixPtr() {return inverseMatrix;}

 protected:
  virtual int     allocateMemory(int __dim);
  virtual void    deleteMemory();
#ifdef _DEBUG_
  int     checkAllocation();
  int     checkIfRatioCalculated();
#endif

};


// ****************************************************************
// *                      SLATERMATRIX_DET                        *
// ****************************************************************
class SlaterMatrix_Det : public SlaterMatrix {

 protected:
  double trialDeterminant;
  double determinant;

 public:
  SlaterMatrix_Det();
  SlaterMatrix_Det(int __dim);
  virtual         ~SlaterMatrix_Det();
  inline void      importNewColumn();
  void             setMatrixToUnity();
  inline void      calcRatio();

  virtual void     summary();

  double&          getDet()         {return determinant;}
  double&          getTrialDet()    {return trialDeterminant;}
};

/*
class SlaterMatrix_Explicit : public SlaterMatrix {

 protected:
  double* matrix;
  double* matrixCurrentColumn;

  double* _matrix;

 public:
  SlaterMatrix_Explicit();
  SlaterMatrix_Explicit(int __dim);
  virtual         ~SlaterMatrix_Explicit();
  inline void      setCurrentColumn(int __currentColumn);
  inline void      setToNextColumn();
  inline void      resetPtr();
  inline void      importNewColumn();
  void             setMatrixToUnity();

  virtual void     summary();

  double*          getMatrixPtr()              {return matrix;}
  double**         getMatrixCurrentColumnPtr() {return &matrixCurrentColumn;}

 protected:
  int              allocateMemory(int __dim);
  void             deleteMemory();

};
*/

#endif
