#define TRUE 1
#define FALSE 0
#define NULL 0

#include "SlaterMatrix.h"

// ****************************************************************
// *                        SLATERMATRIX                          *
// ****************************************************************
//
// ************************ SlaterMatrix **************************
SlaterMatrix::SlaterMatrix() {
#ifdef _DEBUG_
  ratioIsCalculated = FALSE;
#endif

  isAllocated       = FALSE;
}

// ************************ SlaterMatrix **************************
SlaterMatrix::SlaterMatrix(int __dim) {
  isAllocated       = FALSE;
#ifdef _DEBUG_
  ratioIsCalculated = FALSE;
#endif

  redim(__dim);
}


// ************************ ~SlaterMatrix *************************
SlaterMatrix::~SlaterMatrix() {
  deleteMemory();
}


// **************************** init ******************************
void SlaterMatrix::init() {
  setMatrixToUnity();
  setCurrentColumn(0);
}


// **************************** redim *****************************
void SlaterMatrix::redim(int __dim) {
  deleteMemory();
  allocateMemory(__dim);
  init();
}


// ********************** setPtrToNewColumn ***********************
void SlaterMatrix::setPtrToNewColumn(double* __newColumn) {
#ifdef _DEBUG_
  if (newColumn != __newColumn)
    ratioIsCalculated = FALSE;
#endif

  newColumn = __newColumn;
}


// ********************** setCurrentColumn ************************
void SlaterMatrix::setCurrentColumn(int __currentColumn) {
#ifdef _DEBUG_
  if (checkAllocation() == 1) exit(1);
  if (currentColumn != __currentColumn)
    ratioIsCalculated = FALSE;
#endif
  
  currentColumn = __currentColumn;
  inverseMatrixCurrentRow = (inverseMatrix + dim*currentColumn);
}


// *********************** setToNextColumn ************************
void SlaterMatrix::setToNextColumn() {
#ifdef _DEBUG_
  if (checkAllocation() == 1) exit(1);
#endif
  
  if (++currentColumn == dim) {
    currentColumn = 0;
    resetPtr();
  }
  else
    inverseMatrixCurrentRow += dim;
  
#ifdef _DEBUG_
  ratioIsCalculated = FALSE;
#endif
}


// *************************** resetPtr ***************************
void SlaterMatrix::resetPtr() {
  inverseMatrixCurrentRow = inverseMatrix;
}


// ************************** calcRatio ***************************
void SlaterMatrix::calcRatio() {
#ifdef _DEBUG_
  if (checkAllocation() == 1) {
    cerr << "\nCannot calculate determinant ratio.";
    exit(1);
  }
#endif
  
  _newColumn               = newColumn;
  _inverseMatrixCurrentRow = inverseMatrixCurrentRow;
  
  ratio = 0;
  for (int i = 0; i < dim; i++)
    ratio += ( (*_newColumn++) * (*_inverseMatrixCurrentRow++));
  
#ifdef _DEBUG_
  if (ratio == 0)
    cerr << "\nRatio = 0.\n If new column is imported, slater matrix will not beinvertible.\n";
  ratioIsCalculated = TRUE;
#endif
}


// *********************** importNewColumn ************************
void SlaterMatrix::importNewColumn() {
#ifdef _DEBUG_
  if (checkAllocation() == 1) return;
  if (checkIfRatioCalculated() == 1) {
    cerr << "\nCannot update inverse matrix. Determinant ratio not calculated.";
    exit(1);
  }
#endif
  
  _inverseMatrix           = inverseMatrix;
  _inverseMatrixCurrentRow = inverseMatrixCurrentRow + dim;
  invRatio = 1/ratio;
  
  for (int i = 0; i < dim; i++)
    if (i != currentColumn) {
      result     = 0;
      _newColumn = newColumn;
      
      for (int j = 0; j < dim; j++)
        result += ((*_newColumn++) * (*_inverseMatrix++));
      
      result *= invRatio;
      
      _inverseMatrix           -= dim;
      _inverseMatrixCurrentRow -= dim;
      for (int j = 0; j < dim; j++)
        (*_inverseMatrix++) -= ((*_inverseMatrixCurrentRow++) * result);
    }
    else
      _inverseMatrix += dim;
  
  _inverseMatrixCurrentRow -= dim;
  for (int i = 0; i < dim; i++)
    (*_inverseMatrixCurrentRow++) *= invRatio;
  
#ifdef _DEBUG_
  ratioIsCalculated = FALSE;
#endif
  
}


// ****************** importNewMatrixColumnwise *******************
void SlaterMatrix::importNewMatrixColumnwise(double* __matrix) {
#ifdef _DEBUG_
  if (checkAllocation() == 1) {
    cerr << "\nCannot import new matrix.";
    exit(1);
  }
#endif
  
  double* __newColumn     = newColumn;
  int     __currentColumn = currentColumn;
  
  setCurrentColumn(0);
  for (int i = 0; i < dim; i++) {
    setPtrToNewColumn(__matrix);
    calcRatio();
    importNewColumn();
    setToNextColumn();
    __matrix += dim;
  }
  
  setPtrToNewColumn(__newColumn);
  setCurrentColumn(__currentColumn);
  
#ifdef _DEBUG_
  ratioIsCalculated = FALSE;
#endif
}


// ********************** setMatrixToUnity ************************
void SlaterMatrix::setMatrixToUnity() {
#ifdef _DEBUG_
  if (checkAllocation() == 1) {
    cerr << "\nCannot set matrix to unity.";
    exit(1);
  }
#endif

  _inverseMatrix = inverseMatrix;
  
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      *_inverseMatrix++ = ( j==i );
  
#ifdef _DEBUG_
  ratioIsCalculated = FALSE;
#endif
}


// *************************** summary ****************************
void SlaterMatrix::summary() {
  string result;
  char   buffer[50];
  
#ifdef _DEBUG_
  if (checkAllocation() == 1)
    exit(1);
#endif
  
  cerr << "\n--------------------------------------------------";
  cerr << "\nDimension: " << dim << "\n";
  cerr << "\nInverse slater matrix:\n\n";
  
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      snprintf(buffer, 50, "%10.5f", inverseMatrix[i*dim + j]);
      result += buffer;
      if (j != (dim-1))
        result += ", ";
    }
    result += "\n";
  }
  
  cerr << result;
  
}


// *********************** allocateMemory *************************
int SlaterMatrix::allocateMemory(int __dim) {
  if (isAllocated) return 1;
  
  inverseMatrix = new double[__dim*__dim+1];
  dim = __dim;
  
  isAllocated = TRUE;
  
  return 0;
}


// ************************ deleteMemory **************************
void SlaterMatrix::deleteMemory() {
  if (isAllocated)
    delete[] inverseMatrix;
  isAllocated = FALSE;
}


#ifdef _DEBUG_


// *********************** checkAllocation ************************
int SlaterMatrix::checkAllocation() {
  if (!isAllocated) {
    cerr << "\nNo memory allocated for inverse slater matrix.";
    return 1;
  }
  
  return 0;
}


// ******************** checkIfRatioCalculated ********************
int SlaterMatrix::checkIfRatioCalculated() {
  if (!ratioIsCalculated) {
    cerr << "\nDeterminant ratio not yet calculated.";
    return 1;
  }
  
  return 0;
}

#endif


// ****************************************************************
// *                      SLATERMATRIX_DET                        *
// ****************************************************************
//
// ********************** SlaterMatrix_Det ************************
SlaterMatrix_Det::SlaterMatrix_Det() : SlaterMatrix() {
}


// ********************** SlaterMatrix_Det ************************
SlaterMatrix_Det::SlaterMatrix_Det(int __dim) : SlaterMatrix(__dim){
}


// ********************** ~SlaterMatrix_Det ***********************
SlaterMatrix_Det::~SlaterMatrix_Det() {
  deleteMemory();
}


// ************************** calcRatio ***************************
void SlaterMatrix_Det::calcRatio() {
  SlaterMatrix::calcRatio();
  trialDeterminant = ratio*determinant;
}


// *********************** importNewColumn ************************
void SlaterMatrix_Det::importNewColumn() {
  SlaterMatrix::importNewColumn();
  determinant = trialDeterminant;
}


// ********************** setMatrixToUnity ************************
void SlaterMatrix_Det::setMatrixToUnity() {
  SlaterMatrix::setMatrixToUnity();
  determinant = 1.;
}


// *************************** summary ****************************
void SlaterMatrix_Det::summary() {
  SlaterMatrix::summary();
  cerr << "\nDeterminant: " << determinant << endl;
}
