#define TRUE 1
#define FALSE 0
#define NULL 0

#include "SlaterMatrix.h"

SlaterMatrix::SlaterMatrix() {

  isAllocated       = FALSE;
#ifdef _DEBUG_
  ratioIsCalculated = FALSE;
#endif

}

SlaterMatrix::SlaterMatrix(int __dim) {
  isAllocated       = FALSE;
#ifdef _DEBUG_
  ratioIsCalculated = FALSE;
#endif

  redim(__dim);
}


SlaterMatrix::~SlaterMatrix() {
  deleteMemory();
}


void SlaterMatrix::init() {
  setMatrixToUnity();
  setCurrentColumn(0);
}


void SlaterMatrix::redim(int __dim) {
  deleteMemory();
  allocateMemory(__dim);
  init();
}


void SlaterMatrix::setPtrToNewColumn(double* __newColumn) {
#ifdef _DEBUG_
  if (newColumn != __newColumn)
    ratioIsCalculated = FALSE;
#endif

  newColumn = __newColumn;
}


void SlaterMatrix::setCurrentColumn(int __currentColumn) {
#ifdef _DEBUG_
  if (checkAllocation() == 1) exit(1);
  if (currentColumn != __currentColumn)
    ratioIsCalculated = FALSE;
#endif
  
  currentColumn = __currentColumn;
  inverseMatrixCurrentRow = (inverseMatrix + dim*currentColumn);
}


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


void SlaterMatrix::resetPtr() {
  inverseMatrixCurrentRow = inverseMatrix;
}


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


int SlaterMatrix::allocateMemory(int __dim) {
  if (isAllocated) return 1;
  
  inverseMatrix = new double[__dim*__dim+1];
  dim = __dim;
  
  isAllocated = TRUE;
  
  return 0;
}


void SlaterMatrix::deleteMemory() {
  if (isAllocated)
    delete[] inverseMatrix;
  isAllocated = FALSE;
}


#ifdef _DEBUG_

int SlaterMatrix::checkAllocation() {
  if (!isAllocated) {
    cerr << "\nNo memory allocated for inverse slater matrix.";
    return 1;
  }
  
  return 0;
}


int SlaterMatrix::checkIfRatioCalculated() {
  if (!ratioIsCalculated) {
    cerr << "\nDeterminant ratio not yet calculated.";
    return 1;
  }
  
  return 0;
}

#endif


/************************************************************/


SlaterMatrix_Det::SlaterMatrix_Det() : SlaterMatrix() {
}


SlaterMatrix_Det::SlaterMatrix_Det(int __dim) : SlaterMatrix(__dim){
}


SlaterMatrix_Det::~SlaterMatrix_Det() {
  deleteMemory();
}


void SlaterMatrix_Det::calcRatio() {
  SlaterMatrix::calcRatio();
  trialDeterminant = ratio*determinant;
}


void SlaterMatrix_Det::importNewColumn() {
  SlaterMatrix::importNewColumn();
  determinant = trialDeterminant;
}


void SlaterMatrix_Det::setMatrixToUnity() {
  SlaterMatrix::setMatrixToUnity();
  determinant = 1.;
}


void SlaterMatrix_Det::summary() {
  SlaterMatrix::summary();
  cerr << "\nDeterminant: " << determinant << endl;
}


/************************************************************/

/*
SlaterMatrix_Explicit::SlaterMatrix_Explicit() : SlaterMatrix() {
}


SlaterMatrix_Explicit::SlaterMatrix_Explicit(int __dim) : SlaterMatrix(__dim){
}


SlaterMatrix_Explicit::~SlaterMatrix_Explicit() {
  deleteMemory();
}


void SlaterMatrix_Explicit::setCurrentColumn(int __currentColumn) {
#ifdef _DEBUG_
  if (checkAllocation() == 1) exit(1);
  if (currentColumn != __currentColumn)
    ratioIsCalculated = FALSE;
#endif
  
  currentColumn = __currentColumn;
  inverseMatrixCurrentRow = (inverseMatrix + dim*currentColumn);
  matrixCurrentColumn     = (matrix        + dim*currentColumn);
}


void SlaterMatrix_Explicit::setToNextColumn() {
#ifdef _DEBUG_
  if (checkAllocation() == 1) exit(1);
#endif
  
  if (++currentColumn == dim) {
    currentColumn = 0;
    resetPtr();
  }
  else {
    inverseMatrixCurrentRow += dim;
    matrixCurrentColumn     += dim;
  }
  
#ifdef _DEBUG_
  ratioIsCalculated = FALSE;
#endif
}


void SlaterMatrix_Explicit::resetPtr() {
  SlaterMatrix::resetPtr();
  matrixCurrentColumn = matrix;
}


void SlaterMatrix_Explicit::importNewColumn() {
  
  SlaterMatrix::importNewColumn();
  
#ifdef _DEBUG_
  if (checkAllocation() == 1) {
    cerr << "\nCannot import new column.";
    exit(1);
  }
#endif
  
  _matrix    = matrixCurrentColumn;
  _newColumn = newColumn;
  for (int i = 0; i < dim; i++)
    *_matrix++ = *_newColumn++;
  
}


void SlaterMatrix_Explicit::setMatrixToUnity() {
#ifdef _DEBUG_
  if (checkAllocation() == 1) {
    cerr << "\nCannot set matrix to unity";
    exit(1);
  }
#endif
  
  _inverseMatrix = inverseMatrix;
  _matrix        = matrix;
  
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      *_inverseMatrix++  =  *_matrix++  =  ( j==i );
  
#ifdef _DEBUG_
  ratioIsCalculated = FALSE;
#endif
}


void SlaterMatrix_Explicit::summary() {
  SlaterMatrix::summary();
  
  string result;
  char   buffer[50];
  
#ifdef _DEBUG_
  if (checkAllocation() == 1) exit(1);
#endif
  
  cerr << "\nSlater matrix:\n\n";
  
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      snprintf(buffer, 50, "%10.5f", matrix[i + j*dim]);
      result += buffer;
      if (j != (dim-1))
        result += ", ";
    }
    result += "\n";
  }
  
  cerr << result;
  
}


int SlaterMatrix_Explicit::allocateMemory(int __dim) {
  if (isAllocated) return 1;
  
  inverseMatrix = new double[__dim*__dim+1];
  matrix = new double[__dim*__dim+1];
  
  dim = __dim;
  
  isAllocated = TRUE;
  
  return 0;
}


void SlaterMatrix_Explicit::deleteMemory() {
  if (isAllocated) {
    delete[] inverseMatrix;
    delete[] matrix;
  }
  isAllocated = FALSE;
}
*/
