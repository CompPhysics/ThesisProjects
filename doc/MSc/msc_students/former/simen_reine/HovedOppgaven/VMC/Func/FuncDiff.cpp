#ifndef FuncCPP_IS_INCLUDED
#define FuncCPP_IS_INCLUDED

#include "Func.h"


template <class Param>
Func<Param>::Func(Param& _coordinate) {
  setCoordinate(_coordinate);
  calcH();
  numDimensions = 1;
}


template <class Param>
Func<Param>::Func() {
  calcH();
  numDimensions = 1;
}

template <class Param>
Func<Param>::~Func() {
  delete[] function[0];
  //delete[] function;
  delete[] values;
}


template <class Param>
void Func<Param>::init(Param& _coordinate) {
  setCoordinate(_coordinate);
  init();
}


template <class Param>
void Func<Param>::init() {
  function = new SingleParticleFunc<Param>*;
  values   = new double[3];
}


template <class Param>
void Func<Param>::setCoordinate(Param& _coordinate) {
  coordinate = _coordinate;
}


template <class Param>
void Func<Param>::calcH() {
  h           = 1e-4;
  h_half      = 0.5*h;
  h_double    = 2.0*h;
  h_inv       = 1./h;
  h_inv2      = 1./(h*h);
  h_invdouble = 1./(2.*h);
}


template <class Param>
void Func<Param>::calcValueCenter() {
  (*values) = (**function)(coordinate());
}


template <class Param>
void Func<Param>::calcValueCenter(Param& _coordinate) {
  (*values) = (**function)(_coordinate);
}


template <class Param>
void Func<Param>::calcValueSides() {
  values[1] = (**function)(coordinate(), 0);
  values[2] = (**function)(coordinate(), 1);
}


template <class Param>
void Func<Param>::calcValueSides(Param& _coordinate) {
  values[1] = (**function)(_coordinate, 0);
  values[2] = (**function)(_coordinate, 1);
}


template <class Param>
double Func<Param>::valuePt() {
  return (*result = values[0]);
}


template <class Param>
double Func<Param>::diff() {
  *diffResult = (values[1] - values[2])*h_invdouble;
  return *diffResult;
}


template <class Param>
double Func<Param>::ddiff() {
  *ddiffResult = (values[1] + values[2] - 2*values[0])*h_inv2;
  return *ddiffResult;
}


template <class Param>
void Func<Param>::attachResult(double* __result) {
  result = __result;
}


template <class Param>
void Func<Param>::attachDiffResult(double* __diffResult) {
  diffResult = __diffResult;
}


template <class Param>
void Func<Param>::attachDdiffResult(double* __ddiffResult) {
  ddiffResult = __ddiffResult;
}


template <class Param>
void Func<Param>::summary() {
  cerr << "\n\n--------------------------------------------------\n";
  cerr << "Values:\n";
  cerr << " Center: " << values[0] << endl;
  cerr << " Minus:  " << values[1] << endl;
  cerr << " Plus:   " << values[2] << endl;
  cerr << "result:      "  << *result << endl;
  cerr << "diffResult:  "  << *diffResult << endl;
  cerr << "ddiffResult: "  << *ddiffResult << endl;  
}


/************************************************************/


template <class Param>
FuncSet<Param>::FuncSet(Param& _coordinate, int _len) : Func<Param>(_coordinate), len(_len) {}


template <class Param>
FuncSet<Param>::FuncSet() : Func<Param>() {}

template <class Param>
FuncSet<Param>::~FuncSet() {
  for (int i=0; i<len; i++)
    delete[] function[i];
  //delete[] function;
  delete[] values;
}

template <class Param>
void FuncSet<Param>::init(Param& _coordinate, int _len) {
  len = _len;
  Func<Param>::init(_coordinate);
}


template <class Param>
void FuncSet<Param>::init() {
  values   = new double[3*len];
  valuesP  = values + len;
  valuesM  = values + len*2;
  function = new SingleParticleFunc<Param>*[len];
}


template <class Param>
void FuncSet<Param>::calcValueCenter() {
  for (int i = 0; i < len; i++)
    values[i] = function[i][0](coordinate());
}


template <class Param>
void FuncSet<Param>::calcValueCenter(Param& _coordinate) {
  for (int i = 0; i < len; i++)
    values[i] = function[i][0](_coordinate);
}


template <class Param>
void FuncSet<Param>::calcValueSides() {

  for (int i = 0; i < len; i++) {
    valuesP[i] = function[i][0](coordinate(), 0);
    valuesM[i] = function[i][0](coordinate(), numDimensions );
  }
}


template <class Param>
void FuncSet<Param>::calcValueSides(Param& _coordinate) {

  for (int i = 0; i < len; i++) {
    valuesP[i] = function[i][0](_coordinate, 0);
    valuesM[i] = function[i][0](_coordinate, numDimensions );
  }
}


template <class Param>
double FuncSet<Param>::valuePt() {
  for (int i = 0; i < len; i++)
    result[i] = values[i];

  return 0;
}


template <class Param>
double FuncSet<Param>::diff() {
  for (int i = 0; i < len; i++)
    diffResult[i] = (valuesP[i]-valuesM[i])*h_invdouble;
  return 0;
}


template <class Param>
double FuncSet<Param>::ddiff() {
  for (int i = 0; i < len; i++)
    ddiffResult[i] = (valuesM[i]+valuesP[i]-2*values[i])*h_inv2;

  return 0;
}


template <class Param>
void FuncSet<Param>::summary() {
  string returnStr;
  char   buffer[50];

  cerr << "\n\n--------------------------------------------------\n";
  cerr << "Values:\n\n";
  for (int i = 0; i < len; i++) {
    snprintf(buffer, 50, "%5.5f  %5.5f  %5.5f\n", values[i], values[i+len], values[i+2*len]);
    returnStr += buffer;
  }
  cerr << endl;

  returnStr = "";
  cerr << "Result:\n";
  for (int i = 0; i < len; i++) {
    snprintf(buffer, 50, "%5.5f", result[i]);
    returnStr += buffer;
    returnStr += "\n";
  }
  cerr << returnStr;
}


/************************************************************/


template <class Param>
FuncSetMultivar<Param>::FuncSetMultivar(Param& _coordinate, int _len) : FuncSet<Param>(_coordinate, _len) {
  numVar     = coordinate.getLen();
  len_double = 2*_len;
}


template <class Param>
FuncSetMultivar<Param>::FuncSetMultivar() : FuncSet<Param>() {
}

template <class Param>
FuncSetMultivar<Param>::~FuncSetMultivar() {
  for (int i=0; i<len; i++)
    delete[] function[i];
  //delete[] function;
  delete[] values;
}

template <class Param>
void FuncSetMultivar<Param>::init() {
  values   = new double[(1+2*numVar)*len];
  valuesP  = values + len;
  valuesM  = values + len*2;
  function = new SingleParticleFunc<Param>*[len];
}


template <class Param>
void FuncSetMultivar<Param>::init(Param& _coordinate, int _len) {
  numVar     = _coordinate.getLen();
  len_double = 2*_len;
  FuncSet<Param>::init(_coordinate, _len);
}


template <class Param>
void FuncSetMultivar<Param>::calcValueSides() {
  //  for (int i = (numVar-1); i >= 0; i--) {
  for (int i = 0; i < numVar; i++) {

    int iNum = i + numVar;
    for (int j = 0; j < len; j++) {
      valuesP[j] = function[j][0](coordinate(), i);
      valuesM[j] = function[j][0](coordinate(), iNum );
    }
    valuesP += (len_double);
    valuesM += (len_double);
  }
  valuesP = values + len;
  valuesM = values + len_double;
}


template <class Param>
void FuncSetMultivar<Param>::calcValueSides(Param& _coordinate) {
  //  for (int i = (numVar-1); i >= 0; i--) {
  for (int i = 0; i < numVar; i++) {
    int iNum = i + numVar;
    for (int j = 0; j < len; j++) {
      valuesP[j] = function[j][0](_coordinate, i);
      valuesM[j] = function[j][0](_coordinate, iNum );
    }
    valuesP += (len_double);
    valuesM += (len_double);
  }
  valuesP = values + len;
  valuesM = values + len_double;
}


template <class Param>
double FuncSetMultivar<Param>::diff() {
  for (int v = (numVar-1); v > 0; v--) {
    FuncSet<Param>::diff();
    diffResult += len;
    valuesM    += (len_double);
    valuesP    += (len_double);
  }
  FuncSet<Param>::diff();

  diffResult = _diffResult;
  valuesP = values + len;
  valuesM = values + len_double;  


  return 0;
}


template <class Param>
double FuncSetMultivar<Param>::diff(int v) {
  //Calculate the first derivative with respect to the variable indexed
  //v of all the functions and place the result in the result array (!!!)
#ifdef _DEBUG_
  if (v >= numVar) {
    cerr << "\nCannot calculate first derivative. Variable index out of bound.";
    exit(1);
  }
#endif

  diffResult  = result;
  valuesM    += (len_double*v);
  valuesP    += (len_double*v);
  FuncSet<Param>::diff();

  diffResult = _diffResult;
  valuesP    = values + len;
  valuesM    = values + len_double;  

  return 0;
}


template <class Param>
double FuncSetMultivar<Param>::ddiff() {
  for (int i = 0; i < len; i++)
    ddiffResult[i] = -(2*numVar*values[i]);
  for (int v = (numVar-1); v > 0; v--) {
    for (int i = 0; i < len; i++) 
      ddiffResult[i] += (valuesM[i] + valuesP[i]);
    valuesP += (len_double);
    valuesM += (len_double);
  }
  for (int i = 0; i < len; i++) {
    ddiffResult[i] += (valuesM[i] + valuesP[i]);  
    ddiffResult[i] *= h_inv2;
  }

  valuesP = values + len;
  valuesM = values + len_double;

  return 0;
}


template <class Param>
void FuncSetMultivar<Param>::attachDiffResult(double* __diffResult) {
  FuncSet<Param>::attachDiffResult(__diffResult);
  _diffResult = __diffResult;
}


template <class Param>
void FuncSetMultivar<Param>::summary() {
  string returnStr;
  char   buffer[50];

  cerr << "\n\n--------------------------------------------------\n";
  cerr << "Values:\n\n";
  for (int i = 0; i < len; i++) {
    snprintf(buffer, 50, "%5.5f  ", values[i]);
    returnStr += buffer;
    for (int j = 0; j < numVar; j++) {
      snprintf(buffer, 50, "  %5.5f  %5.5f", valuesM[i+j*2*len], valuesP[i+j*2*len]);
      returnStr += buffer;
    }
    returnStr += "\n";
  }
  cerr << returnStr << endl;

  returnStr = "";
  cerr << "Result:\n";
  for (int i = 0; i < len; i++) {
    snprintf(buffer, 50, "%5.5f", result[i]);
    returnStr += buffer;
    returnStr += "\n";
  }

  cerr << returnStr;

}

/************************************************************/
/*

template <class Param>
FuncDiff<Param>::FuncDiff(Param& _coordinate, int _len) : FuncSetMultivar<Param>(_coordinate, _len) {
}


template <class Param>
FuncDiff<Param>::FuncDiff() : FuncSetMultivar<Param>() {
}

template <class Param>
FuncDiff<Param>::~FuncDiff() {
  for (int i=0; i<len; i++)
    delete[] function[i];
  //delete[] function;
  //delete[] values;
  for (int i=0; i<len; i++)
    delete[] dfunction[i];
  //delete[] dfunction;
  for (int i=0; i<len; i++)
    delete[] ddfunction[i];
  //delete[] ddfunction;
}

template <class Param>
void FuncDiff<Param>::init() {
  function   = new SingleParticleFunc<Param>*[len];
  dfunction  = new SingleParticleFunc<Param>*[len*numVar];
  ddfunction = new SingleParticleFunc<Param>[len];
}


template <class Param>
void FuncDiff<Param>::init(Param& _coordinate, int _len) {
  FuncSetMultivar<Param>::init(_coordinate, _len);
}


template <class Param>
void FuncDiff<Param>::calcValueCenter() {
}


template <class Param>
void FuncDiff<Param>::calcValueCenter(Param& _coordinate) {
}


template <class Param>
void FuncDiff<Param>::calcValueSides() {
}


template <class Param>
void FuncDiff<Param>::calcValueSides(Param& _coordinate) {
}


template <class Param>
double FuncDiff<Param>::valuePt() {
  for (int i = 0; i < len; i++)
    result[i] = function[i][0](coordinate);

  return 0;
}


template <class Param>
double FuncDiff<Param>::valuePt(Param& _coordinate) {
  for (int i = 0; i < len; i++)
    result[i] = function[i][0](_coordinate);

  return 0;
}


template <class Param>
double FuncDiff<Param>::diff() {
  int iMax = numVar*len;
  for (int i = 0; i < iMax; i++)
    diffResult[i] = dfunction[i][0](coordinate());

  return 0;
}


template <class Param>
double FuncDiff<Param>::diff(Param& _coordinate) {
  int iMax = numVar*len;
  for (int i = 0; i < iMax; i++)
    diffResult[i] = dfunction[i][0](_coordinate);

  return 0;
}


template <class Param>
double FuncDiff<Param>::ddiff() {
  for (int i = 0; i < len; i++)
    ddiffResult[i] = ddfunction[i][0](coordinate());

  return 0;
}


template <class Param>
double FuncDiff<Param>::ddiff(Param& _coordinate) {
  for (int i = 0; i < len; i++)
    ddiffResult[i] = ddfunction[i][0](_coordinate);

  return 0;
}
*/


#endif
