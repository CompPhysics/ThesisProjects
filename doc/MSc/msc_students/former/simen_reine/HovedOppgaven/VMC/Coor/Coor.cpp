#include "Coor.h"

// ****************************************************************
// *                            COOR                              *
// ****************************************************************
Coor::Coor(double* __x, int _len) : x(__x), len(_len) {
  resetPtr();
}

void Coor::attach(double* __x, int _len) {
  x   = __x;
  len = _len;
  resetPtr();
}


// ****************************************************************
// *                           COORR                              *
// ****************************************************************
CoorR::CoorR(double* __x, int _len) : Coor(__x, _len) {
  rIsCalculated = 0;
}

CoorR::CoorR() : Coor() {
  rIsCalculated = 0;
}

void CoorR::attach(double* __x, int _len) {
  Coor::attach(__x, _len);
}

void CoorR::calculateR() {
  _r = x[0]*x[0];
  for (int i = 1; i < len; i++)
    _r += (x[i]*x[i]);
  _r = sqrt(_r);
}


// ****************************************************************
// *                          COORSPIN                            *
// ****************************************************************
CoorSpin::CoorSpin(double* __x, int _len, int* __spin) : CoorR(__x, _len) {
  _spin = __spin;
}

void CoorSpin::attach(double* __x, int _len, int* __spin) {
  CoorR::attach(__x, _len);
  _spin = __spin;
}


// ****************************************************************
// *                        COORSPINDIFF                          *
// ****************************************************************
CoorSpinDiff::CoorSpinDiff(double* __x, int _len, int* __spin, double _h) : CoorSpin(__x, _len, __spin) {
  h = _h;
  twoh = 2*h;
  rDiff = new double[len*2];
}

void CoorSpinDiff::attach(double* __x, int _len, int* __spin, double _h) {
  CoorSpin::attach(__x, _len, __spin);
  h = _h;
  twoh = 2*h;
  twoLen = 2*len;
  rDiff = new double[twoLen + 1];
}

void CoorSpinDiff::calculateR() {
  _r = x[0]*x[0];
  for (int i = 1; i < len; i++)
    _r += (x[i]*x[i]);
  _r = sqrt(_r);
  rDiff[twoLen] = _r;
}

void CoorSpinDiff::calculateDiffs() {
  double rTemp;
  for (int i=0; i<len; i++) {

    x[i] += h;
    rTemp = x[0]*x[0];
    for (int j = 1; j < len; j++)
      rTemp += (x[j]*x[j]);
    rDiff[i] = sqrt(rTemp);

    x[i] -= twoh;
    rTemp = x[0]*x[0];
    for (int j = 1; j < len; j++)
      rTemp += (x[j]*x[j]);
    rDiff[i+len] = sqrt(rTemp);

    x[i] += h;
  }
}

double& CoorSpinDiff::operator()(int num) {
  return x[num];
}

double CoorSpinDiff::operator()(int num, int currentNumber) {
  return x[num] + ((currentNumber == num) - (currentNumber == num+len))*h; 
}

void CoorSpinDiff::copy(CoorSpinDiff* copyCoor) {
  _r = copyCoor->r();
  for (int i=0; i <= twoLen; i++)
    rDiff[i] = copyCoor->rdiff(i);
}

