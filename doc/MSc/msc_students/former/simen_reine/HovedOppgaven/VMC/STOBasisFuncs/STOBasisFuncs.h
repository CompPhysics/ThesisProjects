#ifndef STOBasisFuncs_IS_INCLUDED
#define STOBasisFuncs_IS_INCLUDED

// ****************************************************************
// *                        STOBASISFUNCS                         *
// ****************************************************************
template <class Param>
class STOBasisFuncs {

 protected:
  double constant, expConst, rExp, r;

 public:
  STOBasisFuncs() {}
  virtual void initConsts(int n, int l, double xsi)
    {
      rExp = n-l-1;
      expConst = -xsi;
      double faculty = 1;
      for (int i=1; i<2*n+1; i++) faculty*=i;
      constant = pow(2*xsi, n+0.5)/sqrt(faculty);
    }
  virtual double operator ()(Param& coordinate, int currentNumber)
    { 
      r = coordinate.r(currentNumber);
      return constant*pow(r,rExp)*exp(expConst*r);
    }

};

#endif
