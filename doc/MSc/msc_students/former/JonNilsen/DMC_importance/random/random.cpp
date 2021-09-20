#include "random.hpp"


// normal random variate generator using box-muller
double Random::gran(double s, double m)
{        /* mean m, standard deviation s */
  double x1, x2, w, y1;
  static double y2;
  static int use_last = 0;

  if (use_last)        /* use value from previous call */
    {
      y1 = y2;
      use_last = 0;
    }
  else
    {
      do {
	x1 = 2.0 * ran1() - 1.0;
	x2 = 2.0 * ran1() - 1.0;
	w = x1 * x1 + x2 * x2;
      } while ( w >= 1.0 );
      w = sqrt( (-2.0 * log( w ) ) / w );
      y1 = x1 * w;
      y2 = x2 * w;
      use_last = 1;
    }
  return( m + y1 * s );
}


/*
** The function
**           ran1()
** is an "Minimal" random number generator of Park and Miller
** (see Numerical recipe page 280) with Bays-Durham shuffle and
** added safeguards. Call with idum a negative integer to initialize;
** thereafter, do not alter idum between sucessive deviates in a
** sequence. RNMX should approximate the largest floating point value
** that is less than 1.
** The function returns a uniform deviate between 0.0 and 1.0
** (exclusive of end-point values).
*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double Random::ran1()
{
   int             j;
   long            k;
   static long     iy=0;
   static long     iv[NTAB];
   double          temp;

   if (idum <= 0 || !iy) {
      if (-(idum) < 1) idum=1;
      else idum = -(idum);
      for(j = NTAB + 7; j >= 0; j--) {
         k     = (idum)/IQ;
         idum = IA*(idum - k*IQ) - IR*k;
         if(idum < 0) idum += IM;
         if(j < NTAB) iv[j] = idum;
      }
      iy = iv[0];
   }
   k     = (idum)/IQ;
   idum = IA*(idum - k*IQ) - IR*k;
   if(idum < 0) idum += IM;
   j     = iy/NDIV;
   iy    = iv[j];
   iv[j] = idum;
   if((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// End: function ran1()

// Long period >2x10^18 random number generator of L'Ecuyer with 
// Bays-Durham shuffle. Return a uniform deviate in the interval (0,1). 
// Call with int variable as argument set to a negative value, and
// don't this variable variable unless you want to reseed the generator

double Random::ran2 () {
    const int IM1 = 2147483563, IM2 = 2147483399;
    const double AM=(1.0/IM1);
    const int IMM1 = IM1-1;
    const int IA1 = 40014, IA2 = 40692, IQ1 = 53668, IQ2 = 52774;
    const int IR1 = 12211, IR2 = 3791, NTAB = 32;
    const int NDIV = 1+IMM1/NTAB;
    const double EPS = 3.0e-16, RNMX = 1.0-EPS;

    int j, k;
    static int idum2=123456789, iy = 0;
    static int iv[NTAB];
    double temp;

    if (idum <= 0) {
        idum = (idum == 0 ? 1 : -idum);
        idum2=idum;
        for (j=NTAB+7;j>=0;j--) {
            k=idum/IQ1;
            idum=IA1*(idum-k*IQ1)-k*IR1;
            if (idum < 0) idum += IM1;
            if (j < NTAB) iv[j] = idum;
        }
        iy=iv[0];
    }
    k=idum/IQ1;
    idum=IA1*(idum-k*IQ1)-k*IR1;
    if (idum < 0) idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

