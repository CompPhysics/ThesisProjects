#include "Random.h"

// ****************************************************************
// *                          RANDOM2                             *
// ****************************************************************
void Random2::demo(int num, FILE* stream) {
  for (int i = 0; i < num; i++)
    fprintf(stream, "%f\n", getNum());
}

// ****************************************************************
// *                            RAN0                              *
// ****************************************************************
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
double Ran0::getNum(void) {
   long     k;
   double ans;

   dummy ^= MASK;
   k      = (dummy)/IQ;
   dummy  = IA*(dummy - k*IQ) - IR*k;
   if(dummy < 0) dummy += IM;
   ans    = AM*(dummy);
   dummy ^= MASK;
   return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

// ****************************************************************
// *                            RAN1                              *
// ****************************************************************
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double Ran1::getNum(void)
{
   int             j;
   long            k;
   static long     iy=0;
   static long     iv[NTAB];
   double          temp;

   if (dummy <= 0 || !iy) {
      if (-(dummy) < 1) dummy=1;
      else dummy = -(dummy);
      for(j = NTAB + 7; j >= 0; j--) {
         k     = (dummy)/IQ;
         dummy = IA*(dummy - k*IQ) - IR*k;
         if(dummy < 0) dummy += IM;
         if(j < NTAB) iv[j] = dummy;
      }
      iy = iv[0];
   }
   k     = (dummy)/IQ;
   dummy = IA*(dummy - k*IQ) - IR*k;
   if(dummy < 0) dummy += IM;
   j     = iy/NDIV;
   iy    = iv[j];
   iv[j] = dummy;
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
