
/*******************   lib.h  *******************/
    /*
     * The definition module                              
     *                      lib.h                    
     * for the library function common for all C programs.
     */

     /* Standard ANSI-C include files */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <malloc.h>
#include <string.h> 
#include <time.h>

#define   NULL_PTR   (void *) 0
#define   ZERO       1.0E-10
#define   INFINITY   1.0E15
#define   UL         unsigned long

         /* a macro used in function pythag() */

static float sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)

      /* Dynamical memory allocation procedure */

#define malloc             "don't call malloc directly!"
#define MALLOC(num_elem, type, func, var)\
               (type *)mem_alloc((num_elem) * sizeof(type), func, var) 
#define calloc             "don't call calloc directly!"
#define CALLOC(num_elem, type, func, var)\
               (type *)mem_calloc((num_elem), sizeof(type), func, var) 
#define realloc             "don't call realloc directly!"
#define REALLOC(ptr, num_elem, type, func, var)\
               (type *)mem_realloc(ptr, (num_elem), sizeof(type), func, var) 

     /* Macro definitions for integer arguments only */

#define   MAX(a,b)   ( ((a) > (b)) ? (a) : (b) )
#define   MIN(a,b)   ( ((a) < (b)) ? (a) : (b) )
#define   PARITY(a)  ( (a) % 2) 
#define   PHASE(a)   (1 - 2 * (abs(a) % 2))
#define   MOD(a,b)   ((a) - ((a)/(b)) * b)
#define   SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))


     /* Function declarations */

void  *mem_alloc(size_t, char *, char *);
void  *mem_calloc(size_t, size_t, char *, char *);
void  *mem_realloc(void *, size_t, size_t, char *, char *);
void  **matrix(int, int, int);
void free_matrix(void **, int);
void rk4(double *, double *, int, double, double, double  *,
	           void (*derivs)(double, double *, double *));
void ludcmp(double **, int, int *, double*);
void lubksb(double **, int, int *, double *);
void tqli(double *, double *, int, double **);
void tred2(double **, int, double *, double *);
double pythag(double, double);
void gauleg(double, double, double *, double *, int);
double simpson ( double, double, int, double (*)(double) );
double trapezoidal_rule ( double, double, int, double (*)(double) );
void spline(double *, double *, int, double, double, double *);
void splint(double *, double *, double *, int, double, double *);
void polint(double *, double *, int, double, double *, double *);
double rtbis(double(*func)(double), double, double, double);
double rtsec(double( *func)(double), double, double, double);
double rtnewt(void ( *funcd)(double, double *, double *), double, double, double);
double zbrent(double( *func)(double), double, double, double);
double ran0(long *);
double ran1(long *);
double ran2(long *);
double ran3(long *);

/****************************  end lib.h  *************/

      
