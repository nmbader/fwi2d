#ifndef SEPFFT_H
#define SEPFFT_H yada
/*this should be changed - don't need to pull in all of the includes that
this implies */
#include<prototypes.h> 
#include <math.h>
#include <complex.h>
#ifndef pi
static double snftEkd=3.14159265358979323846264338327950288419716939937510;
#define pi snftEkd
#endif

#ifndef __cplusplus  /* C++ has its own complex arithmetic class */
#ifndef OLDcomplex
typedef struct { float re, im;} d0u1m2m3y4cmplx;
#define OLDcomplex d0u1m2m3y4cmplx
#endif
#endif





#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern int cefft(float complex*,int ,int ,float);
extern int cvfft(float complex*,int ,int ,int ,double );
extern int refft(register float complex*,int,int,float,int);
extern int rvfft(register float complex*,int,int,int,float,int);
_XFUNCPROTOEND
#else
extern int cefft ();
extern int cvfft ();
extern int refft();
extern int rvfft();
#endif


#endif /*end of  sepfft header defined */
