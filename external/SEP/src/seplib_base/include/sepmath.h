#ifndef SEPMATH_H
#define SEPMATH_H yada

/* Include file <sepmath.h>  This file defines all the functions in
 * the sepmath library. It is include in <seplib.h>. It should also be
 * included in all the c-files in the sepmath library.
 *
 * Author 10/23/91 D.Nichols - split off from seplib.h
 */

#include <prototypes.h>
#include <math.h>
#ifndef pi
#ifndef M_PI
static double snftEkd=3.14159265358979323846264338327950288419716939937510;
#define pi snftEkd
#else
#define pi M_PI
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*typedef struct { float re, im;} complex; */
#ifndef OLDcomplex
typedef struct { float re, im;} d0u1m2m3y4cmplx;
#define OLDcomplex d0u1m2m3y4cmplx
#endif

/* !__cplusplus */
/*#endif */

#ifndef MAX
#define MAX(a,b) ( ((a)>(b)) ? (a):(b) )
#endif
#ifndef MIN
#define MIN(a,b) ( ((a)<(b)) ? (a):(b) )
#endif

/* C++ has its own complex arithmetic class */
#ifdef __cplusplus  
#undef SEP_COMPLEX
#endif

#ifdef SEP_COMPLEX
#ifndef CWP_H
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern float cabs2(complex);
extern complex scadd(complex, complex);
extern float sccabs(complex );
extern complex scdiv(complex , complex );
extern complex scexp(complex);
extern complex sciexp(float);
extern complex scinv(complex);
extern complex sclog(complex);
extern complex scmplx(float , float);
extern complex scmult(complex, complex);
extern complex scneg(complex);
extern complex sconjg(complex);
extern complex scsmult(complex, float);
extern complex scspow(complex,float);
extern complex scsqrt(complex);
extern complex scsub(complex,complex);
extern float sqroot (float, float, float);
_XFUNCPROTOEND
#endif
#endif
#endif

#ifdef __cplusplus
}
#endif
#endif
