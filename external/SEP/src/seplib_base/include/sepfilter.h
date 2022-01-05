#ifndef SEPFILTER_H
#define SEPFILTER_H yada
#include<prototypes.h>

#ifdef __cplusplus
extern "C" {
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int energy(char *intag,char *outtag,float *data, float *data2,int n1, int n2,
 int n3,int lwind,int j1, int n1new,int norm,int ave);
int mksinc (float *sinc,int lsinc,float d,float *space);
int toep (int m,float *r,float *f,float *g,float *a);
_XFUNCPROTOEND
#else
int energy();
#endif
#ifdef __cplusplus
}
#endif
#endif

