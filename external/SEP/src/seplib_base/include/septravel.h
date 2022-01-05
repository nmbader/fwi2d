#ifndef SEPTRAVEL_H
#define SEPTRAVEL_H yada
#include<prototypes.h>

#ifdef __cplusplus
extern "C" {
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
 int fastmarch (int order,float s1, float s2, float s3,
		int b1, int b2, int b3, 
		int nz, int ny, int nx,
		float dz, float dy, float dx, 
		float slow0, float *slow, float *ttime);
 void update (int order, int i1, int i2, int i3, float* tj,unsigned char* mj, float s);
	int hwt_trace_rays(int *n,float *o,float *d,float *shot,float aimi,float aima,
 float oi,float si,int ni, float sa, int na,float st,int nt,float *vel,
 float *xray,float *yray,float *zray,float *vray);
	int hwt_travel_cube(int *n,float *o,float *d,float *shot,float aimi,float aima,
 float oi,float si,int ni, float sa, int na,float st,int nt,float *vel,
	int fast, float *times);
void cvupdate (double cvtime,float* tj,unsigned char* mj);
_XFUNCPROTOEND
#else
int fastmarch();
void update();
void cvupdate();
int hwt_trace_rays();
int hwt_travel_cube();
#endif
#ifdef __cplusplus
}
#endif
#endif

