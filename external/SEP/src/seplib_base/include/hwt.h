#include<prototypes.h>
#ifndef HWT_H
#define HWT_H 12




#define ABSENT     -1999.25
#define PI 3.14159265358979323846264338327950288419716939937510
#define YES 1
#define NO 0
#define NOTIME -1.0
#define sig(x)   x>=0?(x==0?0:1):-1
#ifndef ABS
#define ABS(a) ( ((a)>=(0.0)) ? (a):(-a) )
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double dx,dy,dz;
} v3D;





typedef struct {
  float x,y,z,v,t,l;
} point3D;

typedef struct {
  point3D* point;
} ray3D;

typedef struct {
  float xshot,yshot,zshot; /* source coordinates */
  float aimi,aima;

  int ni; /* number of steps in the aperture direction */
  float oi;     /* origin of the aperture direction = 0.0 */
  float si;     /* aperture step */

  int na;       /* number of steps in the radial direction */
  float oa;     /* origin of the radial direction = 0.0 */  float sa;     /* radial step */

  int nt;
  float ot;
  float st;
  point3D ***c0;
}rays;

typedef struct { /* 1=depth 2=inline 3=xline */
  int   n1,n2,n3;
  float o1,o2,o3;
  float d1,d2,d3;
  float *c0;
}cube;


typedef struct {
  int A_i,A_a; /* ray A */
  int B_i,B_a; /* ray B */
  int C_i,C_a; /* ray C */
} triplet;

typedef struct {
  triplet *R3;
  int nR3;
  ray3D rayA, rayB, rayC;
} triplets;

typedef struct {
  int A_i,A_a,ia; /* index A */
  int B_i,B_a,ib; /* index B */
  int C_i,C_a,ic; /* index C */
  int D_i,D_a,id; /* index D */
} tetrahedron;

typedef struct {
  tetrahedron **T4;
  int nT4;
} tetrahedra;


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
/*hwt_gen functions */
void transfer_rays(rays *pRR, float *xrays, float *yrays, float *zrays, float *vrays);
int init_cube(cube *pVV,float *vel,int *n,float *o,float *d);
void init_rays(rays *pRR,float *shot,int nt,float st,float aimi, float aima,float oi,float si,int ni,float sa,int na);
void allocate_rays(rays* pRR);

/*hwt_ray  */
void tracerays(rays* pRR,cube* pVV);
point3D raytr(rays* pRR,cube* pVV,int it, int ii, int ia);
point3D wfttr(rays* pRR,cube* pVV,int it, int ii, int ia);
int iscusp(rays* pRR,int it, int ii, int ia);
int gottooclose(rays* pRR,cube *pVV,int it,int ii,int ia);
double get_vel_lin(cube* pVV,point3D* p);
double get_vel_pinc(cube* pVV,point3D* p);
void init_wf(rays* pRR,cube* pVV);

/*hwt teselate*/
void tesselate_rays(rays* pRR,cube *pVV,float *times,int fast);

/*hwt_functions*/
double v_l(v3D*);
double v_l2(v3D*);
double v_sina(v3D*,v3D*);
double v_cosa(v3D*,v3D*);
double v_a(v3D*,v3D*);
double v_dp(v3D*,v3D*);
v3D    v_vp(v3D*,v3D*);
double v_mp(v3D*,v3D*,v3D*);
double volume(v3D*,v3D*,v3D*);
float det3(float*);
float det2(float*);
int signum(double);
void set_cube(cube,float);
void tesselate_close(triplets *ptrp, tetrahedra *ptet);
#else
void transfer_rays();
int init_cube();
void init_rays();
void allocate_rays();
void tracerays();
point3D raytr();
point3D wfttr();
int iscusp();
int gottooclose();
double get_vel_lin();
double get_vel_pinc();
void init_wf();
double v_l();
double v_l2();
double v_sina();
double v_cosa();
double v_a();
double v_dp();
v3D    v_vp();
double v_mp();
double volume();
float det3();
float det2();
int signum();
void measure_rays();
void set_cube();
void tesselate_rays();
#endif







#ifdef __cplusplus
}
#endif





#endif
