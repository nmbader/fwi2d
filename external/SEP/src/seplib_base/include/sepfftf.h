#include<sepfft.h>
#include<cfortran.h>

/*C CALLS TO FORTRAN */
#ifdef GNU
#define FT1AXIS(adj, sign1, n1,n2, cx) CCALLSFSUB5(FT1AXIS_,ft1axis_,INT,FLOAT,INT,INT,FLOATVV,adj, sign1, n1,n2, cx)
#define FT2AXIS(adj, sign1, n1,n2, cx) CCALLSFSUB5(FT2AXIS_,ft2axis_,INT,FLOAT,INT,INT,FLOATVV,adj, sign1, n1,n2, cx)
#define FTH( adj,sign, m1, n12, cx) CCALLSFSUB5(FTH, fth, INT, FLOAT, INT, INT, FLOATVV,adj,sign,m1,n12,cx)
#define FTU(signi, nx, cx ) CCALLSFSUB3(FTU_, ftu_, FLOAT,INT,FLOATV,signi,nx,cx)
#define CFFT1(a, npts, iop ) CCALLSFSUB3(CFFT1_, cfft1_, FLOATV,INT,INT,a,npts,iop)
#define FT3D(n1, n2, n3, cx, sign1, sign2, sign3, center1, center2, center3)  CCALLSFSUB10(FT3D_SUB_,ft3d_sub_,INT,INT,INT,FLOATV,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,n1, n2, n3, cx, sign1, sign2, sign3, center1, center2, center3)
#define ROWCC( N1, N2, CX, SIGN2, SCALE ) \
  CCALLSFSUB5(ROWCC,rowcc,INT,INT,FLOATV,FLOAT,FLOAT, \
        N1, N2, CX, SIGN2, SCALE )
#else
#define FT1AXIS(adj, sign1, n1,n2, cx) CCALLSFSUB5(FT1AXIS,ft1axis,INT,FLOAT,INT,INT,FLOATVV,adj, sign1, n1,n2, cx)
#define FT2AXIS(adj, sign1, n1,n2, cx) CCALLSFSUB5(FT2AXIS,ft2axis,INT,FLOAT,INT,INT,FLOATVV,adj, sign1, n1,n2, cx)
#define FTH( adj,sign, m1, n12, cx) CCALLSFSUB5(FTH, fth, INT, FLOAT, INT, INT, FLOATVV,adj,sign,m1,n12,cx)
#define FTU(signi, nx, cx ) CCALLSFSUB3(FTU, ftu, FLOAT,INT,FLOATV,signi,nx,cx)
#define CFFT1(a, npts, iop ) CCALLSFSUB3(CFFT1, cfft1, FLOATV,INT,INT,a,npts,iop)
#define FT3D(n1, n2, n3, cx, sign1, sign2, sign3, center1, center2, center3)  CCALLSFSUB10(FT3D_SUB,ft3d_sub,INT,INT,INT,FLOATV,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,n1, n2, n3, cx, sign1, sign2, sign3, center1, center2, center3)
#define ROWCC( N1, N2, CX, SIGN2, SCALE ) \
  CCALLSFSUB5(ROWCC,rowcc,INT,INT,FLOATV,FLOAT,FLOAT, \
        N1, N2, CX, SIGN2, SCALE )
#endif
