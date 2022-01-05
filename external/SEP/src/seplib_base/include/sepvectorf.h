#include <sepvector.h>

/*C CALLS TO FORTRAN */
#define APROD(mode,m,n,x,y,leniw,iw) CCALLSFSUB7(APROD, aprod, INT, INT, INT, FLOATV, FLOATV, INT, INTV, mode,m,n,x,y,leniw,iw)
#define ADJNULL(adj,add,x,nx,y,ny) CCALLSFSUB6(ADJNULL, adjnull, INT, INT, FLOATV, INT, FLOATV, INT, adj,add,x,nx,y,ny)
#define SCALE(factor,n,data) CCALLSFSUB3(SCALE, scale, FLOAT, INT, FLOATV, factor,n,data)
#define COPY(n,xx,yy) CCALLSFSUB3(COPY, copy, INT, FLOATV, FLOATV, n,xx,yy)
#define INV(input,output,n) CCALLSFSUB3(INV, inv, FLOATV, FLOATV, INT, input,output,n)
#define SQRT(input,output,n) CCALLSFSUB3(SQRT, sqrt, FLOATV, FLOATV, INT, input,output,n)
#define IPOW(input,factor,output,n) CCALLSFSUB4(IPOW, ipow, FLOATV, INT,FLOATV, INT, input,factor,output,n)
