#ifndef SEPVELAN_H
#define SEPVELAN_H
#include<sepvelan.h> 
#define NMO_NEIGHBOR(adj, add, slow,x t0, dt, n, zz, tt) CCALLSFSUB9(NMO_NEIGHBOR,nmo_neighbor,INT,INT,FLOATV,FLOAT,FLOAT, FLOAT, INT, FLOATV, FLOATV, adj, add, slow,    x, t0, dt, n,zz,  tt )
#define NMO_NEIGHBOR_W(adj, add, slow,x t0, dt, n, zz, tt) CCALLSFSUB9(NMO_NEIGHBOR_W,nmo_neighbor_w,INT,INT,FLOATV,FLOAT,FLOAT, FLOAT, INT, FLOATV, FLOATV, adj, add, slow,    x, t0, dt, n,zz,  tt )
#define RMS2INT(inverse, dt, vrms, nt,  vint ) CCALLSFSUB5(RMS2INT, rms2int, INT, FLOAT, FLOATV, INT, FLOATV, inverse, dt, vrms, nt,  vint )
#define VELSIMP(adj,add, t0,dt,x0,dx,s0,ds, nt,nx,ns, modl, data) CCALLSFSUB13(VELSIMP, velsimp, INT, INT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, INT, INT, INT, FLOATVV, FLOATVV, adj,add, t0,dt,x0,dx,s0,ds, nt,nx,ns, modl, data)
#define VELTRAN(adj,add,psun,s02,anti,t0,dt,x0,dx,s0,ds,nt,nx,ns,model,data) CCALLSFSUB16(VELTRAN, veltran, INT, INT, INT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, INT, INT, INT, FLOATVV, FLOATVV, adj,add,psun,s02,anti,t0,dt,x0,dx,s0,ds,nt,nx,ns,model,data)
#define VINT2RMS(inverse, vminallow, dt, vint, nt,  vrms)  CCALLSFSUB6(VINT2RMS,
 vint2rms, INT, FLOAT, FLOAT, FLOATV, INT, FLOATV, inverse, vminallow, dt, vint,
 nt,  vrms)
#endif
