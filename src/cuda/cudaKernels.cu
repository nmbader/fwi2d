#include <algorithm>
#include <stdio.h>
#include "cudaMisc.h"
#include "cudaCppWrappers.h"
#include "injector.hpp"
#include "we_op.hpp"

#define BLOCK_SIZE 16

__constant__ static data_t bnd_coef1[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
__constant__ static data_t coef[10] = {-3.0/4, -5.0/6, -1.0/24, 1.0/6, 1.0/2, 1.0/2, 1.0/6, -1.0/8, 1.0/6, -1.0/8};
__constant__ static data_t bnd_coef2[384] = {920.0/289,-59.0/68,-7549318.0/34157643,-17440994.0/184112825,0,0,0,0,-1740.0/289,0,295314153.0/366719282,262545878.0/1218454145,0,0,0,0,1128.0/289,59.0/68,-19250923.0/26254840,-12192537.0/324892213,0,0,0,0,-308.0/289,0,42283069.0/254173229,-43013531.0/427521546,0,0,0,0,0,0,-18700293.0/1355757959,18700293.0/1355757959,0,0,0,0,0,0,-3.0/833,3.0/833,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    12.0/17,0,89562243.0/385991318,16914691.0/272440171,0,0,0,0,-59.0/68,0,-47979680.0/48852831,-18034913.0/120051851,0,0,0,0,2.0/17,0,299262497.0/373256703,16156647.0/200473586,0,0,0,0,3.0/68,0,-14723969.0/177744748,22633571.0/584543661,0,0,0,0,0,0,46802031.0/1628311862,-46802031.0/1628311862,0,0,0,0,0,0,441.0/181507,-441.0/181507,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    -96.0/731,59.0/172,-47632751.0/164317197,-5570723.0/375470930,0,0,0,0,118.0/731,0,54282598.0/49343777,23802793.0/215253532,0,0,0,0,-16.0/731,-59.0/172,-39119273.0/25083370,-35971870.0/61324629,-26254.0/557679,0,0,0,-6.0/731,0,360454121.0/368940022,17254963.0/80047776,1500708.0/7993399,0,0,0,0,0,-18024731.0/79673021,24178273.0/88099647,-26254.0/185893,0,0,0,0,0,-870707.0/620833782,960119.0/1147305747,13564.0/23980197,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    -36.0/833,0,54741803.0/948483020,-13602043.0/389676498,0,0,0,0,177.0/3332,0,-35820026.0/359121865,24921773.0/534548210,0,0,0,0,-6.0/833,0,813284372.0/948584067,30057666.0/158897885,1500708.0/9108757,0,0,0,-9.0/3332,0,-95056924.0/90903639,-23417695.0/47008088,-7476412.0/9108757,-2.0/49,0,0,0,0,23159719.0/99948527,110687545.0/265515812,4502124.0/9108757,8.0/49,0,0,0,0,-3671038.0/2687426923,-1063649.0/8893843,1473580.0/9108757,-6.0/49,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,-9437957.0/1931986386,9437957.0/1931986386,0,0,0,0,0,0,17289851.0/489388053,-17289851.0/489388053,0,0,0,0,0,0,-66355919.0/327412264,24178273.0/98343792,-564461.0/4461432,0,0,0,0,0,17638343.0/74566894,103749401.0/243793650,375177.0/743572,1.0/6,0,0,0,0,-19321801.0/295845927,-50677283.0/62943042,-280535.0/371786,-5.0/6,-1.0/24,0,0,0,5130905.0/5183662662,35039615.0/213452232,1118749.0/2230716,1.0/2,1.0/6,0,0,0,0,0,-1.0/8,1.0/6,-1.0/8,0,0,0,0,0,0,0,0,0,
    0,0,-1.0/784,1.0/784,0,0,0,0,0,0,8673.0/2904112,-8673.0/2904112,0,0,0,0,0,0,-403062.0/320810033,960119.0/1280713392,3391.0/6692148,0,0,0,0,0,-1920494.0/1377228165,-1063649.0/8712336,368395.0/2230716,-1.0/8,0,0,0,0,5130905.0/5183662662,35039615.0/213452232,1118749.0/2230716,1.0/2,1.0/6,0,0,0,-117380.0/2351569839,-3290636.0/80044587,-5580181.0/6692148,-3.0/4,-5.0/6,-1.0/24,0,0,0,0,1.0/6,1.0/2,1.0/2,1.0/6,0,0,0,0,0,-1.0/8,1.0/6,-1.0/8
};

__global__ void cudaDz_interior(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d){
	
    int nc1=4;

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE + nc1;
    int blockx = BLOCK_SIZE;
    int blockz = BLOCK_SIZE+4;  // shared block size with halo
    int ixL = threadIdx.x;
    int izL = threadIdx.y;
    int i1=ix*nz;

    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    data_t  scale=1.0/d;

    if (blockIdx.x == gridDim.x-1) blockx = nx-BLOCK_SIZE*(gridDim.x-1); // rightmost block
    if (blockIdx.y == gridDim.y-1) blockz = nz-2*nc1-BLOCK_SIZE*(gridDim.y-1)+4; // bottommost block

    __shared__ data_t sh_in[BLOCK_SIZE][BLOCK_SIZE+4]; // allocate shared memory: BLOCKSIZE + halos

    if (ixL<blockx && izL<blockz-4) sh_in[ixL][izL+2] = in[i1+iz]; // copy the stencil tile
    if (ixL<blockx && izL<2) {
        sh_in[ixL][izL] = in[i1+iz-2]; // top halo
        sh_in[ixL][blockz-1-izL] = in[i1+iz+blockz-2-1-2*izL]; // bottom halo
    }
    __syncthreads();

    if (ix < nx && iz < nz - nc1)
        out[ i1+ iz] = add*out[ i1+ iz] + scale * (coef0 * (sh_in[ixL][izL+3] - sh_in[ixL][izL+1]) + coef1 * (sh_in[ixL][izL+4] - sh_in[ixL][izL]));
}

__global__ void cudaDz_bnd(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d){
	
    int nc2=6;

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;
    int iz = threadIdx.y;
    int blockx = BLOCK_SIZE;
    int blockz = 6;  // shared block size with halo
    int ixL = threadIdx.x;
    int i1=ix*nz;

    if (blockIdx.x == gridDim.x-1) blockx = nx-BLOCK_SIZE*(gridDim.x-1); // rightmost block

    __shared__ data_t sh_in[BLOCK_SIZE][6]; // allocate shared memory: BLOCKSIZE + halos

    if (blockIdx.y == 0) // top boundary
    {
        if (ixL<blockx) {
            sh_in[ixL][iz] = in[i1+iz]; // copy the stencil tile
            if (iz<2) sh_in[ixL][blockz-1-iz] = in[i1+blockz-1-iz]; // bottom halo
        }
        __syncthreads();

        if (ix<nx) out[i1+iz] = add* out[i1+iz] + (bnd_coef1[iz*nc2] * sh_in[ixL][0] + bnd_coef1[iz*nc2+1] * sh_in[ixL][1] + bnd_coef1[iz*nc2+2] * sh_in[ixL][2] + bnd_coef1[iz*nc2+3] * sh_in[ixL][3] + bnd_coef1[iz*nc2+4] * sh_in[ixL][4] + bnd_coef1[iz*nc2+5] * sh_in[ixL][5]) / d;
    }
    else  // bottom boundary
    {
        if (ixL<blockx) {
            sh_in[ixL][iz] = in[i1+nz-1-iz]; // copy the stencil tile
            if (iz<2) sh_in[ixL][blockz-1-iz] = in[i1+nz-blockz+iz]; // top halo
        }
        __syncthreads();

        if (ix<nx) out[i1+nz-1-iz] = add* out[i1+nz-1-iz] + (-bnd_coef1[iz*nc2] * sh_in[ixL][0] - bnd_coef1[iz*nc2+1] * sh_in[ixL][1] - bnd_coef1[iz*nc2+2] * sh_in[ixL][2] - bnd_coef1[iz*nc2+3] * sh_in[ixL][3] - bnd_coef1[iz*nc2+4] * sh_in[ixL][4] - bnd_coef1[iz*nc2+5] * sh_in[ixL][5]) / d;
    }
}

__global__ void cudaDx_interior(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d){
	
    int nc1=4;

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE + nc1;
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;
    int blockx = BLOCK_SIZE + 4; // shared block size with halo
    int blockz = BLOCK_SIZE;  
    int ixL = threadIdx.x;
    int izL = threadIdx.y;
    int i1=ix*nz;

    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    data_t  scale=1.0/d;

    if (blockIdx.x == gridDim.x-1) blockx = nx-2*nc1-BLOCK_SIZE*(gridDim.x-1)+4; // rightmost block
    if (blockIdx.y == gridDim.y-1) blockz = nz-BLOCK_SIZE*(gridDim.y-1); // bottommost block

    __shared__ data_t sh_in[BLOCK_SIZE+4][BLOCK_SIZE]; // allocate shared memory: BLOCKSIZE + halos

    if (ixL<blockx-4 && izL<blockz) sh_in[ixL+2][izL] = in[i1+iz]; // copy the stencil tile

    if (izL<blockz && ixL<2) {
        sh_in[ixL][izL] = in[(ix-2)*nz+iz]; // left halo
        sh_in[blockx-1-ixL][izL] = in[(ix+blockx-2-1-2*ixL)*nz+iz]; // right halo
    }
    __syncthreads();

    if (ix < nx-nc1 && iz < nz)
        out[ i1+ iz] = add*out[ i1+ iz] + scale * (coef0 * (sh_in[ixL+3][izL] - sh_in[ixL+1][izL]) + coef1 * (sh_in[ixL+4][izL] - sh_in[ixL][izL]));
}

__global__ void cudaDx_bnd(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d){
	
    int nc2=6;

    int ix = threadIdx.x;
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;
    int blockx = 6; // shared block size with halo
    int blockz = BLOCK_SIZE;  
    int izL = threadIdx.y;
    int i1=ix*nz;

    if (blockIdx.y == gridDim.y-1) blockz = nz-BLOCK_SIZE*(gridDim.y-1); // bottommost block

    __shared__ data_t sh_in[6][BLOCK_SIZE]; // allocate shared memory: BLOCKSIZE + halos

    if (blockIdx.x == 0) // left boundary
    {
        if (izL<blockz) {
            sh_in[ix][izL] = in[i1+iz]; // copy the stencil tile
            if (ix<2) sh_in[blockx-1-ix][izL] = in[(blockx-1-ix)*nz+iz]; // right halo
        }
        __syncthreads();

        if (iz<nz) out[i1+iz] = add* out[i1+iz] + (bnd_coef1[ix*nc2] * sh_in[0][izL] + bnd_coef1[ix*nc2+1] * sh_in[1][izL] + bnd_coef1[ix*nc2+2] * sh_in[2][izL] + bnd_coef1[ix*nc2+3] * sh_in[3][izL] + bnd_coef1[ix*nc2+4] * sh_in[4][izL] + bnd_coef1[ix*nc2+5] * sh_in[5][izL]) / d; 
    }
    else  // right boundary
    {
        if (izL<blockz) {
            sh_in[ix][izL] = in[(nx-1-ix)*nz+iz]; // copy the stencil tile
            if (ix<2) sh_in[blockx-1-ix][izL] = in[(nx-blockx+ix)*nz+iz]; // left halo
        }
        __syncthreads();

        i1=(nx-1-ix)*nz;
        if (iz<nz) out[i1+iz] = add* out[i1+iz] + (-bnd_coef1[ix*nc2] * sh_in[0][izL] - bnd_coef1[ix*nc2+1] * sh_in[1][izL] - bnd_coef1[ix*nc2+2] * sh_in[2][izL] - bnd_coef1[ix*nc2+3] * sh_in[3][izL] - bnd_coef1[ix*nc2+4] * sh_in[4][izL] - bnd_coef1[ix*nc2+5] * sh_in[5][izL]) / d;
    }
}

__global__ void cudaMultDz_interior(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t* par, data_t a){
	
    int nc1=4;

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE + nc1;
    int blockx = BLOCK_SIZE;
    int blockz = BLOCK_SIZE+4;  // shared block size with halo
    int ixL = threadIdx.x;
    int izL = threadIdx.y;
    int i1=ix*nz;

    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    data_t  scale=a/d;

    if (blockIdx.x == gridDim.x-1) blockx = nx-BLOCK_SIZE*(gridDim.x-1); // rightmost block
    if (blockIdx.y == gridDim.y-1) blockz = nz-2*nc1-BLOCK_SIZE*(gridDim.y-1)+4; // bottommost block

    __shared__ data_t sh_in[BLOCK_SIZE][BLOCK_SIZE+4]; // allocate shared memory: BLOCKSIZE + halos
    __shared__ data_t sh_par[BLOCK_SIZE][BLOCK_SIZE+4];
    if (ixL<blockx && izL<blockz-4) {
        sh_in[ixL][izL+2] = in[i1+iz]; // copy the stencil tile
        sh_par[ixL][izL+2] = par[i1+iz];
    }
    if (ixL<blockx && izL<2) {
        sh_in[ixL][izL] = in[i1+iz-2]; // top halo
        sh_par[ixL][izL] = par[i1+iz-2]; 
        sh_in[ixL][blockz-1-izL] = in[i1+iz+blockz-2-1-2*izL]; // bottom halo
        sh_par[ixL][blockz-1-izL] = par[i1+iz+blockz-2-1-2*izL];
    }
    __syncthreads();

    if (ix < nx && iz < nz - nc1)
        out[ i1+ iz] = add*out[ i1+ iz] + scale * (coef0 * (sh_in[ixL][izL+3]*sh_par[ixL][izL+3] - sh_in[ixL][izL+1]*sh_par[ixL][izL+1]) + coef1 * (sh_in[ixL][izL+4]*sh_par[ixL][izL+4] - sh_in[ixL][izL]*sh_par[ixL][izL]));
}

__global__ void cudaMultDz_bnd(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t* par, data_t a){
	
    int nc2=6;

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;
    int iz = threadIdx.y;
    int blockx = BLOCK_SIZE;
    int blockz = 6;  // shared block size with halo
    int ixL = threadIdx.x;
    int i1=ix*nz;

    if (blockIdx.x == gridDim.x-1) blockx = nx-BLOCK_SIZE*(gridDim.x-1); // rightmost block

    __shared__ data_t sh_in[BLOCK_SIZE][6]; // allocate shared memory: BLOCKSIZE + halos
    __shared__ data_t sh_par[BLOCK_SIZE][6];

    if (blockIdx.y == 0) // top boundary
    {
        if (ixL<blockx) {
            sh_in[ixL][iz] = in[i1+iz]; // copy the stencil tile
            sh_par[ixL][iz] = par[i1+iz];
            if (iz<2) {
                sh_in[ixL][blockz-1-iz] = in[i1+blockz-1-iz]; // bottom halo
                sh_par[ixL][blockz-1-iz] = par[i1+blockz-1-iz]; 
            }
        }
        __syncthreads();

        if (ix<nx) out[i1+iz] = add* out[i1+iz] + a/d * (bnd_coef1[iz*nc2] * sh_in[ixL][0]*sh_par[ixL][0] + bnd_coef1[iz*nc2+1] * sh_in[ixL][1]*sh_par[ixL][1] + bnd_coef1[iz*nc2+2] * sh_in[ixL][2]*sh_par[ixL][2] + bnd_coef1[iz*nc2+3] * sh_in[ixL][3]*sh_par[ixL][3] + bnd_coef1[iz*nc2+4] * sh_in[ixL][4]*sh_par[ixL][4] + bnd_coef1[iz*nc2+5] * sh_in[ixL][5]*sh_par[ixL][5]);
    }
    else  // bottom boundary
    {
        if (ixL<blockx) {
            sh_in[ixL][iz] = in[i1+nz-1-iz]; // copy the stencil tile
            sh_par[ixL][iz] = par[i1+nz-1-iz];
            if (iz<2) {
                sh_in[ixL][blockz-1-iz] = in[i1+nz-blockz+iz]; // top halo
                sh_par[ixL][blockz-1-iz] = par[i1+nz-blockz+iz];
            }
        }
        __syncthreads();

        if (ix<nx) out[i1+nz-1-iz] = add* out[i1+nz-1-iz] + a/d * (-bnd_coef1[iz*nc2] * sh_in[ixL][0]*sh_par[ixL][0] - bnd_coef1[iz*nc2+1] * sh_in[ixL][1]*sh_par[ixL][1] - bnd_coef1[iz*nc2+2] * sh_in[ixL][2]*sh_par[ixL][2] - bnd_coef1[iz*nc2+3] * sh_in[ixL][3]*sh_par[ixL][3] - bnd_coef1[iz*nc2+4] * sh_in[ixL][4]*sh_par[ixL][4] - bnd_coef1[iz*nc2+5] * sh_in[ixL][5]*sh_par[ixL][5]);
       
    }
}

__global__ void cudaMultDx_interior(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t* par, data_t a){
	
    int nc1=4;

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE + nc1;
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;
    int blockx = BLOCK_SIZE + 4; // shared block size with halo
    int blockz = BLOCK_SIZE;  
    int ixL = threadIdx.x;
    int izL = threadIdx.y;
    int i1=ix*nz;

    data_t coef0=2.0/3;
    data_t coef1=-1.0/12;
    data_t  scale=a/d;

    if (blockIdx.x == gridDim.x-1) blockx = nx-2*nc1-BLOCK_SIZE*(gridDim.x-1)+4; // rightmost block
    if (blockIdx.y == gridDim.y-1) blockz = nz-BLOCK_SIZE*(gridDim.y-1); // bottommost block

    __shared__ data_t sh_in[BLOCK_SIZE+4][BLOCK_SIZE]; // allocate shared memory: BLOCKSIZE + halos
    __shared__ data_t sh_par[BLOCK_SIZE+4][BLOCK_SIZE];
    if (ixL<blockx-4 && izL<blockz) {
        sh_in[ixL+2][izL] = in[i1+iz]; // copy the stencil tile
        sh_par[ixL+2][izL] = par[i1+iz];
    }

    if (izL<blockz && ixL<2) {
        sh_in[ixL][izL] = in[(ix-2)*nz+iz]; // left halo
        sh_par[ixL][izL] = par[(ix-2)*nz+iz];
        sh_in[blockx-1-ixL][izL] = in[(ix+blockx-2-1-2*ixL)*nz+iz]; // right halo
        sh_par[blockx-1-ixL][izL] = par[(ix+blockx-2-1-2*ixL)*nz+iz];
    }
    __syncthreads();

    if (ix < nx-nc1 && iz < nz)
        out[ i1+ iz] = add*out[ i1+ iz] + scale * (coef0 * (sh_in[ixL+3][izL]*sh_par[ixL+3][izL] - sh_in[ixL+1][izL]*sh_par[ixL+1][izL]) + coef1 * (sh_in[ixL+4][izL]*sh_par[ixL+4][izL] - sh_in[ixL][izL]*sh_par[ixL][izL]));
}

__global__ void cudaMultDx_bnd(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t* par, data_t a){
	
    int nc2=6;

    int ix = threadIdx.x;
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;
    int blockx = 6; // shared block size with halo
    int blockz = BLOCK_SIZE;  
    int izL = threadIdx.y;
    int i1=ix*nz;

    if (blockIdx.y == gridDim.y-1) blockz = nz-BLOCK_SIZE*(gridDim.y-1); // bottommost block

    __shared__ data_t sh_in[6][BLOCK_SIZE]; // allocate shared memory: BLOCKSIZE + halos
    __shared__ data_t sh_par[6][BLOCK_SIZE];

    if (blockIdx.x == 0) // left boundary
    {
        if (izL<blockz) {
            sh_in[ix][izL] = in[i1+iz]; // copy the stencil tile
            sh_par[ix][izL] = par[i1+iz]; 
            if (ix<2){
                sh_in[blockx-1-ix][izL] = in[(blockx-1-ix)*nz+iz]; // right halo
                sh_par[blockx-1-ix][izL] = par[(blockx-1-ix)*nz+iz];
            }
        }
        __syncthreads();

        if (iz<nz) out[i1+iz] = add* out[i1+iz] + a/d * (bnd_coef1[ix*nc2] * sh_in[0][izL]*sh_par[0][izL] + bnd_coef1[ix*nc2+1] * sh_in[1][izL]*sh_par[1][izL] + bnd_coef1[ix*nc2+2] * sh_in[2][izL]*sh_par[2][izL] + bnd_coef1[ix*nc2+3] * sh_in[3][izL]*sh_par[3][izL] + bnd_coef1[ix*nc2+4] * sh_in[4][izL]*sh_par[4][izL] + bnd_coef1[ix*nc2+5] * sh_in[5][izL]*sh_par[5][izL]); 
    }
    else  // right boundary
    {
        if (izL<blockz) {
            sh_in[ix][izL] = in[(nx-1-ix)*nz+iz]; // copy the stencil tile
            sh_par[ix][izL] = par[(nx-1-ix)*nz+iz];
            if (ix<2) {
                sh_in[blockx-1-ix][izL] = in[(nx-blockx+ix)*nz+iz]; // left halo
                sh_par[blockx-1-ix][izL] = par[(nx-blockx+ix)*nz+iz];
            }
        }
        __syncthreads();

        i1=(nx-1-ix)*nz;
        if (iz<nz) out[i1+iz] = add* out[i1+iz] + a/d * (-bnd_coef1[ix*nc2] * sh_in[0][izL]*sh_par[0][izL] - bnd_coef1[ix*nc2+1] * sh_in[1][izL]*sh_par[1][izL] - bnd_coef1[ix*nc2+2] * sh_in[2][izL]*sh_par[2][izL] - bnd_coef1[ix*nc2+3] * sh_in[3][izL]*sh_par[3][izL] - bnd_coef1[ix*nc2+4] * sh_in[4][izL]*sh_par[4][izL] - bnd_coef1[ix*nc2+5] * sh_in[5][izL]*sh_par[5][izL]);
    }
}

__global__ void cudaDzz_var_interior(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t * par, data_t a){
	
    int nc1=6;
    data_t d2 = d*d;

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE + nc1;
    int blockx = BLOCK_SIZE;
    int blockz = BLOCK_SIZE+4;  // shared block size with halo
    int ixL = threadIdx.x;
    int izL = threadIdx.y;
    int i1=ix*nz;

    if (blockIdx.x == gridDim.x-1) blockx = nx-BLOCK_SIZE*(gridDim.x-1); // rightmost block
    if (blockIdx.y == gridDim.y-1) blockz = nz-2*nc1-BLOCK_SIZE*(gridDim.y-1)+4; // bottommost block

    __shared__ data_t sh_in[BLOCK_SIZE][BLOCK_SIZE+4]; // allocate shared memory: BLOCKSIZE + halos
    __shared__ data_t sh_par[BLOCK_SIZE][BLOCK_SIZE+4];
    if (ixL<blockx && izL<blockz-4) {
        sh_in[ixL][izL+2] = in[i1+iz]; // copy the stencil tile
        sh_par[ixL][izL+2] = par[i1+iz];
    }

    if (ixL<blockx && izL<2) {
        sh_in[ixL][izL] = in[i1+iz-2]; // top halo
        sh_par[ixL][izL] = par[i1+iz-2];
        sh_in[ixL][blockz-1-izL] = in[i1+iz+blockz-2-1-2*izL]; // bottom halo
        sh_par[ixL][blockz-1-izL] = par[i1+iz+blockz-2-1-2*izL];
    }
    __syncthreads();

    if (ix < nx && iz < nz - nc1)
        out[i1+iz] = add*out[i1+iz]
                                    + a/d2 * ( (coef[0]*sh_par[ixL][izL+2]+coef[1]*(sh_par[ixL][izL+1]+sh_par[ixL][izL+3])+coef[2]*(sh_par[ixL][izL]+sh_par[ixL][izL+4]))*sh_in[ixL][izL+2] 
                                    + (coef[3]*sh_par[ixL][izL]+coef[4]*sh_par[ixL][izL+1]+coef[5]*sh_par[ixL][izL+2]+coef[6]*sh_par[ixL][izL+3])*sh_in[ixL][izL+1]
                                    + (coef[3]*sh_par[ixL][izL+4]+coef[4]*sh_par[ixL][izL+3]+coef[5]*sh_par[ixL][izL+2]+coef[6]*sh_par[ixL][izL+1])*sh_in[ixL][izL+3]
                                    + (coef[7]*sh_par[ixL][izL]+coef[8]*sh_par[ixL][izL+1]+coef[9]*sh_par[ixL][izL+2])*sh_in[ixL][izL]
                                    + (coef[7]*sh_par[ixL][izL+4]+coef[8]*sh_par[ixL][izL+3]+coef[9]*sh_par[ixL][izL+2])*sh_in[ixL][izL+4]
                                    );

}

__global__ void cudaDzz_var_bnd(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t * par, data_t a){
	
    int nc2=8, nc3=8;
    data_t  val=0;
    data_t d2 = d*d;

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;
    int iz = threadIdx.y;
    int blockx = BLOCK_SIZE;
    int blockz = 8;  // shared block size with halo
    int ixL = threadIdx.x;
    int i1=ix*nz;

    if (blockIdx.x == gridDim.x-1) blockx = nx-BLOCK_SIZE*(gridDim.x-1); // rightmost block

    __shared__ data_t sh_in[BLOCK_SIZE][8]; // allocate shared memory: BLOCKSIZE + halos
    __shared__ data_t sh_par[BLOCK_SIZE][8];

    if (blockIdx.y == 0) // top boundary
    {
        if (ixL<blockx) {
            sh_in[ixL][iz] = in[i1+iz]; // copy the stencil tile
            sh_par[ixL][iz] = par[i1+iz];
            if (iz<2) {
                sh_in[ixL][blockz-1-iz] = in[i1+blockz-1-iz]; // bottom halo
                sh_par[ixL][blockz-1-iz] = par[i1+blockz-1-iz];
            }
        }
        __syncthreads();

        if (ix<nx)
        {
            out[i1+iz] = add*out[i1+iz];
            int i2=iz*nc2*nc3;
            for (int j = 0; j<nc2; j++){
                val=0;
                for (int k = 0; k<nc3; k++){
                        val += bnd_coef2[i2+j*nc3+k] * sh_par[ixL][k];
                }
                out[i1+iz] += a * val * sh_in[ixL][j] / d2;
            }
        }
    }
    else  // bottom boundary
    {
        if (ixL<blockx) {
            sh_in[ixL][iz] = in[i1+nz-1-iz]; // copy the stencil tile
            sh_par[ixL][iz] = par[i1+nz-1-iz];
            if (iz<2) {
                sh_in[ixL][blockz-1-iz] = in[i1+nz-blockz+iz]; // top halo
                sh_par[ixL][blockz-1-iz] = par[i1+nz-blockz+iz];
            }
        }
        __syncthreads();
        if (ix<nx){
            out[i1+nz-1-iz] = add*out[i1+nz-1-iz];
            int i2=iz*nc2*nc3;
            for (int j = 0; j<nc2; j++){
                    val=0;
                for (int k = 0; k<nc3; k++){
                        val += bnd_coef2[i2+j*nc3+k] * sh_par[ixL][k];
                }
                out[i1+nz-1-iz] +=  a * val * sh_in[ixL][j]/d2;
            }
        }
    }
}

__global__ void cudaDxx_var_interior(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t * par, data_t a){
	
    int nc1=6;
    data_t d2 = d*d;

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE + nc1;
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;
    int blockx = BLOCK_SIZE + 4; // shared block size with halo
    int blockz = BLOCK_SIZE;  
    int ixL = threadIdx.x;
    int izL = threadIdx.y;
    int i1=ix*nz;

    if (blockIdx.x == gridDim.x-1) blockx = nx-2*nc1-BLOCK_SIZE*(gridDim.x-1)+4; // rightmost block
    if (blockIdx.y == gridDim.y-1) blockz = nz-BLOCK_SIZE*(gridDim.y-1); // bottommost block

    __shared__ data_t sh_in[BLOCK_SIZE+4][BLOCK_SIZE]; // allocate shared memory: BLOCKSIZE + halos
    __shared__ data_t sh_par[BLOCK_SIZE+4][BLOCK_SIZE];
    if (ixL<blockx-4 && izL<blockz) {
        sh_in[ixL+2][izL] = in[i1+iz]; // copy the stencil tile
        sh_par[ixL+2][izL] = par[i1+iz];
    }

    if (izL<blockz && ixL<2) {
        sh_in[ixL][izL] = in[(ix-2)*nz+iz]; // left halo
        sh_par[ixL][izL] = par[(ix-2)*nz+iz];
        sh_in[blockx-1-ixL][izL] = in[(ix+blockx-2-1-2*ixL)*nz+iz]; // right halo
        sh_par[blockx-1-ixL][izL] = par[(ix+blockx-2-1-2*ixL)*nz+iz];
    }
    __syncthreads();

    if (ix < nx-nc1 && iz < nz)
        out[ ix* nz+ iz] = add*out[ ix* nz+ iz] + a/d2 * ( (coef[0]* sh_par[ixL+2][izL]+coef[1]*( sh_par[ixL+1][izL]+ sh_par[ixL+3][izL])+coef[2]*( sh_par[ixL][izL]+ sh_par[ixL+4][izL]))*sh_in[ixL+2][izL] 
                        + (coef[3]* sh_par[ixL][izL]+coef[4]* sh_par[ixL+1][izL]+coef[5]* sh_par[ixL+2][izL]+coef[6]* sh_par[ixL+3][izL])*sh_in[ixL+1][izL]
                        + (coef[3]* sh_par[ixL+4][izL]+coef[4]* sh_par[ixL+3][izL]+coef[5]* sh_par[ixL+2][izL]+coef[6]* sh_par[ixL+1][izL])*sh_in[ixL+3][izL]
                        + (coef[7]* sh_par[ixL][izL]+coef[8]* sh_par[ixL+1][izL]+coef[9]* sh_par[ixL+2][izL])*sh_in[ixL][izL]
                        + (coef[7]* sh_par[ixL+4][izL]+coef[8]* sh_par[ixL+3][izL]+coef[9]* sh_par[ixL+2][izL])*sh_in[ixL+4][izL]
                        );
}

__global__ void cudaDxx_var_bnd(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t * par, data_t a){
	
    int nc2=8, nc3=8;
    data_t  val=0;
    data_t d2 = d*d;

    int ix = threadIdx.x;
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;
    int blockx = 8; // shared block size with halo
    int blockz = BLOCK_SIZE;  
    int izL = threadIdx.y;
    int i1=ix*nz;

    if (blockIdx.y == gridDim.y-1) blockz = nz-BLOCK_SIZE*(gridDim.y-1); // bottommost block

    __shared__ data_t sh_in[8][BLOCK_SIZE]; // allocate shared memory: BLOCKSIZE + halos
    __shared__ data_t sh_par[8][BLOCK_SIZE];

    if (blockIdx.x == 0) // left boundary
    {
        if (izL<blockz) {
            sh_in[ix][izL] = in[i1+iz]; // copy the stencil tile
            sh_par[ix][izL] = par[i1+iz];
            if (ix<2) {
                sh_in[blockx-1-ix][izL] = in[(blockx-1-ix)*nz+iz]; // right halo
                sh_par[blockx-1-ix][izL] = par[(blockx-1-ix)*nz+iz];
            }
        }
        __syncthreads();

        if (iz<nz)
        {
            out[i1+iz] = add*out[i1+iz];
            int i2=ix*nc2*nc3;
            for (int j = 0; j<nc2; j++){
                val=0;
                for (int k = 0; k<nc3; k++){
                        val += bnd_coef2[i2+j*nc3+k] * sh_par[k][izL];
                }
                out[i1+iz] += a * val * sh_in[j][izL] / d2;
            }
        }
    }
    else  // right boundary
    {
        if (izL<blockz) {
            sh_in[ix][izL] = in[(nx-1-ix)*nz+iz]; // copy the stencil tile
            sh_par[ix][izL] = par[(nx-1-ix)*nz+iz];
            if (ix<2) {
                sh_in[blockx-1-ix][izL] = in[(nx-blockx+ix)*nz+iz]; // left halo
                sh_par[blockx-1-ix][izL] = par[(nx-blockx+ix)*nz+iz];
            }
        }
        __syncthreads();

        if (iz<nz)
        {
            out[(nx-1-ix)*nz+iz] = add*out[(nx-1-ix)*nz+iz];
            int i2=ix*nc2*nc3;
            for (int j = 0; j<nc2; j++){
                val=0;
                for (int k = 0; k<nc3; k++){
                        val += bnd_coef2[i2+j*nc3+k] * sh_par[k][izL];
                }
                out[(nx-1-ix)*nz+iz] +=  a * val * sh_in[j][izL]/d2;
            }
        }
    }
}

__global__ void cudaEsatNeumannTop(bool add, const data_t* in0, const data_t*in1, data_t* out, int nx, int nz, data_t dx, data_t dz, const data_t * par0, const data_t * par1, data_t a0, data_t a1)
{
    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    data_t sumx=0, sumz=0;

    // top left
    if (ix<nc1)
    {
        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx += bnd_coef1[ix*nc2+j] * in0[j*nz];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in1[ix*nz+iz];
        }
        out[ix*nz] = add*out[ix*nz] - 1.0 /(dz * h0) * (-a0*par0[ix*nz]*sumx/dx + a1*par1[ix*nz]*sumz/dz);
    }

    // top middle
    else if (ix<nx-nc1)
    {
        // (Dx.in1)_0
        sumx=0;
        for (int j=1; j<=2; j++){
            sumx += coef[j-1] * (in0[(ix+j)*nz]-in0[(ix-j)*nz]);
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in1[ix*nz+iz];
        }

        out[ix*nz] = add*out[ix*nz] - 1.0 /(dz * h0) * (-a0*par0[ix*nz]*sumx/dx + a1*par1[ix*nz]*sumz/dz);
    }

    // top right
    else if (ix<nx)
    {
        int jx=nx-ix-1;
        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx -= bnd_coef1[jx*nc2+j] * in0[(nx-1-j)*nz];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in1[(nx-1-jx)*nz+iz];
        }
        out[(nx-1-jx)*nz] = add*out[(nx-1-jx)*nz] - 1.0 /(dz * h0) * (-a0*par0[(nx-1-jx)*nz]*sumx/dx + a1*par1[(nx-1-jx)*nz]*sumz/dz);
    }
}

__global__ void cudaEsatNeumannBottom(bool add, const data_t* in0, const data_t*in1, data_t* out, int nx, int nz, data_t dx, data_t dz, const data_t * par0, const data_t * par1, data_t a0, data_t a1)
{
    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    data_t sumx=0, sumz=0;

    // bottom left
    if (ix<nc1)
    {
        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx += bnd_coef1[ix*nc2+j] * in0[j*nz+nz-1];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in1[ix*nz+nz-1-iz];
        }
        out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - 1.0 /(dz * h0) * (a0*par0[ix*nz+nz-1]*sumx/dx + a1*par1[ix*nz+nz-1]*sumz/dz);
    }

    // bottom middle
    else if (ix<nx-nc1)
    {
        // (Dx.in1)_0
        sumx=0;
        for (int j=1; j<=2; j++){
            sumx += coef[j-1] * (in0[(ix+j)*nz+nz-1]-in0[(ix-j)*nz+nz-1]);
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in1[ix*nz+nz-1-iz];
        }

        out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - 1.0 /(dz * h0) * (a0*par0[ix*nz+nz-1]*sumx/dx + a1*par1[ix*nz+nz-1]*sumz/dz);
    }

    // bottom right
    else if (ix<nx)
    {
        int jx=nx-ix-1;
        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx -= bnd_coef1[jx*nc2+j] * in0[(nx-1-j)*nz+nz-1];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in1[(nx-1-jx)*nz+nz-1-iz];
        }
        out[(nx-1-jx)*nz+nz-1] = add*out[(nx-1-jx)*nz+nz-1] - 1.0 /(dz * h0) * (a0*par0[(nx-1-jx)*nz+nz-1]*sumx/dx + a1*par1[(nx-1-jx)*nz+nz-1]*sumz/dz);
    }
}

__global__ void cudaEsatNeumannLeft(bool add, const data_t* in0, const data_t*in1, data_t* out, int nx, int nz, data_t dx, data_t dz, const data_t * par0, const data_t * par1, data_t a0, data_t a1)
{
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    data_t sumx=0, sumz=0;

    // left top
    if (iz<nc1)
    {
        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz += bnd_coef1[iz*nc2+j] * in0[j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in1[ix*nz+iz];
        }
        out[iz] = add*out[iz] - 1.0 /(dx * h0) * (-a0*par0[iz]*sumz/dz + a1*par1[iz]*sumx/dx);
    }

    // left middle
    else if (iz<nz-nc1)
    {
        // (Dz.in1)_0
        sumz=0;
        for (int j=1; j<=2; j++){
            sumz += coef[j-1] * (in0[iz+j]-in0[iz-j]);
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in1[ix*nz+iz];
        }

        out[iz] = add*out[iz] - 1.0 /(dx * h0) * (-a0*par0[iz]*sumz/dz + a1*par1[iz]*sumx/dx);
    }

    // left bottom
    else if (iz<nz)
    {
        int jz=nz-iz-1;
        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz -= bnd_coef1[jz*nc2+j] * in0[nz-1-j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in1[ix*nz+nz-1-jz];
        }
        out[nz-1-jz] = add*out[nz-1-jz] - 1.0 /(dx * h0) * (-a0*par0[nz-1-jz]*sumz/dz + a1*par1[nz-1-jz]*sumx/dx);
    }
}

__global__ void cudaEsatNeumannRight(bool add, const data_t* in0, const data_t*in1, data_t* out, int nx, int nz, data_t dx, data_t dz, const data_t * par0, const data_t * par1, data_t a0, data_t a1)
{
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    data_t sumx=0, sumz=0;

    // right top
    if (iz<nc1)
    {
        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz += bnd_coef1[iz*nc2+j] * in0[(nx-1)*nz+j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in1[(nx-1-ix)*nz+iz];
        }
        out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - 1.0 /(dx * h0) * (a0*par0[(nx-1)*nz+iz]*sumz/dz + a1*par1[(nx-1)*nz+iz]*sumx/dx);
    }

    // right middle
    else if (iz<nz-nc1)
    {
        // (Dz.in1)_0
        sumz=0;
        for (int j=1; j<=2; j++){
            sumz += coef[j-1] * (in0[(nx-1)*nz+iz+j]-in0[(nx-1)*nz+iz-j]);
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in1[(nx-1-ix)*nz+iz];
        }

        out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - 1.0 /(dx * h0) * (a0*par0[(nx-1)*nz+iz]*sumz/dz + a1*par1[(nx-1)*nz+iz]*sumx/dx);
    }

    // right bottom
    else if (iz<nz)
    {
        int jz=nz-1-iz;
        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz -= bnd_coef1[jz*nc2+j] * in0[(nx-1)*nz+nz-1-j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in1[(nx-1-ix)*nz+nz-1-jz];
        }
        out[(nx-1)*nz+nz-1-jz] = add*out[(nx-1)*nz+nz-1-jz] - 1.0 /(dx * h0) * (a0*par0[(nx-1)*nz+nz-1-jz]*sumz/dz + a1*par1[(nx-1)*nz+nz-1-jz]*sumx/dx);
    }
}

__global__ void cudaEsatAbsorbingTop(bool add, const data_t* in0, const data_t* in1, const data_t* in2, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, const data_t * par0,  const data_t * par1, const data_t * par2, data_t a0, data_t a1, data_t a2)
{
    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    data_t sumx=0, sumz=0;

    // SAT = - H-1 (-f1.Dx.in1 + f2.Sz.in2 -f3.in3/dt)_0
    // f3 is often P or S impedance (rho*Vp or rho*Vs) divided by 2 ; f3=par[2] is a boundary array only

    // top left
    if (ix<nc1)
    {
        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx += bnd_coef1[ix*nc2+j] * in0[j*nz];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in1[ix*nz+iz];
        }
        out[ix*nz] = add*out[ix*nz] - 1.0 /(dz * h0) * (-a0*par0[ix*nz]*sumx/dx + a1*par1[ix*nz]*sumz/dz - a2*par2[ix]*in2[ix*nz]/dt);
    }

    // top middle
    else if (ix<nx-nc1)
    {
        // (Dx.in1)_0
        sumx=0;
        for (int j=1; j<=2; j++){
            sumx += coef[j-1] * (in0[(ix+j)*nz]-in0[(ix-j)*nz]);
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in1[ix*nz+iz];
        }

        out[ix*nz] = add*out[ix*nz] - 1.0 /(dz * h0) * (-a0*par0[ix*nz]*sumx/dx + a1*par1[ix*nz]*sumz/dz - a2*par2[ix]*in2[ix*nz]/dt);
    }

    // top right
    else if (ix<nx)
    {
        int jx=nx-ix-1;
        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx -= bnd_coef1[jx*nc2+j] * in0[(nx-1-j)*nz];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in1[(nx-1-jx)*nz+iz];
        }
        out[(nx-1-jx)*nz] = add*out[(nx-1-jx)*nz] - 1.0 /(dz * h0) * (-a0*par0[(nx-1-jx)*nz]*sumx/dx + a1*par1[(nx-1-jx)*nz]*sumz/dz  - a2*par2[(nx-1-jx)]*in2[(nx-1-jx)*nz]/dt);
    }
}

__global__ void cudaEsatAbsorbingBottom(bool add, const data_t* in0, const data_t* in1, const data_t* in2, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, const data_t * par0,  const data_t * par1, const data_t * par2, data_t a0, data_t a1, data_t a2)
{
    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    data_t sumx=0, sumz=0;

    // SAT = - H-1 (f1.Dx.in1 + f2.Sz.in2 - f3.in3/dt)_0

    // bottom left
    if (ix<nc1)
    {
        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx += bnd_coef1[ix*nc2+j] * in0[j*nz+nz-1];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in1[ix*nz+nz-1-iz];
        }
        out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - 1.0 /(dz * h0) * (a0*par0[ix*nz+nz-1]*sumx/dx + a1*par1[ix*nz+nz-1]*sumz/dz - a2*par2[ix]*in2[ix*nz+nz-1]/dt);
    }

    // bottom middle
    else if (ix<nx-nc1)
    {
        // (Dx.in1)_0
        sumx=0;
        for (int j=1; j<=2; j++){
            sumx += coef[j-1] * (in0[(ix+j)*nz+nz-1]-in0[(ix-j)*nz+nz-1]);
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in1[ix*nz+nz-1-iz];
        }

        out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - 1.0 /(dz * h0) * (a0*par0[ix*nz+nz-1]*sumx/dx + a1*par1[ix*nz+nz-1]*sumz/dz - a2*par2[ix]*in2[ix*nz+nz-1]/dt);
    }

    // bottom right
    else if (ix<nx)
    {
        int jx=nx-ix-1;
        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx -= bnd_coef1[jx*nc2+j] * in0[(nx-1-j)*nz+nz-1];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in1[(nx-1-jx)*nz+nz-1-iz];
        }
        out[(nx-1-jx)*nz+nz-1] = add*out[(nx-1-jx)*nz+nz-1] - 1.0 /(dz * h0) * (a0*par0[(nx-1-jx)*nz+nz-1]*sumx/dx + a1*par1[(nx-1-jx)*nz+nz-1]*sumz/dz - a2*par2[nx-1-jx]*in2[(nx-1-jx)*nz+nz-1]/dt);
    }
}

__global__ void cudaEsatAbsorbingLeft(bool add, const data_t* in0, const data_t* in1, const data_t* in2, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, const data_t * par0,  const data_t * par1, const data_t * par2, data_t a0, data_t a1, data_t a2)
{
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    data_t sumx=0, sumz=0;

    // SAT = - H-1 (-f1.Dz.in1 + f2.Sx.in2 - f3.in3/dt)_0

    // left top
    if (iz<nc1)
    {
        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz += bnd_coef1[iz*nc2+j] * in0[j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in1[ix*nz+iz];
        }
        out[iz] = add*out[iz] - 1.0 /(dx * h0) * (-a0*par0[iz]*sumz/dz + a1*par1[iz]*sumx/dx - a2*par2[iz]*in2[iz]/dt);
    }

    // left middle
    else if (iz<nz-nc1)
    {
        // (Dz.in1)_0
        sumz=0;
        for (int j=1; j<=2; j++){
            sumz += coef[j-1] * (in0[iz+j]-in0[iz-j]);
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in1[ix*nz+iz];
        }

        out[iz] = add*out[iz] - 1.0 /(dx * h0) * (-a0*par0[iz]*sumz/dz + a1*par1[iz]*sumx/dx - a2*par2[iz]*in2[iz]/dt);
    }

    // left bottom
    else if (iz<nz)
    {
        int jz=nz-iz-1;
        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz -= bnd_coef1[jz*nc2+j] * in0[nz-1-j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in1[ix*nz+nz-1-jz];
        }
        out[nz-1-jz] = add*out[nz-1-jz] - 1.0 /(dx * h0) * (-a0*par0[nz-1-jz]*sumz/dz + a1*par1[nz-1-jz]*sumx/dx - a2*par2[nz-1-jz]*in2[nz-1-jz]/dt);
    }
}

__global__ void cudaEsatAbsorbingRight(bool add, const data_t* in0, const data_t* in1, const data_t* in2, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, const data_t * par0,  const data_t * par1, const data_t * par2, data_t a0, data_t a1, data_t a2)
{
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    data_t sumx=0, sumz=0;

    // SAT = - H-1 (f1.Dz.in1 + f2.Sx.in2 - f3.in3/dt)_0

    // right top
    if (iz<nc1)
    {
        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz += bnd_coef1[iz*nc2+j] * in0[(nx-1)*nz+j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in1[(nx-1-ix)*nz+iz];
        }
        out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - 1.0 /(dx * h0) * (a0*par0[(nx-1)*nz+iz]*sumz/dz + a1*par1[(nx-1)*nz+iz]*sumx/dx - a2*par2[iz]*in2[(nx-1)*nz+iz]/dt);
    }

    // right middle
    else if (iz<nz-nc1)
    {
        // (Dz.in1)_0
        sumz=0;
        for (int j=1; j<=2; j++){
            sumz += coef[j-1] * (in0[(nx-1)*nz+iz+j]-in0[(nx-1)*nz+iz-j]);
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in1[(nx-1-ix)*nz+iz];
        }

        out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - 1.0 /(dx * h0) * (a0*par0[(nx-1)*nz+iz]*sumz/dz + a1*par1[(nx-1)*nz+iz]*sumx/dx - a2*par2[iz]*in2[(nx-1)*nz+iz]/dt);
    }

    // right bottom
    else if (iz<nz)
    {
        int jz=nz-1-iz;
        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz -= bnd_coef1[jz*nc2+j] * in0[(nx-1)*nz+nz-1-j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in1[(nx-1-ix)*nz+nz-1-jz];
        }
        out[(nx-1)*nz+nz-1-jz] = add*out[(nx-1)*nz+nz-1-jz] - 1.0 /(dx * h0) * (a0*par0[(nx-1)*nz+nz-1-jz]*sumz/dz + a1*par1[(nx-1)*nz+nz-1-jz]*sumx/dx - a2*par2[nz-1-jz]*in2[(nx-1)*nz+nz-1-jz]/dt);
    }
}

__global__ void cudaDzTop(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, const data_t * par, data_t a){
    
    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;

    data_t h0 = 17.0/48;
    int nc2=6;
    data_t sc = a /(d*d*h0);
    data_t sum=0;

    // additional SAT = - H-1 (-f.Dz.in)_0
    if (ix<nx)
    {
        sum=0;
        for (int j=0; j<nc2; j++){
            sum -= bnd_coef1[j] * in[ix*nz+j];
        }
        out[ix*nz] = add*out[ix*nz] - sc * par[ix*nz]*sum;
    }
}
__global__ void cudaDzBottom(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, const data_t * par, data_t a){
    
    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;

    data_t h0 = 17.0/48;
    int nc2=6;
    data_t sc = a /(d*d*h0);
    data_t sum=0;

    // additional SAT = - H-1 (f.Dz.in)_0
    if (ix<nx)
    {
        sum=0;
        for (int j=0; j<nc2; j++){
            sum -= bnd_coef1[j] * in[ix*nz+nz-1-j];
        }
        out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - sc * par[ix*nz+nz-1]*sum;
    }
}
__global__ void cudaDxLeft(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, const data_t * par, data_t a){
    
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;

    data_t h0 = 17.0/48;
    int nc2=6;
    data_t sc = a /(d*d*h0);
    data_t sum=0;

    // additional SAT = - H-1 (-f.Dx.in)_0
    if (iz<nz)
    {
        sum=0;
        for (int j=0; j<nc2; j++){
            sum -= bnd_coef1[j] * in[j*nz+iz];
        }
        out[iz] = add*out[iz] - sc * par[iz]*sum;
    }
}
__global__ void cudaDxRight(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, const data_t * par, data_t a){
    
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;

    data_t h0 = 17.0/48;
    int nc2=6;
    data_t sc = a /(d*d*h0);
    data_t sum=0;

    // additional SAT = - H-1 (f.Dx.in)_0
    if (iz<nz)
    {
        sum=0;
        for (int j=0; j<nc2; j++){
            sum -= bnd_coef1[j] * in[(nx-1-j)*nz+iz];
        }
        out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - sc * par[(nx-1)*nz+iz]*sum;
    }
}

__global__ void cudaScaleTopBottom(data_t* in, int nx, int nz, data_t dx, data_t dz, const data_t* par, data_t dt, bool top, bool bottom){

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE + 1;

    data_t h0 = 17.0/48;
    data_t scx=0, scz=0;
    const data_t * par0 = par, *par1=par+nx*nz, *par2=par+2*nx*nz;
    data_t * in0 = in, *in1=in+nx*nz;

    if (blockIdx.y==0) // top
    {
        if (ix<nx-1)
        {
            scx = top*sqrt(par1[ix*nz]/par2[ix*nz])*dt / (2 * dz * h0);
            scz = scx*sqrt((par0[ix*nz]+2*par1[ix*nz])/par1[ix*nz]);
            in0[ix*nz] /= (1+scx);
            in1[ix*nz] /= (1+scz);
        }
    }
    else // bottom
    {
        if (ix<nx-1)
        {
            scx = bottom*sqrt(par1[ix*nz+nz-1]/par2[ix*nz+nz-1])*dt / (2 * dz * h0);
            scz = scx*sqrt((par0[ix*nz+nz-1]+2*par1[ix*nz+nz-1])/par1[ix*nz+nz-1]);
            in0[ix*nz+nz-1] /= (1+scx);
            in1[ix*nz+nz-1] /= (1+scz);
        }
    }
}

__global__ void cudaScaleLeftRight(data_t* in, int nx, int nz, data_t dx, data_t dz, const data_t* par, data_t dt, bool top, bool bottom, bool left, bool right){

    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;

    data_t h0 = 17.0/48;
    data_t scx=0, scz=0;
    const data_t * par0 = par, *par1=par+nx*nz, *par2=par+2*nx*nz;
    data_t * in0 = in, *in1=in+nx*nz;

    if (blockIdx.x==0) // left
    {
        if (iz<nz)
        {
            scx = left*sqrt((par0[iz]+2*par1[iz])/par2[iz])*dt / (2 * dx * h0);
            scz = scx*sqrt(par1[iz]/(par0[iz]+2*par1[iz]));
            if (iz==0) { // left top
                scx += top*sqrt(par1[0]/par2[0])*dt / (2  * dz * h0);
                scz += top*sqrt((par0[0]+2*par1[0])/par2[0])*dt / (2  * dz * h0);
            }
            else if (iz==nz-1) { // left bottom
                scx += bottom*sqrt(par1[nz-1]/par2[nz-1])*dt / (2 * dz * h0);
                scz += bottom*sqrt((par0[nz-1]+2*par1[nz-1])/par2[nz-1])*dt / (2 * dz * h0);
            }
            in0[iz] = in0[iz] / (1+scx);
            in1[iz] = in1[iz] / (1+scz);
        }
    }
    else // right
    {
        if (iz<nz)
        {
            scx = right*sqrt((par0[(nx-1)*nz+iz]+2*par1[(nx-1)*nz+iz])/par2[(nx-1)*nz+iz])*dt / (2 * dx * h0);
            scz = scx*sqrt(par1[(nx-1)*nz+iz]/(par0[(nx-1)*nz+iz]+2*par1[(nx-1)*nz+iz]));
            if (iz==0) { // right top
                scx += top*sqrt(par1[(nx-1)*nz]/par2[(nx-1)*nz])*dt / (2 * dz * h0);
                scz += top*sqrt((par0[(nx-1)*nz]+2*par1[(nx-1)*nz])/par2[(nx-1)*nz])*dt / (2 * dz * h0);
            }
            else if (iz==nz-1) { // right bottom
                scx += bottom*sqrt(par1[(nx-1)*nz+nz-1]/par2[(nx-1)*nz+nz-1])*dt / (2 * dz * h0);
                scz += bottom*sqrt((par0[(nx-1)*nz+nz-1]+2*par1[(nx-1)*nz+nz-1])/par2[(nx-1)*nz+nz-1])*dt / (2 * dz * h0);
            }
            in0[(nx-1)*nz+iz] = in0[(nx-1)*nz+iz] / (1+scx);
            in1[(nx-1)*nz+iz] = in1[(nx-1)*nz+iz] / (1+scz);
        }
    }
}

__global__ void cudaScaleLeftRightVTI(data_t* in, int nx, int nz, data_t dx, data_t dz, const data_t* par, data_t dt, bool top, bool bottom, bool left, bool right){

    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;

    data_t h0 = 17.0/48;
    data_t scx=0, scz=0;
    const data_t * par0 = par, *par1=par+nx*nz, *par2=par+2*nx*nz, *par4=par+4*nx*nz;
    data_t * in0 = in, *in1=in+nx*nz;

    if (blockIdx.x==0) // left
    {
        if (iz<nz)
        {
            scx = left*sqrt((1+2*par4[iz])*(par0[iz]+2*par1[iz])/par2[iz])*dt / (2 * dx * h0);
            scz = left*sqrt(par1[iz]/par2[iz])*dt / (2 * dx * h0);
            if (iz==0) { // left top
                scx += top*sqrt(par1[0]/par2[0])*dt / (2  * dz * h0);
                scz += top*sqrt((par0[0]+2*par1[0])/par2[0])*dt / (2  * dz * h0);
            }
            else if (iz==nz-1) { // left bottom
                scx += bottom*sqrt(par1[nz-1]/par2[nz-1])*dt / (2 * dz * h0);
                scz += bottom*sqrt((par0[nz-1]+2*par1[nz-1])/par2[nz-1])*dt / (2 * dz * h0);
            }
            in0[iz] = in0[iz] / (1+scx);
            in1[iz] = in1[iz] / (1+scz);
        }
    }
    else // right
    {
        if (iz<nz)
        {
            scx = right*sqrt((1+2*par4[(nx-1)*nz+iz])*(par0[(nx-1)*nz+iz]+2*par1[(nx-1)*nz+iz])/par2[(nx-1)*nz+iz])*dt / (2 * dx * h0);
            scz = right*sqrt(par1[(nx-1)*nz+iz]/par2[(nx-1)*nz+iz])*dt / (2 * dx * h0);
            if (iz==0) { // right top
                scx += top*sqrt(par1[(nx-1)*nz]/par2[(nx-1)*nz])*dt / (2 * dz * h0);
                scz += top*sqrt((par0[(nx-1)*nz]+2*par1[(nx-1)*nz])/par2[(nx-1)*nz])*dt / (2 * dz * h0);
            }
            else if (iz==nz-1) { // right bottom
                scx += bottom*sqrt(par1[(nx-1)*nz+nz-1]/par2[(nx-1)*nz+nz-1])*dt / (2 * dz * h0);
                scz += bottom*sqrt((par0[(nx-1)*nz+nz-1]+2*par1[(nx-1)*nz+nz-1])/par2[(nx-1)*nz+nz-1])*dt / (2 * dz * h0);
            }
            in0[(nx-1)*nz+iz] = in0[(nx-1)*nz+iz] / (1+scx);
            in1[(nx-1)*nz+iz] = in1[(nx-1)*nz+iz] / (1+scz);
        }
    }
}

__global__ void cudaTaperTop(data_t* in, int nx, int nz, int taper, data_t a){

    int ix = threadIdx.x + blockIdx.x*blockDim.x;
    int iz = threadIdx.y;

    data_t val = cos(a*0.5*M_PI*(taper-1-iz)/taper);
    if (ix<nx) in[ix*nz+iz] *= val*val;
}
__global__ void cudaTaperBottom(data_t* in, int nx, int nz, int taper, data_t a){

    int ix = threadIdx.x + blockIdx.x*blockDim.x;
    int iz = threadIdx.y;

    data_t val = cos(a*0.5*M_PI*(taper-1-iz)/taper);
    if (ix<nx) in[ix*nz+nz-1-iz] *= val*val;
}
__global__ void cudaTaperLeft(data_t* in, int nx, int nz, int taper, data_t a){

    int iz = threadIdx.y + blockIdx.y*blockDim.y;
    int ix = threadIdx.x;

    data_t val = cos(a*0.5*M_PI*(taper-1-ix)/taper);
    if (iz<nz) in[ix*nz+iz] *= val*val;
}
__global__ void cudaTaperRight(data_t* in, int nx, int nz, int taper, data_t a){

    int iz = threadIdx.y + blockIdx.y*blockDim.y;
    int ix = threadIdx.x;

    data_t val = cos(a*0.5*M_PI*(taper-1-ix)/taper);
    if (iz<nz) in[(nx-1-ix)*nz+iz] *= val*val;
}

__global__ void cudaTimeStep(const data_t * prev, const data_t * curr, data_t * next, const data_t * par, int nx, int nz, data_t dt){

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;

    if (ix<nx && iz<nz) next[ix*nz+iz] = dt*dt*next[ix*nz+iz]/par[ix*nz+iz] + 2*curr[ix*nz+iz] - prev[ix*nz+iz];
}

__global__ void cudaInjectDM3(const data_t ** in, data_t * out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw)
{
    int ix = threadIdx.x;
    int iz = threadIdx.y;
    int itr = blockIdx.x+itr_min;

    const int (* p_xind) [3] = (const int (*) [3]) xind;
    const int (* p_zind) [3] = (const int (*) [3]) zind;
    const data_t (* p_xw) [2][3] = (const data_t (*) [2][3]) xw;
    const data_t (* p_zw) [2][3] = (const data_t (*) [2][3]) zw;

    data_t val = p_xw[itr][0][ix] * p_zw[itr][0][iz] * in[0][itr*nt+it];
    atomicAdd(out + p_xind[itr][ix]*nz+p_zind[itr][iz], val);
}

__global__ void cudaExtractDM3(const data_t * in, data_t ** out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw)
{
    int ix = threadIdx.x;
    int iz = threadIdx.y;
    int itr = blockIdx.x+itr_min;

    const int (* p_xind) [3] = (const int (*) [3]) xind;
    const int (* p_zind) [3] = (const int (*) [3]) zind;
    const data_t (* p_xw) [2][3] = (const data_t (*) [2][3]) xw;
    const data_t (* p_zw) [2][3] = (const data_t (*) [2][3]) zw;

    data_t val = p_xw[itr][1][ix] * p_zw[itr][1][iz] * in[p_xind[itr][ix]*nz+p_zind[itr][iz]];
    atomicAdd(out[0]+itr*nt+it, val);
}

__global__ void cudaInjectDDM3(const data_t ** in, data_t * out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw)
{
    int ix = threadIdx.x;
    int iz = threadIdx.y;
    int itr = blockIdx.x+itr_min;

    const int (* p_xind) [3] = (const int (*) [3]) xind;
    const int (* p_zind) [3] = (const int (*) [3]) zind;
    const data_t (* p_xw) [2][6] = (const data_t (*) [2][6]) xw;
    const data_t (* p_zw) [2][6] = (const data_t (*) [2][6]) zw;

    data_t val = (-p_xw[itr][0][ix+3] * p_zw[itr][0][iz] * in[0][itr*nt+it] - p_xw[itr][0][ix] * p_zw[itr][0][iz+3] * in[1][itr*nt+it]);
    atomicAdd(out + p_xind[itr][ix]*nz+p_zind[itr][iz], val);
}

__global__ void cudaExtractDDM3(const data_t * in, data_t ** out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw)
{
    int ix = threadIdx.x;
    int iz = threadIdx.y;
    int itr = blockIdx.x+itr_min;

    const int (* p_xind) [3] = (const int (*) [3]) xind;
    const int (* p_zind) [3] = (const int (*) [3]) zind;
    const data_t (* p_xw) [2][6] = (const data_t (*) [2][6]) xw;
    const data_t (* p_zw) [2][6] = (const data_t (*) [2][6]) zw;

    data_t val = - p_xw[itr][1][ix+3] * p_zw[itr][1][iz] * in[p_xind[itr][ix]*nz+p_zind[itr][iz]];
    atomicAdd(out[0]+itr*nt+it, val);
    val = - p_xw[itr][1][ix] * p_zw[itr][1][iz+3] * in[p_xind[itr][ix]*nz+p_zind[itr][iz]];
    atomicAdd(out[1]+itr*nt+it, val);
}

__global__ void cudaInjectDIPM3(const data_t ** in, data_t * out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw)
{
    int ix = threadIdx.x;
    int iz = threadIdx.y;
    int itr = blockIdx.x+itr_min;

    const int (* p_xind) [6] = (const int (*) [6]) xind;
    const int (* p_zind) [6] = (const int (*) [6]) zind;
    const data_t (* p_xw) [2][6] = (const data_t (*) [2][6]) xw;
    const data_t (* p_zw) [2][6] = (const data_t (*) [2][6]) zw;

    data_t val = -p_xw[itr][0][ix] * p_zw[itr][0][iz] * in[0][itr*nt+it];
    atomicAdd(out + p_xind[itr][ix]*nz+p_zind[itr][iz], val);
    val = p_xw[itr][0][ix+3] * p_zw[itr][0][iz+3] * in[0][itr*nt+it];
    atomicAdd(out + p_xind[itr][ix+3]*nz+p_zind[itr][iz+3], val);
}

__global__ void cudaExtractDIPM3(const data_t * in, data_t ** out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw)
{
    int ix = threadIdx.x;
    int iz = threadIdx.y;
    int itr = blockIdx.x+itr_min;

    const int (* p_xind) [6] = (const int (*) [6]) xind;
    const int (* p_zind) [6] = (const int (*) [6]) zind;
    const data_t (* p_xw) [2][6] = (const data_t (*) [2][6]) xw;
    const data_t (* p_zw) [2][6] = (const data_t (*) [2][6]) zw;

    data_t val = p_xw[itr][1][ix+3] * p_zw[itr][1][iz+3] * in[p_xind[itr][ix+3]*nz+p_zind[itr][iz+3]]
    - p_xw[itr][1][ix] * p_zw[itr][1][iz] * in[p_xind[itr][ix]*nz+p_zind[itr][iz]];
    atomicAdd(out[0]+itr*nt+it, val);
}

__global__ void cudaComputeGradients(const data_t * model, const data_t * u_for, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, int nx, int nz, int nt, int sub, int it, data_t dx, data_t dz, data_t dt){

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;

    int nxz=nx*nz;
    const data_t * pfor0x=u_for;
    const data_t * pfor0z=u_for+nxz;
    const data_t * pfor1x=u_for+2*nxz;
    const data_t * pfor1z=u_for+3*nxz;
    const data_t * pfor2x=u_for+4*nxz;
    const data_t * pfor2z=u_for+5*nxz;
    const data_t *padjx = curr, *padjz=curr+nxz;
    const data_t *padjx_x = u_x, *padjz_x=u_x+nxz;
    const data_t *padjx_z = u_z, *padjz_z=u_z+nxz;
    data_t *pforx_x=tmp, *pforz_z=tmp+nxz, *pforz_x=tmp+2*nxz, *pforx_z=tmp+3*nxz;
    data_t *gla = grad, *gmu=grad+nxz, *grho=grad+2*nxz;

    int i=0;
    if (ix<nx && iz<nz)
    {
        i=ix*nz+iz;
        if (nt/sub+1-it-1>0){
            gla[i] += dt*(padjx_x[i] + padjz_z[i])*(pforx_x[i] + pforz_z[i]);
            gmu[i] += dt*((padjx_z[i] + padjz_x[i])*(pforz_x[i] + pforx_z[i]) + 2*padjx_x[i]*pforx_x[i] + 2*padjz_z[i]*pforz_z[i]);
            grho[i] += 1.0/dt*(padjx[i]*(pfor2x[i]-2*pfor1x[i]+pfor0x[i]) + padjz[i]*(pfor2z[i]-2*pfor1z[i]+pfor0z[i]));
        }
        else{
            gla[i] += 0.5*dt*(padjx_x[i] + padjz_z[i])*(pforx_x[i] + pforz_z[i]);
            gmu[i] += 0.5*dt*((padjx_z[i] + padjz_x[i])*(pforz_x[i] + pforx_z[i]) + 2*padjx_x[i]*pforx_x[i] + 2*padjz_z[i]*pforz_z[i]);
            grho[i] += 1.0/dt*(padjx[i]*(pfor2x[i]-pfor1x[i]) + padjz[i]*(pfor2z[i]-pfor1z[i]));
        }
    }
}

__global__ void cudaComputeGradientsVTI(const data_t * model, const data_t * u_for, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, int nx, int nz, int nt, int sub, int it, data_t dx, data_t dz, data_t dt){

    int ix = threadIdx.x + blockIdx.x*BLOCK_SIZE;
    int iz = threadIdx.y + blockIdx.y*BLOCK_SIZE;

    int nxz=nx*nz;
    const data_t * pfor0x=u_for;
    const data_t * pfor0z=u_for+nxz;
    const data_t * pfor1x=u_for+2*nxz;
    const data_t * pfor1z=u_for+3*nxz;
    const data_t * pfor2x=u_for+4*nxz;
    const data_t * pfor2z=u_for+5*nxz;
    const data_t *padjx = curr, *padjz=curr+nxz;
    const data_t *padjx_x = u_x, *padjz_x=u_x+nxz;
    const data_t *padjx_z = u_z, *padjz_z=u_z+nxz;
    data_t *pforx_x=tmp, *pforz_z=tmp+nxz, *pforz_x=tmp+2*nxz, *pforx_z=tmp+3*nxz;
    data_t *gla = grad, *gmu=grad+nxz, *grho=grad+2*nxz, *gc13=grad+3*nxz, *geps=grad+4*nxz;
    const data_t *pm0 = model, *pm1=model+nxz, *pm3=model+3*nxz, *pm4=model+4*nxz;

    data_t val1=0, val2=0, val3=0, del=0;
    int i=0;
    if (ix<nx && iz<nz)
    {
        i=ix*nz+iz;
        del = ((pm3[i]+pm1[i])*(pm3[i]+pm1[i]) - (pm0[i]+pm1[i])*(pm0[i]+pm1[i])) / (2*(pm0[i]+pm1[i])*(pm0[i]+2*pm1[i]));((pm3[i]+pm1[i])*(pm3[i]+pm1[i]) - (pm0[i]+pm1[i])*(pm0[i]+pm1[i])) / (2*(pm0[i]+pm1[i])*(pm0[i]+2*pm1[i]));
        val1 = sqrt(2*(pm0[i]+2*pm1[i])*(pm0[i]+pm1[i])*del + (pm0[i]+pm1[i])*(pm0[i]+pm1[i]));
        val2 = ((1+2*del)*pm0[i] + (1+3*del)*pm1[i])/val1; // d(C13)/d(lambda)
        val3 = ((1+3*del)*pm0[i] + (1+4*del)*pm1[i])/val1 - 1; // d(C13)/d(mu)
        if (nt/sub+1-it-1>0){
            gla[i] += dt*((1+2*pm4[i])*padjx_x[i]*pforx_x[i] + padjz_z[i]*pforz_z[i] + val2*(padjx_x[i]*pforz_z[i] + padjz_z[i]*pforx_x[i])); // lambda gradient
            gmu[i] += dt*((padjx_z[i] + padjz_x[i])*(pforz_x[i] + pforx_z[i]) + 2*(1+2*pm4[i])*padjx_x[i]*pforx_x[i] + 2*padjz_z[i]*pforz_z[i] + val3*(padjx_x[i]*pforz_z[i] + padjz_z[i]*pforx_x[i])); // mu gradient
            grho[i] += 1.0/dt*(padjx[i]*(pfor2x[i]-2*pfor1x[i]+pfor0x[i]) + padjz[i]*(pfor2z[i]-2*pfor1z[i]+pfor0z[i])); // rho gradient
            gc13[i] = 0;
            geps[i] = 0;
        }
        else{
            gla[i] += 0.5*dt*((1+2*pm4[i])*padjx_x[i]*pforx_x[i] + padjz_z[i]*pforz_z[i] + val2*(padjx_x[i]*pforz_z[i] + padjz_z[i]*pforx_x[i])); // lambda gradient
            gmu[i] += 0.5*dt*((padjx_z[i] + padjz_x[i])*(pforz_x[i] + pforx_z[i]) + 2*(1+2*pm4[i])*padjx_x[i]*pforx_x[i] + 2*padjz_z[i]*pforz_z[i] + val3*(padjx_x[i]*pforz_z[i] + padjz_z[i]*pforx_x[i])); // mu gradient
            grho[i] += 1.0/dt*(padjx[i]*(pfor2x[i]-pfor1x[i]) + padjz[i]*(pfor2z[i]-pfor1z[i])); // rho gradient
            gc13[i] = 0;
            geps[i] = 0;
        }
    }
}

// ################################ C++ wrappers ##################################

cudaStream_t streams[5];

void Dz_gpu(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int stream1, int stream2){
    
    dim3 threads1(BLOCK_SIZE,BLOCK_SIZE);
    dim3 blocks1((nx+BLOCK_SIZE-1)/BLOCK_SIZE,(nz-8+BLOCK_SIZE-1)/BLOCK_SIZE);
    dim3 threads2(BLOCK_SIZE,4);
    dim3 blocks2((nx+BLOCK_SIZE-1)/BLOCK_SIZE,2);

    cudaDz_interior<<<blocks1,threads1,0,streams[stream1]>>>(add, in, out, nx, nz, d);
    cudaKernelError();
    cudaDz_bnd<<<blocks2,threads2,0,streams[stream2]>>>(add, in, out, nx, nz, d);
    cudaKernelError();
}

void Dx_gpu(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int stream1, int stream2){

    dim3 threads1(BLOCK_SIZE,BLOCK_SIZE);
    dim3 blocks1((nx-8+BLOCK_SIZE-1)/BLOCK_SIZE,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);
    dim3 threads2(4,BLOCK_SIZE);
    dim3 blocks2(2,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);

    cudaDx_interior<<<blocks1,threads1,0,streams[stream1]>>>(add, in, out, nx, nz, d);
    cudaKernelError();
    cudaDx_bnd<<<blocks2,threads2,0,streams[stream2]>>>(add, in, out, nx, nz, d);
    cudaKernelError();
}

void mult_Dz_gpu(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t* par, data_t a, int stream1, int stream2){
    
    dim3 threads1(BLOCK_SIZE,BLOCK_SIZE);
    dim3 blocks1((nx+BLOCK_SIZE-1)/BLOCK_SIZE,(nz-8+BLOCK_SIZE-1)/BLOCK_SIZE);
    dim3 threads2(BLOCK_SIZE,4);
    dim3 blocks2((nx+BLOCK_SIZE-1)/BLOCK_SIZE,2);

    cudaMultDz_interior<<<blocks1,threads1,0,streams[stream1]>>>(add, in, out, nx, nz, d, par, a);
    cudaKernelError();
    cudaMultDz_bnd<<<blocks2,threads2,0,streams[stream2]>>>(add, in, out, nx, nz, d, par, a);
    cudaKernelError();
}

void mult_Dx_gpu(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t * par, data_t a, int stream1, int stream2){

    dim3 threads1(BLOCK_SIZE,BLOCK_SIZE);
    dim3 blocks1((nx-8+BLOCK_SIZE-1)/BLOCK_SIZE,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);
    dim3 threads2(4,BLOCK_SIZE);
    dim3 blocks2(2,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);

    cudaMultDx_interior<<<blocks1,threads1,0,streams[stream1]>>>(add, in, out, nx, nz, d, par, a);
    cudaKernelError();
    cudaMultDx_bnd<<<blocks2,threads2,0,streams[stream2]>>>(add, in, out, nx, nz, d, par, a);
    cudaKernelError();
}

void Dzz_var_gpu(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t * par, data_t a, int stream1, int stream2){

    dim3 threads1(BLOCK_SIZE,BLOCK_SIZE);
    dim3 blocks1((nx+BLOCK_SIZE-1)/BLOCK_SIZE,(nz-12+BLOCK_SIZE-1)/BLOCK_SIZE);
    dim3 threads2(BLOCK_SIZE,6);
    dim3 blocks2((nx+BLOCK_SIZE-1)/BLOCK_SIZE,2);

    cudaDzz_var_interior<<<blocks1,threads1,0,streams[stream1]>>>(add, in, out, nx, nz, d, par, a);
    cudaKernelError();
    cudaDzz_var_bnd<<<blocks2,threads2,0,streams[stream2]>>>(add, in, out, nx, nz, d, par, a);
    cudaKernelError();
}

void Dxx_var_gpu(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, const data_t * par, data_t a, int stream1, int stream2){

    dim3 threads1(BLOCK_SIZE,BLOCK_SIZE);
    dim3 blocks1((nx-12+BLOCK_SIZE-1)/BLOCK_SIZE,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);
    dim3 threads2(6,BLOCK_SIZE);
    dim3 blocks2(2,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);

    cudaDxx_var_interior<<<blocks1,threads1,0,streams[stream1]>>>(add, in, out, nx, nz, d, par, a);
    cudaKernelError();
    cudaDxx_var_bnd<<<blocks2,threads2,0,streams[stream2]>>>(add, in, out, nx, nz, d, par, a);
    cudaKernelError();
}

void esat_neumann_top_gpu(bool add, const data_t* in0, const data_t*in1, data_t* out, int nx, int nz, data_t dx, data_t dz, const data_t * par0, const data_t * par1, data_t a0, data_t a1, int stream){
    
    dim3 threads(BLOCK_SIZE,1);
    dim3 blocks((nx+BLOCK_SIZE-1)/BLOCK_SIZE,1);

    cudaEsatNeumannTop<<<blocks,threads,0,streams[stream]>>>(add, in0, in1, out, nx, nz, dx, dz, par0, par1, a0, a1);
    cudaKernelError();
}
void esat_neumann_bottom_gpu(bool add, const data_t* in0, const data_t*in1, data_t* out, int nx, int nz, data_t dx, data_t dz, const data_t * par0, const data_t * par1, data_t a0, data_t a1, int stream){
    
    dim3 threads(BLOCK_SIZE,1);
    dim3 blocks((nx+BLOCK_SIZE-1)/BLOCK_SIZE,1);

    cudaEsatNeumannBottom<<<blocks,threads,0,streams[stream]>>>(add, in0, in1, out, nx, nz, dx, dz, par0, par1, a0, a1);
    cudaKernelError();
}
void esat_neumann_left_gpu(bool add, const data_t* in0, const data_t*in1, data_t* out, int nx, int nz, data_t dx, data_t dz, const data_t * par0, const data_t * par1, data_t a0, data_t a1, int stream){
    
    dim3 threads(1,BLOCK_SIZE);
    dim3 blocks(1,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);

    cudaEsatNeumannLeft<<<blocks,threads,0,streams[stream]>>>(add, in0, in1, out, nx, nz, dx, dz, par0, par1, a0, a1);
    cudaKernelError();
}
void esat_neumann_right_gpu(bool add, const data_t* in0, const data_t*in1, data_t* out, int nx, int nz, data_t dx, data_t dz, const data_t * par0, const data_t * par1, data_t a0, data_t a1, int stream){
    
    dim3 threads(1,BLOCK_SIZE);
    dim3 blocks(1,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);

    cudaEsatNeumannRight<<<blocks,threads,0,streams[stream]>>>(add, in0, in1, out, nx, nz, dx, dz, par0, par1, a0, a1);
    cudaKernelError();
}

void esat_absorbing_top_gpu(bool add, const data_t* in0, const data_t* in1, const data_t* in2, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, const data_t * par0,  const data_t * par1, const data_t * par2, data_t a0, data_t a1, data_t a2, int stream){
    
    dim3 threads(BLOCK_SIZE,1);
    dim3 blocks((nx+BLOCK_SIZE-1)/BLOCK_SIZE,1);

    cudaEsatAbsorbingTop<<<blocks,threads,0,streams[stream]>>>(add, in0, in1, in2, out, nx, nz, dx, dz, dt, par0, par1, par2, a0, a1, a2);
    cudaKernelError();
}
void esat_absorbing_bottom_gpu(bool add, const data_t* in0, const data_t* in1, const data_t* in2, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, const data_t * par0,  const data_t * par1, const data_t * par2, data_t a0, data_t a1, data_t a2, int stream){
    
    dim3 threads(BLOCK_SIZE,1);
    dim3 blocks((nx+BLOCK_SIZE-1)/BLOCK_SIZE,1);

    cudaEsatAbsorbingBottom<<<blocks,threads,0,streams[stream]>>>(add, in0, in1, in2, out, nx, nz, dx, dz, dt, par0, par1, par2, a0, a1, a2);
    cudaKernelError();
}
void esat_absorbing_left_gpu(bool add, const data_t* in0, const data_t* in1, const data_t* in2, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, const data_t * par0,  const data_t * par1, const data_t * par2, data_t a0, data_t a1, data_t a2, int stream){
    
    dim3 threads(1,BLOCK_SIZE);
    dim3 blocks(1,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);

    cudaEsatAbsorbingLeft<<<blocks,threads,0,streams[stream]>>>(add, in0, in1, in2, out, nx, nz, dx, dz, dt, par0, par1, par2, a0, a1, a2);
    cudaKernelError();
}
void esat_absorbing_right_gpu(bool add, const data_t* in0, const data_t* in1, const data_t* in2, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, const data_t * par0,  const data_t * par1, const data_t * par2, data_t a0, data_t a1, data_t a2, int stream){
    
    dim3 threads(1,BLOCK_SIZE);
    dim3 blocks(1,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);

    cudaEsatAbsorbingRight<<<blocks,threads,0,streams[stream]>>>(add, in0, in1, in2, out, nx, nz, dx, dz, dt, par0, par1, par2, a0, a1, a2);
    cudaKernelError();
}

void esat_Dz_top_gpu(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, const data_t * par, data_t a, int stream){ 
    dim3 threads(BLOCK_SIZE,1);
    dim3 blocks((nx+BLOCK_SIZE-1)/BLOCK_SIZE,1);

    cudaDzTop<<<blocks,threads,0,streams[stream]>>>(add, in, out, nx, nz, d, par, a);
    cudaKernelError();
}
void esat_Dz_bottom_gpu(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, const data_t * par, data_t a, int stream){
    dim3 threads(BLOCK_SIZE,1);
    dim3 blocks((nx+BLOCK_SIZE-1)/BLOCK_SIZE,1);

    cudaDzBottom<<<blocks,threads,0,streams[stream]>>>(add, in, out, nx, nz, d, par, a);
    cudaKernelError();
}
void esat_Dx_left_gpu(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, const data_t * par, data_t a, int stream){
    dim3 threads(1,BLOCK_SIZE);
    dim3 blocks(1,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);

    cudaDxLeft<<<blocks,threads,0,streams[stream]>>>(add, in, out, nx, nz, d, par, a);
    cudaKernelError();
}
void esat_Dx_right_gpu(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, const data_t * par, data_t a, int stream){
    dim3 threads(1,BLOCK_SIZE);
    dim3 blocks(1,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);

    cudaDxRight<<<blocks,threads,0,streams[stream]>>>(add, in, out, nx, nz, d, par, a);
    cudaKernelError();
}

void esat_scale_boundaries_gpu(data_t* in, int nx, int nz, data_t dx, data_t dz, const data_t* par, data_t dt, bool top, bool bottom, bool left, bool right, int stream1, int stream2){
    
    dim3 threads1(BLOCK_SIZE,1); // top-bottom excluding the corners
    dim3 blocks1((nx-2+BLOCK_SIZE-1)/BLOCK_SIZE,2);
    dim3 threads2(1,BLOCK_SIZE); // left-right 
    dim3 blocks2(2,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);
    
    if  (top || bottom) {
        cudaScaleTopBottom<<<blocks1,threads1,0,streams[stream1]>>>(in, nx, nz, dx, dz, par, dt, top, bottom);
        cudaKernelError();
    }
    if  (top || bottom || left || right) {
        cudaScaleLeftRight<<<blocks2,threads2,0,streams[stream2]>>>(in, nx, nz, dx, dz, par, dt, top, bottom, left, right);
        cudaKernelError();
    }
}
void vtisat_scale_boundaries_gpu(data_t* in, int nx, int nz, data_t dx, data_t dz, const data_t* par, data_t dt, bool top, bool bottom, bool left, bool right, int stream1, int stream2){

    dim3 threads1(BLOCK_SIZE,1); // top-bottom excluding the corners
    dim3 blocks1((nx-2+BLOCK_SIZE-1)/BLOCK_SIZE,2);
    dim3 threads2(1,BLOCK_SIZE); // left-right 
    dim3 blocks2(2,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);
    
    if  (top || bottom) {
        cudaScaleTopBottom<<<blocks1,threads1,0,streams[stream1]>>>(in, nx, nz, dx, dz, par, dt, top, bottom);
        cudaKernelError();
    }
    if  (top || bottom || left || right) {
        cudaScaleLeftRightVTI<<<blocks2,threads2,0,streams[stream2]>>>(in, nx, nz, dx, dz, par, dt, top, bottom, left, right);
        cudaKernelError();
    }
}

void taper_top_gpu(data_t* in, int nx, int nz, int taper, data_t a, int stream){

    dim3 threads(BLOCK_SIZE/2,taper);
    dim3 blocks(2*(nx+BLOCK_SIZE/2-1)/BLOCK_SIZE,1);
    cudaTaperTop<<<blocks,threads,0,streams[stream]>>>(in,nx,nz,taper,a);
    cudaKernelError();
}

void taper_bottom_gpu(data_t* in, int nx, int nz, int taper, data_t a, int stream){

    dim3 threads(BLOCK_SIZE/2,taper);
    dim3 blocks(2*(nx+BLOCK_SIZE/2-1)/BLOCK_SIZE,1);
    cudaTaperBottom<<<blocks,threads,0,streams[stream]>>>(in,nx,nz,taper,a);
    cudaKernelError();
}

void taper_left_gpu(data_t* in, int nx, int nz, int taper, data_t a, int stream){

    dim3 threads(taper,BLOCK_SIZE/2);
    dim3 blocks(1,2*(nx+BLOCK_SIZE/2-1)/BLOCK_SIZE);
    cudaTaperLeft<<<blocks,threads,0,streams[stream]>>>(in,nx,nz,taper,a);
    cudaKernelError();
}

void taper_right_gpu(data_t* in, int nx, int nz, int taper, data_t a, int stream){

    dim3 threads(taper,BLOCK_SIZE/2);
    dim3 blocks(1,2*(nx+BLOCK_SIZE/2-1)/BLOCK_SIZE);
    cudaTaperRight<<<blocks,threads,0,streams[stream]>>>(in,nx,nz,taper,a);
    cudaKernelError();
}

void time_step_gpu(const data_t * prev, const data_t * curr, data_t * next, const data_t * par, int nx, int nz, data_t dt, int stream){

    dim3 threads(BLOCK_SIZE,BLOCK_SIZE);
    dim3 blocks((nx+BLOCK_SIZE-1)/BLOCK_SIZE,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);
    cudaTimeStep<<<blocks,threads,0,streams[stream]>>>(prev,curr,next,par,nx,nz,dt);
    cudaKernelError();
}

void delta_m3::inject_gpu(bool add, const data_t ** in, data_t * out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw) const{

    if (add==false) cudaMemset(out,0,nx*nz*sizeof(data_t));
    dim3 threads(3,3);
    dim3 blocks(itr_max-itr_min,1);
    cudaInjectDM3<<<blocks,threads,0,streams[1]>>>(in, out, nx, nz, nt, ntr, it, itr_min, itr_max, xind, zind, xw, zw);
    cudaKernelError();
}

void delta_m3::extract_gpu(bool add, const data_t * in, data_t ** out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw) const{

    if (add==false) cudaMemset(out[0]+itr_min,0,nt*(itr_max-itr_min)*sizeof(data_t));
    dim3 threads(3,3);
    dim3 blocks(itr_max-itr_min,1);
    cudaExtractDM3<<<blocks,threads,0,streams[1]>>>(in, out, nx, nz, nt, ntr, it, itr_min, itr_max, xind, zind, xw, zw);
    cudaKernelError();
}

void ddelta_m3::inject_gpu(bool add, const data_t ** in, data_t * out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw) const{

    if (add==false) cudaMemset(out,0,nx*nz*sizeof(data_t));
    dim3 threads(3,3);
    dim3 blocks(itr_max-itr_min,1);
    cudaInjectDDM3<<<blocks,threads,0,streams[1]>>>(in, out, nx, nz, nt, ntr, it, itr_min, itr_max, xind, zind, xw, zw);
    cudaKernelError();
}

void ddelta_m3::extract_gpu(bool add, const data_t * in, data_t ** out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw) const{

    if (add==false) {
        cudaMemset(out[0]+itr_min,0,nt*(itr_max-itr_min)*sizeof(data_t));
        cudaMemset(out[1]+itr_min,0,nt*(itr_max-itr_min)*sizeof(data_t));
    }
    dim3 threads(3,3);
    dim3 blocks(itr_max-itr_min,1);
    cudaExtractDDM3<<<blocks,threads,0,streams[1]>>>(in, out, nx, nz, nt, ntr, it, itr_min, itr_max, xind, zind, xw, zw);
    cudaKernelError();
}

void dipole_m3::inject_gpu(bool add, const data_t ** in, data_t * out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw) const{

    if (add==false) cudaMemset(out,0,nx*nz*sizeof(data_t));
    dim3 threads(3,3);
    dim3 blocks(itr_max-itr_min,1);
    cudaInjectDIPM3<<<blocks,threads,0,streams[1]>>>(in, out, nx, nz, nt, ntr, it, itr_min, itr_max, xind, zind, xw, zw);
    cudaKernelError();
}

void dipole_m3::extract_gpu(bool add, const data_t * in, data_t ** out, int nx, int nz, int nt, int ntr, int it, int itr_min, int itr_max, const int * xind, const int * zind, const data_t * xw, const data_t * zw) const{

    if (add==false) cudaMemset(out[0]+itr_min,0,nt*(itr_max-itr_min)*sizeof(data_t));
    dim3 threads(3,3);
    dim3 blocks(itr_max-itr_min,1);
    cudaExtractDIPM3<<<blocks,threads,0,streams[1]>>>(in, out, nx, nz, nt, ntr, it, itr_min, itr_max, xind, zind, xw, zw);
    cudaKernelError();
}

void nl_we_op_e::compute_gradients_gpu(const data_t * model, const data_t * u_for, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const {
    
    int nxz = nx*nz;
    const data_t * pfor1x=u_for+2*nxz;
    const data_t * pfor1z=u_for+3*nxz;

    Dx_gpu(false, pfor1x, tmp, nx, nz, dx, 3, 4); // forwardx_x
    Dz_gpu(false, pfor1z, tmp+nxz, nx, nz, dz, 4, 3); // forwardz_z
    Dx_gpu(false, pfor1z, tmp+2*nxz, nx, nz, dx, 3, 4); // forwardz_x
    Dz_gpu(false, pfor1x, tmp+3*nxz, nx, nz, dz, 4, 3); // forwardx_z
    cudaStreamSynchronize(streams[3]);
    cudaStreamSynchronize(streams[4]);

    dim3 threads(BLOCK_SIZE,BLOCK_SIZE);
    dim3 blocks((nx+BLOCK_SIZE-1)/BLOCK_SIZE,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);
    cudaComputeGradients<<<blocks,threads,0,streams[3]>>>(model, u_for, curr, u_x, u_z, tmp, grad, nx, nz, par.nt, par.sub, it, dx, dz, dt);
    cudaKernelError();
}

void nl_we_op_vti::compute_gradients_gpu(const data_t * model, const data_t * u_for, const data_t * curr, const data_t * u_x, const data_t * u_z, data_t * tmp, data_t * grad, const param &par, int nx, int nz, int it, data_t dx, data_t dz, data_t dt) const {
    
    int nxz = nx*nz;
    const data_t * pfor1x=u_for+2*nxz;
    const data_t * pfor1z=u_for+3*nxz;

    Dx_gpu(false, pfor1x, tmp, nx, nz, dx, 3, 4); // forwardx_x
    Dz_gpu(false, pfor1z, tmp+nxz, nx, nz, dz, 4, 3); // forwardz_z
    Dx_gpu(false, pfor1z, tmp+2*nxz, nx, nz, dx, 3, 4); // forwardz_x
    Dz_gpu(false, pfor1x, tmp+3*nxz, nx, nz, dz, 4, 3); // forwardx_z
    cudaStreamSynchronize(streams[3]);
    cudaStreamSynchronize(streams[4]);

    dim3 threads(BLOCK_SIZE,BLOCK_SIZE);
    dim3 blocks((nx+BLOCK_SIZE-1)/BLOCK_SIZE,(nz+BLOCK_SIZE-1)/BLOCK_SIZE);
    cudaComputeGradientsVTI<<<blocks,threads,0,streams[3]>>>(model, u_for, curr, u_x, u_z, tmp, grad, nx, nz, par.nt, par.sub, it, dx, dz, dt);
    cudaKernelError();
}


#undef BLOCK_SIZE