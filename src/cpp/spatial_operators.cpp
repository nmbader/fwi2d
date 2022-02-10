#include <string.h>
#include <algorithm>
#include <omp.h>
#include "spatial_operators.hpp"

void applyHz(bool inv, bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax){

    int izminb=std::max(4,izmin);
    int izmaxb=std::min(nz - 4, izmax);

    data_t coef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};

    if (!inv){
        #pragma omp parallel for
        for (int ix=ixmin; ix<ixmax; ix++){
            int i1=ix*nz;
            for (int iz = izmin; iz<4; iz++){
                out[i1+iz] = add*out[i1+iz] + in[i1+iz]*coef[iz]*d;
            }
            for (int iz = nz-izmax; iz<4; iz++){
                out[i1+nz-1-iz] = add*out[i1+nz-1-iz] + in[i1+nz-1-iz]*coef[iz]*d;
            }
            for (int iz=izminb; iz<izmaxb; iz++){
                out[i1+iz] = add*out[i1+iz] + in[i1+iz]*d;
            }
        }
    }
    else{
        #pragma omp parallel for
        for (int ix=ixmin; ix<ixmax; ix++){
            int i1=ix*nz;
            for (int iz = izmin; iz<4; iz++){
                out[i1+iz] = add*out[i1+iz] + in[i1+iz]/(coef[iz]*d);
            }
            for (int iz = nz-izmax; iz<4; iz++){
                out[i1+nz-1-iz] = add*out[i1+nz-1-iz] + in[i1+nz-1-iz]/(coef[iz]*d);
            }
            for (int iz=izminb; iz<izmaxb; iz++){
                out[i1+iz] = add*out[i1+iz] + in[i1+iz]/d;
            }
        }
    }
}

void applyHx(bool inv, bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax){

    int ixminb=std::max(4,ixmin);
    int ixmaxb=std::min(nx - 4, ixmax);

    data_t coef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};

    if (!inv){
        #pragma omp parallel for
        for (int ix = ixmin; ix<4; ix++){
            int i1=ix*nz;
            for (int iz=izmin; iz<izmax; iz++){
                out[i1+iz] = add*out[i1+iz] + in[i1+iz]*coef[ix]*d;
            }
        }
        #pragma omp parallel for
        for (int ix = nx-ixmax; ix<4; ix++){
            int i1=(nx-1-ix)*nz;
            for (int iz=izmin; iz<izmax; iz++){
                out[i1+iz] = add*out[i1+iz] + in[i1+iz]*coef[ix]*d;
            }
        }
        #pragma omp parallel for
        for (int ix=ixminb; ix<ixmaxb; ix++){
            int i1=ix*nz;
            for (int iz=izmin; iz<izmax; iz++){
                out[i1+iz] = add*out[i1+iz] + in[i1+iz]*d;
            }
        }
    }
    else{
        #pragma omp parallel for
        for (int ix = ixmin; ix<4; ix++){
            int i1=ix*nz;
            for (int iz=izmin; iz<izmax; iz++){
                out[i1+iz] = add*out[i1+iz] + in[i1+iz]/(coef[ix]*d);
            }
        }
        #pragma omp parallel for
        for (int ix = nx-ixmax; ix<4; ix++){
            int i1=(nx-1-ix)*nz;
            for (int iz=izmin; iz<izmax; iz++){
                out[i1+iz] = add*out[i1+iz] + in[i1+iz]/(coef[ix]*d);
            }
        }
        #pragma omp parallel for
        for (int ix=ixminb; ix<ixmaxb; ix++){
            int i1=ix*nz;
            for (int iz=izmin; iz<izmax; iz++){
                out[i1+iz] = add*out[i1+iz] + in[i1+iz]/d;
            }
        }
    }
}

void Dz(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax){

    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    
    int nc1=4, nc2=6;
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);

    #pragma omp parallel for
    for (int ix = ixmin; ix < ixmax; ix++){
        int i1=ix*nz;

        // apply the operator near the top boundary if included
        for (int iz=izmin; iz<nc1; iz++){
             out[i1+iz] = add* out[i1+iz] + (bnd_coef[iz*nc2] * in[i1] + bnd_coef[iz*nc2+1] * in[i1+1] + bnd_coef[iz*nc2+2] * in[i1+2] + bnd_coef[iz*nc2+3] * in[i1+3] + bnd_coef[iz*nc2+4] * in[i1+4] + bnd_coef[iz*nc2+5] * in[i1+5]) / d;
        }

        // apply the operator near the bottom boundary if included
        for (int iz=nz-izmax; iz<nc1; iz++){
             out[i1+nz-1-iz] = add* out[i1+nz-1-iz] + (-bnd_coef[iz*nc2] * in[i1+nz-1] - bnd_coef[iz*nc2+1] * in[i1+nz-2] - bnd_coef[iz*nc2+2] * in[i1+nz-3] - bnd_coef[iz*nc2+3] * in[i1+nz-4] - bnd_coef[iz*nc2+4] * in[i1+nz-5] - bnd_coef[iz*nc2+5] * in[i1+nz-6]) / d;
        }

        // apply the operator in the interior
        Dz_ispc(add, in, out, 1.0/d,  i1,  izminb,  izmaxb);
    }
}

void Dx(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax){
    
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    
    int nc1=4, nc2=6;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);

    // apply the operator near the left boundary if included
    #pragma omp parallel for
    for (int ix=ixmin; ix<nc1; ix++){
        int i1=ix*nz;
        for (int iz = izmin; iz < izmax; iz++) out[i1+iz] = add*out[i1+iz] + (bnd_coef[ix*nc2] * in[iz] + bnd_coef[ix*nc2+1] * in[nz+iz] + bnd_coef[ix*nc2+2] * in[2*nz+iz] + bnd_coef[ix*nc2+3] * in[3*nz+iz] + bnd_coef[ix*nc2+4] * in[4*nz+iz] + bnd_coef[ix*nc2+5] * in[5*nz+iz]) / d; 
    }

    // apply the operator near the right boundary if included
    #pragma omp parallel for
    for (int ix=nx-ixmax; ix<nc1; ix++){
        int i1=(nx-1-ix)*nz;
        for (int iz = izmin; iz < izmax; iz++) out[i1+iz] = add*out[i1+iz] + (-bnd_coef[ix*nc2] * in[(nx-1)*nz+iz] - bnd_coef[ix*nc2+1] * in[(nx-2)*nz+iz] - bnd_coef[ix*nc2+2] * in[(nx-3)*nz+iz] - bnd_coef[ix*nc2+3] * in[(nx-4)*nz+iz] - bnd_coef[ix*nc2+4] * in[(nx-5)*nz+iz] - bnd_coef[ix*nc2+5] * in[(nx-6)*nz+iz]) / d; 
    }

    // apply the operator in the interior
    #pragma omp parallel for
    for (int ix=ixminb; ix<ixmaxb; ix++){
        Dx_ispc(add, in, out, 1.0/d,  ix, nz, izmin,  izmax);
    }
}

void mult_Dz(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax, const data_t * par, data_t a){
    
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    
    int nc1=4, nc2=6;
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);

    #pragma omp parallel for
    for (int ix = ixmin; ix < ixmax; ix++){
        int i1=ix*nz;

        // apply the operator near the top boundary if included
        for (int iz=izmin; iz<nc1; iz++){
             out[i1+iz] = add* out[i1+iz] + a/d * (bnd_coef[iz*nc2] * in[i1]*par[i1] + bnd_coef[iz*nc2+1] * in[i1+1]*par[i1+1] + bnd_coef[iz*nc2+2] * in[i1+2]*par[i1+2] + bnd_coef[iz*nc2+3] * in[i1+3]*par[i1+3] + bnd_coef[iz*nc2+4] * in[i1+4]*par[i1+4] + bnd_coef[iz*nc2+5] * in[i1+5]*par[i1+5]);
        }

        // apply the operator near the bottom boundary if included
        for (int iz=nz-izmax; iz<nc1; iz++){
             out[i1+nz-1-iz] = add* out[i1+nz-1-iz] + a/d * (-bnd_coef[iz*nc2] * in[i1+nz-1]*par[i1+nz-1] - bnd_coef[iz*nc2+1] * in[i1+nz-2]*par[i1+nz-2] - bnd_coef[iz*nc2+2] * in[i1+nz-3]*par[i1+nz-3] - bnd_coef[iz*nc2+3] * in[i1+nz-4]*par[i1+nz-4] - bnd_coef[iz*nc2+4] * in[i1+nz-5]*par[i1+nz-5] - bnd_coef[iz*nc2+5] * in[i1+nz-6]*par[i1+nz-6]);
        }

        // apply the operator in the interior
        mult_Dz_ispc(add, in, out,  a/d, i1,  izminb,  izmaxb, par);
    }
}

void mult_Dx(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax, const data_t * par, data_t a){
    
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};

    int nc1=4, nc2=6;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);

    // apply the operator near the left boundary if included
    #pragma omp parallel for
    for (int ix=ixmin; ix<nc1; ix++){
        int i1=ix*nz;
        for (int iz = izmin; iz < izmax; iz++) out[i1+iz] = add*out[i1+iz] + a/d * (bnd_coef[ix*nc2] * in[iz]*par[iz] + bnd_coef[ix*nc2+1] * in[nz+iz]*par[nz+iz] + bnd_coef[ix*nc2+2] * in[2*nz+iz]*par[2*nz+iz] + bnd_coef[ix*nc2+3] * in[3*nz+iz]*par[3*nz+iz] + bnd_coef[ix*nc2+4] * in[4*nz+iz]*par[4*nz+iz] + bnd_coef[ix*nc2+5] * in[5*nz+iz]*par[5*nz+iz]); 
    }

    // apply the operator near the right boundary if included
    #pragma omp parallel for
    for (int ix=nx-ixmax; ix<nc1; ix++){
        int i1=(nx-1-ix);
        for (int iz = izmin; iz < izmax; iz++) out[i1*nz+iz] = add*out[i1*nz+iz] + a/d * (-bnd_coef[ix*nc2] * in[(nx-1)*nz+iz]*par[(nx-1)*nz+iz] - bnd_coef[ix*nc2+1] * in[(nx-2)*nz+iz]*par[(nx-2)*nz+iz] - bnd_coef[ix*nc2+2] * in[(nx-3)*nz+iz]*par[(nx-3)*nz+iz] - bnd_coef[ix*nc2+3] * in[(nx-4)*nz+iz]*par[(nx-4)*nz+iz] - bnd_coef[ix*nc2+4] * in[(nx-5)*nz+iz]*par[(nx-5)*nz+iz] - bnd_coef[ix*nc2+5] * in[(nx-6)*nz+iz]*par[(nx-6)*nz+iz]); 
    }

    // apply the operator in the interior
    #pragma omp parallel for
    for (int ix=ixminb; ix<ixmaxb; ix++){
        mult_Dx_ispc(add, in, out, a/d, ix, nz, izmin,  izmax, par);
    }
}

void esat_scale_boundaries(data_t** in, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, int izmin, int izmax, const data_t** par, data_t dt, bool top, bool bottom, bool left, bool right){

    // par must at least contain lambda, mu, rho in this order
    data_t h0 = 17.0/48;
    data_t scalerx, scalerz;
    data_t sc1x=0, sc1z=0, sc2x=0, sc2z=0, sc3x=0, sc3z=0, sc4x=0, sc4z=0;

    // scale the top boundary without the corners
    if (top){
        sc1x = sqrt(par[1][0]/par[2][0])*dt / (2  * dz * h0); // top left
        sc1z = sc1x*sqrt((par[0][0]+2*par[1][0])/par[1][0]);
        sc2x = sqrt(par[1][(nx-1)*nz]/par[2][(nx-1)*nz])*dt / (2 * dz * h0); // top right
        sc2z = sc1x*sqrt((par[0][(nx-1)*nz]+2*par[1][(nx-1)*nz])/par[1][(nx-1)*nz]);
        if (izmin==0){
            for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                scalerx = sqrt(par[1][ix*nz]/par[2][ix*nz])*dt / (2 * dz * h0);
                scalerz = scalerx*sqrt((par[0][ix*nz]+2*par[1][ix*nz])/par[1][ix*nz]);
                in[0][ix*nz] = in[0][ix*nz] / (1+scalerx);
                in[1][ix*nz] = in[1][ix*nz] / (1+scalerz);
            }
        }
    }

    // scale the left boundary
    if (left){
        sc1x += sqrt((par[0][0]+2*par[1][0])/par[2][0])*dt / (2 * dx * h0); // top left
        sc1z += sqrt(par[1][0]/par[2][0])*dt / (2 * dx * h0);
        in[0][0] = in[0][0] / (1+sc1x);
        in[1][0] = in[1][0] / (1+sc1z);
        sc1x = sqrt((par[0][nz-1]+2*par[1][nz-1])/par[2][nz-1])*dt / (2 * dx * h0); // bottom left
        sc1z = sc1x*sqrt(par[1][nz-1]/(par[0][nz-1]+2*par[1][nz-1]));
        if (ixmin==0){
            for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                scalerx = sqrt((par[0][iz]+2*par[1][iz])/par[2][iz])*dt / (2 * dx * h0);
                scalerz = scalerx*sqrt(par[1][iz]/(par[0][iz]+2*par[1][iz]));
                in[0][iz] = in[0][iz] / (1+scalerx);
                in[1][iz] = in[1][iz] / (1+scalerz);
            }
        }
    }
    else{
        if (izmin==0 && ixmin==0){
            in[0][0] = in[0][0] / (1+sc1x); // top left
            in[1][0] = in[1][0] / (1+sc1z);
        }
        sc1x = 0; // bottom left
        sc1z = 0;
    }

    // scale the bottom boundary
    if (bottom){
        sc1x += sqrt(par[1][nz-1]/par[2][nz-1])*dt / (2 * dz * h0); // bottom left
        sc1z += sqrt((par[0][nz-1]+2*par[1][nz-1])/par[2][nz-1])*dt / (2 * dz * h0);
        in[0][nz-1] = in[0][nz-1] / (1+sc1x);
        in[1][nz-1] = in[1][nz-1] / (1+sc1z);
        sc1x = sqrt(par[1][(nx-1)*nz+nz-1]/par[2][(nx-1)*nz+nz-1])*dt / (2 * dz * h0); // bottom right
        sc1z = sc1x*sqrt((par[0][(nx-1)*nz+nz-1]+2*par[1][(nx-1)*nz+nz-1])/par[1][(nx-1)*nz+nz-1]);
        if (izmax==nz){
            for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                scalerx = sqrt(par[1][ix*nz+nz-1]/par[2][ix*nz+nz-1])*dt / (2 * dz * h0);
                scalerz = scalerx*sqrt((par[0][ix*nz+nz-1]+2*par[1][ix*nz+nz-1])/par[1][ix*nz+nz-1]);
                in[0][ix*nz+nz-1] = in[0][ix*nz+nz-1] / (1+scalerx);
                in[1][ix*nz+nz-1] = in[1][ix*nz+nz-1] / (1+scalerz);
            }
        }
    }
    else {
        if (izmax==nz && ixmin==0){
            in[0][nz-1] = in[0][nz-1] / (1+sc1x); // bottom left
            in[1][nz-1] = in[1][nz-1] / (1+sc1z);
        }
        sc1x = 0; // bottom right
        sc1z = 0;
    }

    // scale the right boundary
    if (right){
        sc1x += sqrt((par[0][(nx-1)*nz+nz-1]+2*par[1][(nx-1)*nz+nz-1])/par[2][(nx-1)*nz+nz-1])*dt / (2 * dx * h0); // bottom right
        sc1z += sqrt(par[1][(nx-1)*nz+nz-1]/par[2][(nx-1)*nz+nz-1])*dt / (2 * dx * h0);
        in[0][(nx-1)*nz+nz-1] = in[0][(nx-1)*nz+nz-1] / (1+sc1x);
        in[1][(nx-1)*nz+nz-1] = in[1][(nx-1)*nz+nz-1] / (1+sc1z);
        if (ixmax==nx){
            for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                scalerx = sqrt((par[0][(nx-1)*nz+iz]+2*par[1][(nx-1)*nz+iz])/par[2][(nx-1)*nz+iz])*dt / (2 * dx * h0);
                scalerz = scalerx*sqrt(par[1][(nx-1)*nz+iz]/(par[0][(nx-1)*nz+iz]+2*par[1][(nx-1)*nz+iz]));
                in[0][(nx-1)*nz+iz] = in[0][(nx-1)*nz+iz] / (1+scalerx);
                in[1][(nx-1)*nz+iz] = in[1][(nx-1)*nz+iz] / (1+scalerz);
            }
        }
        sc2x += sqrt((par[0][(nx-1)*nz]+2*par[1][(nx-1)*nz])/par[2][(nx-1)*nz])*dt / (2 * dx *h0); // top right
        sc2z += sqrt(par[1][(nx-1)*nz]/par[2][(nx-1)*nz])*dt / (2 * dx * h0);
        if (izmin==0 && ixmax==nx){
            in[0][(nx-1)*nz] = in[0][(nx-1)*nz] / (1+sc2x);
            in[1][(nx-1)*nz] = in[1][(nx-1)*nz] / (1+sc2z);
        }
    }
    else {
        if (izmax==nz && ixmax==nx){
            in[0][(nx-1)*nz+nz-1] = in[0][(nx-1)*nz+nz-1] / (1+sc1x); // bottom right
            in[1][(nx-1)*nz+nz-1] = in[1][(nx-1)*nz+nz-1] / (1+sc1z);
        }
        if (izmin==0 && ixmax==nx){
            in[0][(nx-1)*nz] = in[0][(nx-1)*nz] / (1+sc2x); // top right
            in[1][(nx-1)*nz] = in[1][(nx-1)*nz] / (1+sc2z);
        }
    }
}

void vtisat_scale_boundaries(data_t** in, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, int izmin, int izmax, const data_t** par, data_t dt, bool top, bool bottom, bool left, bool right){

    // par must at least contain lambda, mu, rho, c13, eps in this order
    data_t h0 = 17.0/48;
    data_t scalerx, scalerz;
    data_t sc1x=0, sc1z=0, sc2x=0, sc2z=0, sc3x=0, sc3z=0, sc4x=0, sc4z=0;

    // scale the top boundary without the corners
    if (top){
        sc1x = sqrt(par[1][0]/par[2][0])*dt / (2  * dz * h0); // top left
        sc1z = sc1x*sqrt((par[0][0]+2*par[1][0])/par[1][0]);
        sc2x = sqrt(par[1][(nx-1)*nz]/par[2][(nx-1)*nz])*dt / (2 * dz * h0); // top right
        sc2z = sc1x*sqrt((par[0][(nx-1)*nz]+2*par[1][(nx-1)*nz])/par[1][(nx-1)*nz]);
        if (izmin==0){
            for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                scalerx = sqrt(par[1][ix*nz]/par[2][ix*nz])*dt / (2 * dz * h0);
                scalerz = scalerx*sqrt((par[0][ix*nz]+2*par[1][ix*nz])/par[1][ix*nz]);
                in[0][ix*nz] = in[0][ix*nz] / (1+scalerx);
                in[1][ix*nz] = in[1][ix*nz] / (1+scalerz);
            }
        }
    }

    // scale the left boundary
    if (left){
        sc1x += sqrt((1+2*par[4][0])*(par[0][0]+2*par[1][0])/par[2][0])*dt / (2 * dx * h0); // top left
        sc1z += sqrt(par[1][0]/par[2][0])*dt / (2 * dx * h0);
        in[0][0] = in[0][0] / (1+sc1x);
        in[1][0] = in[1][0] / (1+sc1z);
        sc1x = sqrt((1+2*par[4][nz-1])*(par[0][nz-1]+2*par[1][nz-1])/par[2][nz-1])*dt / (2 * dx * h0); // bottom left
        sc1z = sc1x*sqrt(par[1][nz-1]/((1+2*par[4][nz-1])*(par[0][nz-1]+2*par[1][nz-1])));
        if (ixmin==0){
            for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                scalerx = sqrt((1+2*par[4][iz])*(par[0][iz]+2*par[1][iz])/par[2][iz])*dt / (2 * dx * h0);
                scalerz = scalerx*sqrt(par[1][iz]/((1+2*par[4][iz])*(par[0][iz]+2*par[1][iz])));
                in[0][iz] = in[0][iz] / (1+scalerx);
                in[1][iz] = in[1][iz] / (1+scalerz);
            }
        }
    }
    else{
        if (izmin==0 && ixmin==0){
            in[0][0] = in[0][0] / (1+sc1x); // top left
            in[1][0] = in[1][0] / (1+sc1z);
        }
        sc1x = 0; // bottom left
        sc1z = 0;
    }

    // scale the bottom boundary
    if (bottom){
        sc1x += sqrt(par[1][nz-1]/par[2][nz-1])*dt / (2 * dz * h0); // bottom left
        sc1z += sqrt((par[0][nz-1]+2*par[1][nz-1])/par[2][nz-1])*dt / (2 * dz * h0);
        in[0][nz-1] = in[0][nz-1] / (1+sc1x);
        in[1][nz-1] = in[1][nz-1] / (1+sc1z);
        sc1x = sqrt(par[1][(nx-1)*nz+nz-1]/par[2][(nx-1)*nz+nz-1])*dt / (2 * dz * h0); // bottom right
        sc1z = sc1x*sqrt((par[0][(nx-1)*nz+nz-1]+2*par[1][(nx-1)*nz+nz-1])/par[1][(nx-1)*nz+nz-1]);
        if (izmax==nz){
            for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                scalerx = sqrt(par[1][ix*nz+nz-1]/par[2][ix*nz+nz-1])*dt / (2 * dz * h0);
                scalerz = scalerx*sqrt((par[0][ix*nz+nz-1]+2*par[1][ix*nz+nz-1])/par[1][ix*nz+nz-1]);
                in[0][ix*nz+nz-1] = in[0][ix*nz+nz-1] / (1+scalerx);
                in[1][ix*nz+nz-1] = in[1][ix*nz+nz-1] / (1+scalerz);
            }
        }
    }
    else {
        if (izmax==nz && ixmin==0){
            in[0][nz-1] = in[0][nz-1] / (1+sc1x); // bottom left
            in[1][nz-1] = in[1][nz-1] / (1+sc1z);
        }
        sc1x = 0; // bottom right
        sc1z = 0;
    }

    // scale the right boundary
    if (right){
        sc1x += sqrt((1+2*par[4][(nx-1)*nz+nz-1])*(par[0][(nx-1)*nz+nz-1]+2*par[1][(nx-1)*nz+nz-1])/par[2][(nx-1)*nz+nz-1])*dt / (2 * dx * h0); // bottom right
        sc1z += sqrt(par[1][(nx-1)*nz+nz-1]/par[2][(nx-1)*nz+nz-1])*dt / (2 * dx * h0);
        in[0][(nx-1)*nz+nz-1] = in[0][(nx-1)*nz+nz-1] / (1+sc1x);
        in[1][(nx-1)*nz+nz-1] = in[1][(nx-1)*nz+nz-1] / (1+sc1z);
        if (ixmax==nx){
            for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                scalerx = sqrt((1+2*par[4][(nx-1)*nz+iz])*(par[0][(nx-1)*nz+iz]+2*par[1][(nx-1)*nz+iz])/par[2][(nx-1)*nz+iz])*dt / (2 * dx * h0);
                scalerz = scalerx*sqrt(par[1][(nx-1)*nz+iz]/((1+2*par[4][(nx-1)*nz+iz])*(par[0][(nx-1)*nz+iz]+2*par[1][(nx-1)*nz+iz])));
                in[0][(nx-1)*nz+iz] = in[0][(nx-1)*nz+iz] / (1+scalerx);
                in[1][(nx-1)*nz+iz] = in[1][(nx-1)*nz+iz] / (1+scalerz);
            }
        }
        sc2x += sqrt((1+2*par[4][(nx-1)*nz])*(par[0][(nx-1)*nz]+2*par[1][(nx-1)*nz])/par[2][(nx-1)*nz])*dt / (2 * dx *h0); // top right
        sc2z += sqrt(par[1][(nx-1)*nz]/par[2][(nx-1)*nz])*dt / (2 * dx * h0);
        if (izmin==0 && ixmax==nx){
            in[0][(nx-1)*nz] = in[0][(nx-1)*nz] / (1+sc2x);
            in[1][(nx-1)*nz] = in[1][(nx-1)*nz] / (1+sc2z);
        }
    }
    else {
        if (izmax==nz && ixmax==nx){
            in[0][(nx-1)*nz+nz-1] = in[0][(nx-1)*nz+nz-1] / (1+sc1x); // bottom right
            in[1][(nx-1)*nz+nz-1] = in[1][(nx-1)*nz+nz-1] / (1+sc1z);
        }
        if (izmin==0 && ixmax==nx){
            in[0][(nx-1)*nz] = in[0][(nx-1)*nz] / (1+sc2x); // top right
            in[1][(nx-1)*nz] = in[1][(nx-1)*nz] / (1+sc2z);
        }
    }
}

void taperz(data_t* in, int nx, int nz, int ixmin, int ixmax, int istart, int iend, data_t a){

    #pragma omp parallel for
    for (int ix = ixmin; ix < ixmax; ix++){
        int i1=ix*nz;
        taperz_ispc(in, i1, istart, iend, a);
    }
}
void taperx(data_t* in, int nx, int nz, int izmin, int izmax, int istart, int iend, data_t a){

    if (iend>istart){
        #pragma omp parallel for
        for (int ix = istart; ix < iend; ix++){
            int i1=ix*nz;
            data_t val = cos(a*0.5*M_PI*(ix-istart)/(iend-istart));
            taperx_ispc(in, i1, nz, izmin, izmax, val);
        }
    }
    else{
        #pragma omp parallel for
        for (int ix = iend; ix < istart; ix++){
            int i1=ix*nz;
            data_t val = cos(a*0.5*M_PI*(ix+1-istart)/(iend-istart));
            taperx_ispc(in, i1, nz, izmin, izmax, val);
        }
    }
}

void asat_dirichlet_top(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, const data_t ** par, data_t a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t hcoef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};
    int nc1=4;
    data_t dz2=dz*dz;

    // SAT = + H-1 (-S'.(f2*in_0)) - H-1(f2*in_0/(h.a))_0  f2=reciprocal of density
    #pragma omp parallel for
    for (int ix=ixmin; ix<ixmax; ix++){
        for (int iz=0; iz<nc1; iz++){
            out[ix*nz+iz] = add*out[ix*nz+iz] + 1.0/(hcoef[iz]*dz2)*scoef[iz]*par[1][ix*nz]*in[0][ix*nz];
        }
        out[ix*nz] -= par[1][ix*nz]*in[0][ix*nz]/(a*hcoef[0]*dz2);
    }
}

void asat_dirichlet_bottom(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, const data_t ** par, data_t a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t hcoef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};
    int nc1=4;
    data_t dz2=dz*dz;

    // SAT = + H-1 (S'.(f2*in_0)) - H-1(f2*in_0/(h.a))_0  f2=reciprocal of density
    #pragma omp parallel for
    for (int ix=ixmin; ix<ixmax; ix++){
        for (int iz=0; iz<nc1; iz++){
            out[ix*nz+nz-1-iz] = add*out[ix*nz+nz-1-iz] + 1.0/(hcoef[iz]*dz2)*scoef[iz]*par[1][ix*nz+nz-1]*in[0][ix*nz+nz-1];
        }
        out[ix*nz+nz-1] -= par[1][ix*nz+nz-1]*in[0][ix*nz+nz-1]/(a*hcoef[0]*dz2);
    }
}

void asat_dirichlet_left(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int izmin, int izmax, const data_t ** par, data_t a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t hcoef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};
    int nc1=4;
    data_t dx2=dx*dx;

    // SAT = + H-1 (-S'.(f2*in_0)) - H-1(f2*in_0/(h.a))_0  f2=reciprocal of density
    #pragma omp parallel for
    for (int iz=izmin; iz<izmax; iz++){
        for (int ix=0; ix<nc1; ix++){
            out[ix*nz+iz] = add*out[ix*nz+iz] + 1.0/(hcoef[ix]*dx2)*scoef[ix]*par[1][iz]*in[0][iz];
        }
    }
    #pragma omp parallel for
    for (int iz=izmin; iz<izmax; iz++) out[iz] -= par[1][iz]*in[0][iz]/(a*hcoef[0]*dx2);
}

void asat_dirichlet_right(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int izmin, int izmax, const data_t ** par, data_t a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t hcoef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};
    int nc1=4;
    data_t dx2=dx*dx;

    // SAT = + H-1 (S'.(f2*in_0)) - H-1(f2*in_0/(h.a))_0  f2=reciprocal of density
    #pragma omp parallel for
    for (int iz=izmin; iz<izmax; iz++){
        for (int ix=0; ix<nc1; ix++){
            out[(nx-1-ix)*nz+iz] = add*out[(nx-1-ix)*nz+iz] + 1.0/(hcoef[ix]*dx2)*scoef[ix]*par[1][(nx-1)*nz+iz]*in[0][(nx-1)*nz+iz];
        }
    }
    #pragma omp parallel for
    for (int iz=izmin; iz<izmax; iz++) out[(nx-1)*nz+iz] -= par[1][(nx-1)*nz+iz]*in[0][(nx-1)*nz+iz]/(a*hcoef[0]*dx2);
}

void asat_neumann_top(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, const data_t ** par, data_t a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4;
    data_t sum=0;

    // SAT = - H-1 (-f2.S.in )_0  f2=reciprocal of density
    #pragma omp parallel for private(sum)
    for (int ix=ixmin; ix<ixmax; ix++){
        sum=0;
        for (int iz=0; iz<nc1; iz++){
            sum += scoef[iz] * in[0][ix*nz+iz];
        }
        out[ix*nz] = add*out[ix*nz] - 1.0/(h0*dz) * par[1][ix*nz]*sum/dz;
    }
}

void asat_neumann_bottom(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, const data_t ** par, data_t a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4;
    data_t sum=0;

    // SAT = - H-1 (f2.S.in )_0  f2=reciprocal of density
    #pragma omp parallel for private(sum)
    for (int ix=ixmin; ix<ixmax; ix++){
        sum=0;
        for (int iz=0; iz<nc1; iz++){
            sum += scoef[iz] * in[0][ix*nz+nz-1-iz];
        }
        out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - 1.0/(h0*dz) * par[1][ix*nz+nz-1]*sum/dz;
    }
}

void asat_neumann_left(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int izmin, int izmax, const data_t ** par, data_t a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4;
    data_t sum=0;

    // SAT = - H-1 (-f2.S.in)_0  f2=reciprocal of density
    #pragma omp parallel for private(sum)
    for (int iz=izmin; iz<izmax; iz++){
        sum=0;
        for (int ix=0; ix<nc1; ix++){
            sum += scoef[ix] * in[0][ix*nz+iz];
        }
        out[iz] = add*out[iz] - 1.0/(h0*dx) * par[1][iz]*sum/dx;
    }
}

void asat_neumann_right(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int izmin, int izmax, const data_t ** par, data_t a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4;
    data_t sum=0;

    // SAT = - H-1 (f2.S.in)_0  f2=reciprocal of density
    #pragma omp parallel for private(sum)
    for (int iz=izmin; iz<izmax; iz++){
        sum=0;
        for (int ix=0; ix<nc1; ix++){
            sum += scoef[ix] * in[0][(nx-1-ix)*nz+iz];
        }
        out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - 1.0/(h0*dx) * par[1][(nx-1)*nz+iz]*sum/dx;
    }
}

void asat_absorbing_top(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int ixmin, int ixmax, const data_t ** par, data_t a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4;
    data_t sum=0;

    // SAT = - H-1 (-f2.S.in - f3.in/(2dt) )_0  f2=reciprocal of density ; f3=1/(rho.v)=1/sqrt(rho.K)
    #pragma omp parallel for private(sum)
    for (int ix=ixmin; ix<ixmax; ix++){
        sum=0;
        for (int iz=0; iz<nc1; iz++){
            sum += scoef[iz] * in[0][ix*nz+iz];
        }
        out[ix*nz] = add*out[ix*nz] - 1.0/(h0*dz) * ( par[1][ix*nz]*sum/dz - in[1][ix*nz]/(2*dt*sqrt(par[0][ix*nz]/par[1][ix*nz])));
    }
}

void asat_absorbing_bottom(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int ixmin, int ixmax, const data_t ** par, data_t a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4;
    data_t sum=0;

    // SAT = - H-1 (f2.S.in - f3.in/(2dt) )_0  f2=reciprocal of density ; f3=1/(rho.v)=1/sqrt(rho.K)
    #pragma omp parallel for private(sum)
    for (int ix=ixmin; ix<ixmax; ix++){
        sum=0;
        for (int iz=0; iz<nc1; iz++){
            sum += scoef[iz] * in[0][ix*nz+nz-1-iz];
        }
        out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - 1.0/(h0*dz) * ( par[1][ix*nz+nz-1]*sum/dz - in[1][ix*nz+nz-1]/(2*dt*sqrt(par[0][ix*nz+nz-1]/par[1][ix*nz+nz-1])));
    }
}

void asat_absorbing_left(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int izmin, int izmax, const data_t ** par, data_t a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4;
    data_t sum=0;

    // SAT = - H-1 (-f2.S.in - f3.in/(2dt) )_0  f2=reciprocal of density ; f3=1/(rho.v)=1/sqrt(rho.K)
    #pragma omp parallel for private(sum)
    for (int iz=izmin; iz<izmax; iz++){
        sum=0;
        for (int ix=0; ix<nc1; ix++){
            sum += scoef[ix] * in[0][ix*nz+iz];
        }
        out[iz] = add*out[iz] - 1.0/(h0*dx) * ( par[1][iz]*sum/dx - in[1][iz]/(2*dt*sqrt(par[0][iz]/par[1][iz])));
    }
}

void asat_absorbing_right(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int izmin, int izmax, const data_t ** par, data_t a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4;
    data_t sum=0;

    // SAT = - H-1 (f2.S.in - f3.in/(2dt) )_0  f2=reciprocal of density ; f3=1/(rho.v)=1/sqrt(rho.K)
    #pragma omp parallel for private(sum)
    for (int iz=izmin; iz<izmax; iz++){
        sum=0;
        for (int ix=0; ix<nc1; ix++){
            sum += scoef[ix] * in[0][(nx-1-ix)*nz+iz];
        }
        out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - 1.0/(h0*dx) * ( par[1][(nx-1)*nz+iz]*sum/dx - in[1][(nx-1)*nz+iz]/(2*dt*sqrt(par[0][(nx-1)*nz+iz]/par[1][(nx-1)*nz+iz])));
    }
}

void asat_scale_boundaries(data_t** in, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, int izmin, int izmax, const data_t** par, data_t dt, bool top, bool bottom, bool left, bool right){

    // par must at least contain K, 1/rho in this order
    data_t h0 = 17.0/48;
    data_t scaler;
    data_t sc1=0, sc2=0, sc3=0, sc4=0;

    // scale the top boundary without the corners
    if (top){
        sc1 = sqrt(par[0][0]*par[1][0])*dt / (2  * dz * h0); // top left
        sc2 = sqrt(par[0][(nx-1)*nz]*par[1][(nx-1)*nz])*dt / (2 * dz * h0); // top right
        if (izmin==0){
            for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                scaler = sqrt(par[0][ix*nz]*par[1][ix*nz])*dt / (2 * dz * h0);
                in[0][ix*nz] = in[0][ix*nz] / (1+scaler);
            }
        }
    }

    // scale the left boundary
    if (left){
        sc1 += sqrt(par[0][0]*par[1][0])*dt / (2 * dx * h0); // top left
        in[0][0] = in[0][0] / (1+sc1);
        sc1 = sqrt(par[0][nz-1]*par[1][nz-1])*dt / (2 * dx * h0); // bottom left
        if (ixmin==0){
            for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                scaler = sqrt(par[0][iz]*par[1][iz])*dt / (2 * dx * h0);
                in[0][iz] = in[0][iz] / (1+scaler);
            }
        }
    }
    else{
        if (izmin==0 && ixmin==0){
            in[0][0] = in[0][0] / (1+sc1); // top left
        }
        sc1 = 0; // bottom left
    }

    // scale the bottom boundary
    if (bottom){
        sc1 += sqrt(par[0][nz-1]*par[1][nz-1])*dt / (2 * dz * h0); // bottom left
        in[0][nz-1] = in[0][nz-1] / (1+sc1);
        sc1 = sqrt(par[0][(nx-1)*nz+nz-1]*par[1][(nx-1)*nz+nz-1])*dt / (2 * dz * h0); // bottom right
        if (izmax==nz){
            for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                scaler = sqrt(par[0][ix*nz+nz-1]*par[1][ix*nz+nz-1])*dt / (2 * dz * h0);
                in[0][ix*nz+nz-1] = in[0][ix*nz+nz-1] / (1+scaler);
            }
        }
    }
    else {
        if (izmax==nz && ixmin==0){
            in[0][nz-1] = in[0][nz-1] / (1+sc1); // bottom left
        }
        sc1 = 0; // bottom right
    }

    // scale the right boundary
    if (right){
        sc1 += sqrt(par[0][(nx-1)*nz+nz-1]*par[1][(nx-1)*nz+nz-1])*dt / (2 * dx * h0); // bottom right
        in[0][(nx-1)*nz+nz-1] = in[0][(nx-1)*nz+nz-1] / (1+sc1);
        if (ixmax==nx){
            for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                scaler = sqrt(par[0][(nx-1)*nz+iz]*par[1][(nx-1)*nz+iz])*dt / (2 * dx * h0);
                in[0][(nx-1)*nz+iz] = in[0][(nx-1)*nz+iz] / (1+scaler);
            }
        }
        sc2 += sqrt(par[0][(nx-1)*nz]*par[1][(nx-1)*nz])*dt / (2 * dx *h0); // top right
        if (izmin==0 && ixmax==nx){
            in[0][(nx-1)*nz] = in[0][(nx-1)*nz] / (1+sc2);
        }
    }
    else {
        if (izmax==nz && ixmax==nx){
            in[0][(nx-1)*nz+nz-1] = in[0][(nx-1)*nz+nz-1] / (1+sc1); // bottom right
        }
        if (izmin==0 && ixmax==nx){
            in[0][(nx-1)*nz] = in[0][(nx-1)*nz] / (1+sc2); // top right
        }
    }
}


void asat_neumann_absorbing_top(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int ixmin, int ixmax, const data_t ** par, const data_t * a)
{

    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t h0 = 17.0/48;
    int nc1=4;
    data_t sum=0;

    // SAT = - H-1 (-f2.S.in - f3.in/(2dt) )_0  f2=reciprocal of density ; f3=a/(rho.v)=a/sqrt(rho.K)  ;  'a' is like a transmission coefficient (0<=a<=1)
    #pragma omp parallel for private(sum)
    for (int ix=ixmin; ix<ixmax; ix++){
        sum=0;
        for (int iz=0; iz<nc1; iz++){
            sum += scoef[iz] * in[0][ix*nz+iz];
        }
        out[ix*nz] = add*out[ix*nz] - 1.0/(h0*dz) * ( par[1][ix*nz]*sum/dz - in[1][ix*nz]*a[ix]/(2*dt*sqrt(par[0][ix*nz]/par[1][ix*nz])));
    }
}

void asat_neumann_dirichlet_top(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, const data_t ** par, data_t a, const data_t * m)
{
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3};
    data_t hcoef[4] = {17.0/48, 59.0/48, 43.0/48, 49.0/48};
    int nc1=4;
    data_t dz2=dz*dz;
    data_t sum=0;

    // SAT = (1-m).[ -H-1 (-f2.S.in )_0 ] + m.[ H-1 (-S'.(f2.in_0)) - H-1(f2.in_0/(h.a))_0 ] f2=reciprocal of density ; 'm' is like a mask (0,1) to flip between the two SATs
    #pragma omp parallel for private(sum)
    for (int ix=ixmin; ix<ixmax; ix++){
        if (m[ix]==0)
        {
            sum=0;
            for (int iz=0; iz<nc1; iz++){
                sum += scoef[iz] * in[0][ix*nz+iz];
            }
            out[ix*nz] = add*out[ix*nz] - 1.0/(hcoef[0]*dz) * par[1][ix*nz]*sum/dz;
        }
        else
        {
            for (int iz=0; iz<nc1; iz++){
                out[ix*nz+iz] = add*out[ix*nz+iz] + 1.0/(hcoef[iz]*dz2)*scoef[iz]*par[1][ix*nz]*in[0][ix*nz];
            }
            out[ix*nz] -= par[1][ix*nz]*in[0][ix*nz]/(a*hcoef[0]*dz2);
        }
    }
}

void asat_scale_boundaries_bis(data_t** in, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, int izmin, int izmax, const data_t** par, const data_t * a, data_t dt, bool top, bool bottom, bool left, bool right){

    // par must at least contain K, 1/rho in this order
    data_t h0 = 17.0/48;
    data_t scaler;
    data_t sc1=0, sc2=0, sc3=0, sc4=0;

    // scale the top boundary without the corners
    if (top){
        sc1 = a[0]*sqrt(par[0][0]*par[1][0])*dt / (2  * dz * h0); // top left
        sc2 = a[nx-1]*sqrt(par[0][(nx-1)*nz]*par[1][(nx-1)*nz])*dt / (2 * dz * h0); // top right
        if (izmin==0){
            for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                scaler = a[ix]*sqrt(par[0][ix*nz]*par[1][ix*nz])*dt / (2 * dz * h0);
                in[0][ix*nz] = in[0][ix*nz] / (1+scaler);
            }
        }
    }

    // scale the left boundary
    if (left){
        sc1 += sqrt(par[0][0]*par[1][0])*dt / (2 * dx * h0); // top left
        in[0][0] = in[0][0] / (1+sc1);
        sc1 = sqrt(par[0][nz-1]*par[1][nz-1])*dt / (2 * dx * h0); // bottom left
        if (ixmin==0){
            for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                scaler = sqrt(par[0][iz]*par[1][iz])*dt / (2 * dx * h0);
                in[0][iz] = in[0][iz] / (1+scaler);
            }
        }
    }
    else{
        if (izmin==0 && ixmin==0){
            in[0][0] = in[0][0] / (1+sc1); // top left
        }
        sc1 = 0; // bottom left
    }

    // scale the bottom boundary
    if (bottom){
        sc1 += sqrt(par[0][nz-1]*par[1][nz-1])*dt / (2 * dz * h0); // bottom left
        in[0][nz-1] = in[0][nz-1] / (1+sc1);
        sc1 = sqrt(par[0][(nx-1)*nz+nz-1]*par[1][(nx-1)*nz+nz-1])*dt / (2 * dz * h0); // bottom right
        if (izmax==nz){
            for (int ix = std::max(1,ixmin); ix < std::min(nx-1,ixmax); ix++){
                scaler = sqrt(par[0][ix*nz+nz-1]*par[1][ix*nz+nz-1])*dt / (2 * dz * h0);
                in[0][ix*nz+nz-1] = in[0][ix*nz+nz-1] / (1+scaler);
            }
        }
    }
    else {
        if (izmax==nz && ixmin==0){
            in[0][nz-1] = in[0][nz-1] / (1+sc1); // bottom left
        }
        sc1 = 0; // bottom right
    }

    // scale the right boundary
    if (right){
        sc1 += sqrt(par[0][(nx-1)*nz+nz-1]*par[1][(nx-1)*nz+nz-1])*dt / (2 * dx * h0); // bottom right
        in[0][(nx-1)*nz+nz-1] = in[0][(nx-1)*nz+nz-1] / (1+sc1);
        if (ixmax==nx){
            for (int iz = std::max(1,izmin); iz < std::min(nz-1,izmax); iz++){
                scaler = sqrt(par[0][(nx-1)*nz+iz]*par[1][(nx-1)*nz+iz])*dt / (2 * dx * h0);
                in[0][(nx-1)*nz+iz] = in[0][(nx-1)*nz+iz] / (1+scaler);
            }
        }
        sc2 += sqrt(par[0][(nx-1)*nz]*par[1][(nx-1)*nz])*dt / (2 * dx *h0); // top right
        if (izmin==0 && ixmax==nx){
            in[0][(nx-1)*nz] = in[0][(nx-1)*nz] / (1+sc2);
        }
    }
    else {
        if (izmax==nz && ixmax==nx){
            in[0][(nx-1)*nz+nz-1] = in[0][(nx-1)*nz+nz-1] / (1+sc1); // bottom right
        }
        if (izmin==0 && ixmax==nx){
            in[0][(nx-1)*nz] = in[0][(nx-1)*nz] / (1+sc2); // top right
        }
    }
}

// variable coefficients expressions for 2nd derivative SBP operators
static inline data_t expr1(const data_t ** par, int i){return par[0][i];} // e.g. = lambda
static inline data_t expr2(const data_t ** par, int i){return par[1][i];} // e.g. = mu
static inline data_t expr3(const data_t ** par, int i){return par[2][i];} // e.g. = rho
static inline data_t expr4a(const data_t ** par, int i){return par[0][i] + 2*par[1][i];} // e.g. = lambda + 2 mu
static inline data_t expr4b(const data_t ** par, int i){return 2*par[1][i];} // e.g. = 2 mu

void pml_top(const data_t * in, data_t * out, // input and output wavefields
            const data_t * u_x, const data_t * u_z, // first derivative wavefields
            data_t * &ac, data_t * &an, data_t * s2, // auxilary arrays to store the split fields s1, t, s3, s4, and s2
            const data_t * model, // elastic parameters lamda, mu, rho
            const data_t * g, const data_t * gprime, // stretching function and its derivative
            int nx, int nz, int l, int ixmin, int ixmax, data_t dx, data_t dz, data_t dt, // sizes and sampling
            int version)
{
    data_t (*pin) [nx*nz] = (data_t (*)[nx*nz]) in;
    data_t (*pout) [nx*nz] = (data_t (*)[nx*nz]) out;
    data_t (*pu_x) [nx*nz] = (data_t (*)[nx*nz]) u_x;
    data_t (*pu_z) [nx*nz] = (data_t (*)[nx*nz]) u_z;
    data_t (*pac) [2*nx*l] = (data_t (*)[2*nx*l]) ac;
    data_t (*pan) [2*nx*l] = (data_t (*)[2*nx*l]) an;
    data_t * bucket;
    const data_t* mod[3] = {model, model+nx*nz, model+2*nx*nz};

    // stretching function g = (p+1)/(2L) * vmax * log(1/R) * (z/L)^p = g0 * (z/L)^p   (eq. 24)

    data_t dt2 = dt*dt;

    // first and second eq. in eq. (21) along with eq. (20)
    Dzz_var<expr2>(false, pin[0], pout[0], nx, nz, dz, ixmin, ixmax, 0, l, mod, 1.0); // d/dz(mu.dux/dz)
    if (version==1) Dzz_var<expr4a>(false, pin[1], pout[1], nx, nz, dz, ixmin, ixmax, 0, l, mod, 1.0); // d/dz((la+2mu).duz/dz)
    else{
        Dzz_var<expr4b>(false, pin[1], pout[1], nx, nz, dz, ixmin, ixmax, 0, l, mod, 1.0); // d/dz((2mu).duz/dz)
        mult_Dz(true, pu_z[1], pout[1], nx, nz, dz, ixmin, ixmax, 0, l, mod[0], 1.0); // d/dz(lambda.duz/dz)
    }
    
    #pragma omp parallel for
    for (int ix=ixmin; ix<ixmax; ix++){
        for (int iz=0; iz<l; iz++){

            int i1 = ix*l + iz;
            int j = ix*nz + iz;
            int i2 = (nx+ix)* l + iz;

            pan[0][i1] = (dt2*pout[0][j]/mod[2][j] - dt2*g[iz]*g[iz]*pac[0][i1] + 2*pac[0][i1]-pan[0][i1] + g[iz]*dt*pan[0][i1])/(1+g[iz]*dt);
            pan[1][i1] = (gprime[iz]*dt2*mod[1][j]*pu_z[0][j]/mod[2][j] - dt2*g[iz]*g[iz]*pac[1][i1] + 2*pac[1][i1]-pan[1][i1]+ g[iz]*dt*pan[1][i1])/(1+g[iz]*dt);
            //s2[i1] = (dt*pan[1][i1] + s2[i1])/(1+g[iz]*dt);
            s2[i1] += dt*(pac[1][i1] - g[iz]*s2[i1]);

            pan[0][i2] = (dt2*pout[1][j]/mod[2][j] - dt2*g[iz]*g[iz]*pac[0][i2] + 2*pac[0][i2]-pan[0][i2] + g[iz]*dt*pan[0][i2])/(1+g[iz]*dt);
            pan[1][i2] = (gprime[iz]*dt2*(mod[0][j]+2*mod[1][j])*pu_z[1][j]/mod[2][j] - dt2*g[iz]*g[iz]*pac[1][i2] + 2*pac[1][i2]-pan[1][i2]+ g[iz]*dt*pan[1][i2])/(1+g[iz]*dt);
            //s2[i2] = (dt*pan[1][i2] + s2[i2])/(1+g[iz]*dt);
            s2[i2] += dt*(pac[1][i2] - g[iz]*s2[i2]);
        }
    }

    // third eq. in eq. (21)
    mult_Dx(false, pu_z[1], pout[0], nx, nz, dx, ixmin, ixmax, 0, l, mod[0], 1.0);
    mult_Dz(true, pu_x[1], pout[0], nx, nz, dz, ixmin, ixmax, 0, l, mod[1], 1.0);
    mult_Dz(false, pu_x[0], pout[1], nx, nz, dz, ixmin, ixmax, 0, l, mod[0], 1.0);
    mult_Dx(true, pu_z[0], pout[1], nx, nz, dx, ixmin, ixmax, 0, l, mod[1], 1.0);

    #pragma omp parallel for
    for (int ix=ixmin; ix<ixmax; ix++){
        for (int iz=0; iz<l; iz++){

            int i1 = ix*l + iz;
            int j = ix*nz + iz;
            int i2 = (nx+ix)* l + iz;

            pan[2][i1] = (dt2*pout[0][j]/mod[2][j] + 2*pac[2][i1]-pan[2][i1] + 0.5*g[iz]*dt*pan[2][i1])/(1+0.5*g[iz]*dt);
            pan[2][i2] = (dt2*pout[1][j]/mod[2][j] + 2*pac[2][i2]-pan[2][i2] + 0.5*g[iz]*dt*pan[2][i2])/(1+0.5*g[iz]*dt);
        }
    }

    // fourth eq. in eq. (21) and sum the split fields
    if (version==1) Dxx_var<expr4a>(false, pin[0], pout[0], nx, nz, dx, ixmin, ixmax, 0, l, mod, 1.0); // d/dx((la+2mu).dux/dx)
    else{
        Dxx_var<expr4b>(false, pin[0], pout[0], nx, nz, dx, ixmin, ixmax, 0, l, mod, 1.0); // d/dx((2mu).dux/dx)
        mult_Dx(true, pu_x[0], pout[0], nx, nz, dx, ixmin, ixmax, 0, l, mod[0], 1.0); // d/dx(lambda.dux/dx)
    }
    Dxx_var<expr2>(false, pin[1], pout[1], nx, nz, dx, ixmin, ixmax, 0, l, mod, 1.0); // d/dx(mu.duz/dx)

    #pragma omp parallel for
    for (int ix=ixmin; ix<ixmax; ix++){
        for (int iz=0; iz<l; iz++){

            int i1 = ix*l + iz;
            int j = ix*nz + iz;
            int i2 = (nx+ix)* l + iz;

            pan[3][i1] = dt2*pout[0][j]/mod[2][j] + 2*pac[3][i1] - pan[3][i1];
            pout[0][j] = pan[0][i1] + s2[i1] + pan[2][i1] + pan[3][i1];

            pan[3][i2] = dt2*pout[1][j]/mod[2][j] + 2*pac[3][i2] - pan[3][i2];
            pout[1][j] = pan[0][i2] + s2[i2] + pan[2][i2] + pan[3][i2];
        }
    }

    bucket = ac;
    ac = an;
    an = bucket;
}

void pml_bottom(const data_t * in, data_t * out, const data_t * u_x, const data_t * u_z,
            data_t * &ac, data_t * &an, data_t * s2, const data_t * model, const data_t * g, const data_t * gprime,
            int nx, int nz, int l, int ixmin, int ixmax, data_t dx, data_t dz, data_t dt, int version)
{
    data_t (*pin) [nx*nz] = (data_t (*)[nx*nz]) in;
    data_t (*pout) [nx*nz] = (data_t (*)[nx*nz]) out;
    data_t (*pu_x) [nx*nz] = (data_t (*)[nx*nz]) u_x;
    data_t (*pu_z) [nx*nz] = (data_t (*)[nx*nz]) u_z;
    data_t (*pac) [2*nx*l] = (data_t (*)[2*nx*l]) ac;
    data_t (*pan) [2*nx*l] = (data_t (*)[2*nx*l]) an;
    data_t * bucket;
    const data_t* mod[3] = {model, model+nx*nz, model+2*nx*nz};

    data_t dt2 = dt*dt;

    Dzz_var<expr2>(false, pin[0], pout[0], nx, nz, dz, ixmin, ixmax, nz-l, nz, mod, 1.0); // d/dz(mu.dux/dz)
    if (version==1) Dzz_var<expr4a>(false, pin[1], pout[1], nx, nz, dz, ixmin, ixmax, nz-l, nz, mod, 1.0); // d/dz((la+2mu).duz/dz)
    else{
        Dzz_var<expr4b>(false, pin[1], pout[1], nx, nz, dz, ixmin, ixmax, nz-l, nz, mod, 1.0); // d/dz((2mu).duz/dz)
        mult_Dz(true, pu_z[1], pout[1], nx, nz, dz, ixmin, ixmax, nz-l, nz, mod[0], 1.0); // d/dz(lambda.duz/dz)
    }
    
    #pragma omp parallel for
    for (int ix=ixmin; ix<ixmax; ix++){
        for (int iz=0; iz<l; iz++){

            int i1 = ix*l + iz;
            int j = ix*nz + nz-1-iz;
            int i2 = (nx+ix)* l + iz;

            pan[0][i1] = (dt2*pout[0][j]/mod[2][j] - dt2*g[iz]*g[iz]*pac[0][i1] + 2*pac[0][i1]-pan[0][i1] + g[iz]*dt*pan[0][i1])/(1+g[iz]*dt);
            pan[1][i1] = (-gprime[iz]*dt2*mod[1][j]*pu_z[0][j]/mod[2][j] - dt2*g[iz]*g[iz]*pac[1][i1] + 2*pac[1][i1]-pan[1][i1]+ g[iz]*dt*pan[1][i1])/(1+g[iz]*dt);
            //s2[i1] = (dt*pan[1][i1] + s2[i1])/(1+g[iz]*dt);
            s2[i1] += dt*(pac[1][i1] - g[iz]*s2[i1]);

            pan[0][i2] = (dt2*pout[1][j]/mod[2][j] - dt2*g[iz]*g[iz]*pac[0][i2] + 2*pac[0][i2]-pan[0][i2] + g[iz]*dt*pan[0][i2])/(1+g[iz]*dt);
            pan[1][i2] = (-gprime[iz]*dt2*(mod[0][j]+2*mod[1][j])*pu_z[1][j]/mod[2][j] - dt2*g[iz]*g[iz]*pac[1][i2] + 2*pac[1][i2]-pan[1][i2]+ g[iz]*dt*pan[1][i2])/(1+g[iz]*dt);
            //s2[i2] = (dt*pan[1][i2] + s2[i2])/(1+g[iz]*dt);
            s2[i2] += dt*(pac[1][i2] - g[iz]*s2[i2]);
        }
    }

    mult_Dx(false, pu_z[1], pout[0], nx, nz, dx, ixmin, ixmax, nz-l, nz, mod[0], 1.0);
    mult_Dz(true, pu_x[1], pout[0], nx, nz, dz, ixmin, ixmax, nz-l, nz, mod[1], 1.0);
    mult_Dz(false, pu_x[0], pout[1], nx, nz, dz, ixmin, ixmax, nz-l, nz, mod[0], 1.0);
    mult_Dx(true, pu_z[0], pout[1], nx, nz, dx, ixmin, ixmax, nz-l, nz, mod[1], 1.0);

    #pragma omp parallel for
    for (int ix=ixmin; ix<ixmax; ix++){
        for (int iz=0; iz<l; iz++){

            int i1 = ix*l + iz;
            int j = ix*nz + nz-1-iz;
            int i2 = (nx+ix)* l + iz;

            pan[2][i1] = (dt2*pout[0][j]/mod[2][j] + 2*pac[2][i1]-pan[2][i1] + 0.5*g[iz]*dt*pan[2][i1])/(1+0.5*g[iz]*dt);
            pan[2][i2] = (dt2*pout[1][j]/mod[2][j] + 2*pac[2][i2]-pan[2][i2] + 0.5*g[iz]*dt*pan[2][i2])/(1+0.5*g[iz]*dt);
        }
    }

    // fourth eq. in eq. (21) and sum the split fields
    if (version==1) Dxx_var<expr4a>(false, pin[0], pout[0], nx, nz, dx, ixmin, ixmax, nz-l, nz, mod, 1.0); // d/dx((la+2mu).dux/dx)
    else{
        Dxx_var<expr4b>(false, pin[0], pout[0], nx, nz, dx, ixmin, ixmax, nz-l, nz, mod, 1.0); // d/dx((2mu).dux/dx)
        mult_Dx(true, pu_x[0], pout[0], nx, nz, dx, ixmin, ixmax, nz-l, nz, mod[0], 1.0); // d/dx(lambda.dux/dx)
    }
    Dxx_var<expr2>(false, pin[1], pout[1], nx, nz, dx, ixmin, ixmax, nz-l, nz, mod, 1.0); // d/dx(mu.duz/dx)

    #pragma omp parallel for
    for (int ix=ixmin; ix<ixmax; ix++){
        for (int iz=0; iz<l; iz++){

            int i1 = ix*l + iz;
            int j = ix*nz + nz-1-iz;
            int i2 = (nx+ix)* l + iz;

            pan[3][i1] = dt2*pout[0][j]/mod[2][j] + 2*pac[3][i1] - pan[3][i1];
            pout[0][j] = pan[0][i1] + s2[i1] + pan[2][i1] + pan[3][i1];

            pan[3][i2] = dt2*pout[1][j]/mod[2][j] + 2*pac[3][i2] - pan[3][i2];
            pout[1][j] = pan[0][i2] + s2[i2] + pan[2][i2] + pan[3][i2];
        }
    }

    bucket = ac;
    ac = an;
    an = bucket;
}

void pml_left(const data_t * in, data_t * out, const data_t * u_x, const data_t * u_z,
            data_t * &ac, data_t * &an, data_t * s2,
            data_t * &t2c, data_t * &t2n, data_t * &t3c, data_t * &t3n, data_t * &s5, data_t * &s6, // additional auxiliary arrays for the top left and bottom left corners
            bool pml_T, bool pml_B, // flags to apply or not the 2D PML in the corners. If true izmin must = 0 and izmax must = nz
            const data_t * model, const data_t * g, const data_t * gprime,
            int nx, int nz, int l, int izmin, int izmax, data_t dx, data_t dz, data_t dt, int version)
{
    data_t (*pin) [nx*nz] = (data_t (*)[nx*nz]) in;
    data_t (*pout) [nx*nz] = (data_t (*)[nx*nz]) out;
    data_t (*pu_x) [nx*nz] = (data_t (*)[nx*nz]) u_x;
    data_t (*pu_z) [nx*nz] = (data_t (*)[nx*nz]) u_z;
    data_t (*pac) [2*nx*l] = (data_t (*)[2*nx*l]) ac;
    data_t (*pan) [2*nx*l] = (data_t (*)[2*nx*l]) an;
    data_t * bucket;
    const data_t* mod[3] = {model, model+nx*nz, model+2*nx*nz};

    data_t dt2 = dt*dt, r=dx/dz;
    int izmin1=izmin, izmax1=izmax;
    if (pml_T) izmin1=l;
    if (pml_B) izmax1=nz-l;

    // first and second eq. in eq. (21) along with eq. (20)
    if (version==1) Dxx_var<expr4a>(false, pin[0], pout[0], nx, nz, dx, 0, l, izmin, izmax, mod, 1.0); // d/dx((la+2mu).dux/dx)
    else{
        Dxx_var<expr4b>(false, pin[0], pout[0], nx, nz, dx, 0, l, izmin, izmax, mod, 1.0); // d/dx((2mu).dux/dx)
        mult_Dx(true, pu_x[0], pout[0], nx, nz, dx, 0, l, izmin, izmax, mod[0], 1.0); // d/dx(lambda.dux/dx)
    }
    Dxx_var<expr2>(false, pin[1], pout[1], nx, nz, dx, 0, l, izmin, izmax, mod, 1.0); // d/dz(mu.duz/dz)

    #pragma omp parallel for
    for (int ix=0; ix<l; ix++){
        for (int iz=izmin; iz<izmax; iz++){

            int i1 = ix*nz + iz;
            int j = i1;
            int i2 = nz*(ix+l) + iz;

            pan[0][i1] = (dt2*pout[0][j]/mod[2][j] - dt2*g[ix]*g[ix]*pac[0][i1] + 2*pac[0][i1]-pan[0][i1] + g[ix]*dt*pan[0][i1])/(1+g[ix]*dt);
            pan[1][i1] = (gprime[ix]*dt2*(mod[0][j]+2*mod[1][j])*pu_x[0][j]/mod[2][j] - dt2*g[ix]*g[ix]*pac[1][i1] + 2*pac[1][i1]-pan[1][i1]+ g[ix]*dt*pan[1][i1])/(1+g[ix]*dt);
            s2[i1] += dt*(pac[1][i1] - g[ix]*s2[i1]);

            pan[0][i2] = (dt2*pout[1][j]/mod[2][j] - dt2*g[ix]*g[ix]*pac[0][i2] + 2*pac[0][i2]-pan[0][i2] + g[ix]*dt*pan[0][i2])/(1+g[ix]*dt);
            pan[1][i2] = (gprime[ix]*dt2*mod[1][j]*pu_x[1][j]/mod[2][j] - dt2*g[ix]*g[ix]*pac[1][i2] + 2*pac[1][i2]-pan[1][i2]+ g[ix]*dt*pan[1][i2])/(1+g[ix]*dt);
            s2[i2] += dt*(pac[1][i2] - g[ix]*s2[i2]);
        }
    }

    // third eq. in eq. (21)
    mult_Dx(false, pu_z[1], pout[0], nx, nz, dx, 0, l, izmin, izmax, mod[0], 1.0);
    mult_Dz(true, pu_x[1], pout[0], nx, nz, dz, 0, l, izmin, izmax, mod[1], 1.0);
    mult_Dz(false, pu_x[0], pout[1], nx, nz, dz, 0, l, izmin, izmax, mod[0], 1.0);
    mult_Dx(true, pu_z[0], pout[1], nx, nz, dx, 0, l, izmin, izmax, mod[1], 1.0);
    #pragma omp parallel for
    for (int ix=0; ix<l; ix++){
        for (int iz=izmin1; iz<izmax1; iz++){

            int i1 = ix*nz + iz;
            int j = i1;
            int i2 = nz*(ix+l) + iz;

            pan[2][i1] = (dt2*pout[0][j]/mod[2][j] + 2*pac[2][i1]-pan[2][i1] + 0.5*g[ix]*dt*pan[2][i1])/(1+0.5*g[ix]*dt);
            pan[2][i2] = (dt2*pout[1][j]/mod[2][j] + 2*pac[2][i2]-pan[2][i2] + 0.5*g[ix]*dt*pan[2][i2])/(1+0.5*g[ix]*dt);
        }
    }

    // apply a 2D PML at the top left corner
    if (pml_T){
        #pragma omp parallel for
        for (int ix=0; ix<l; ix++){
            for (int iz=0; iz<l; iz++){

                int i1 = ix*nz + iz;
                int j = i1;
                int i2 = nz*(ix+l) + iz;

                pan[2][i1] = (dt2*pout[0][j]/mod[2][j] + 2*pac[2][i1]-pan[2][i1] + 0.5*(g[ix]+g[iz]*r)*dt*pan[2][i1] - dt2*g[ix]*g[iz]*r*pac[2][i1])/(1+0.5*(g[ix]+g[iz]*r)*dt);
                pan[2][i2] = (dt2*pout[1][j]/mod[2][j] + 2*pac[2][i2]-pan[2][i2] + 0.5*(g[ix]+g[iz]*r)*dt*pan[2][i2] - dt2*g[ix]*g[iz]*r*pac[2][i2])/(1+0.5*(g[ix]+g[iz]*r)*dt);
            }
        }
    }

    // apply a 2D PML at the bottom left corner
    if (pml_B){
        #pragma omp parallel for
        for (int ix=0; ix<l; ix++){
            for (int iz=0; iz<l; iz++){

                int i1 = ix*nz + nz-1-iz;
                int j = i1;
                int i2 = nz*(ix+l) + nz-1-iz;

                pan[2][i1] = (dt2*pout[0][j]/mod[2][j] + 2*pac[2][i1]-pan[2][i1] + 0.5*(g[ix]+g[iz]*r)*dt*pan[2][i1] - dt2*g[ix]*g[iz]*r*pac[2][i1])/(1+0.5*(g[ix]+g[iz]*r)*dt);
                pan[2][i2] = (dt2*pout[1][j]/mod[2][j] + 2*pac[2][i2]-pan[2][i2] + 0.5*(g[ix]+g[iz]*r)*dt*pan[2][i2] - dt2*g[ix]*g[iz]*r*pac[2][i2])/(1+0.5*(g[ix]+g[iz]*r)*dt);
            }
        }
    }

    // fourth eq. in eq. (21) and sum the split fields
    Dzz_var<expr2>(false, pin[0], pout[0], nx, nz, dz, 0, l, izmin, izmax, mod, 1.0); // d/dz(mu.dux/dz)
    if (version==1) Dzz_var<expr4a>(false, pin[1], pout[1], nx, nz, dz, 0, l, izmin, izmax, mod, 1.0); // d/dz((la+2mu).duz/dz)
    else{
        Dzz_var<expr4b>(false, pin[1], pout[1], nx, nz, dz, 0, l, izmin, izmax, mod, 1.0); // d/dz((2mu).duz/dz)
        mult_Dz(true, pu_z[1], pout[1], nx, nz, dz, 0, l, izmin, izmax, mod[0], 1.0); // d/dz(lambda.duz/dz)
    }

    #pragma omp parallel for
    for (int ix=0; ix<l; ix++){
        for (int iz=izmin1; iz<izmax1; iz++){

            int i1 = ix*nz + iz;
            int j = i1;
            int i2 = nz*(ix+l) + iz;

            pan[3][i1] = dt2*pout[0][j]/mod[2][j] + 2*pac[3][i1] - pan[3][i1];
            pout[0][j] = pan[0][i1] + s2[i1] + pan[2][i1] + pan[3][i1];

            pan[3][i2] = dt2*pout[1][j]/mod[2][j] + 2*pac[3][i2] - pan[3][i2];
            pout[1][j] = pan[0][i2] + s2[i2] + pan[2][i2] + pan[3][i2];
        }
    }

    // Corner PML will split the fourth equation in eq. (21) into two equations
    if (pml_T){
        #pragma omp parallel for
        for (int ix=0; ix<l; ix++){
            for (int iz=0; iz<l; iz++){

                int i1 = ix*nz + iz;
                int k1 = ix*l + iz;
                int j = i1;
                int i2 = nz*(ix+l) + iz;
                int k2 = l*(ix+l) + iz;

                pan[3][i1] = (dt2*pout[0][j]/mod[2][j] + 2*pac[3][i1]-pan[3][i1] - dt2*g[iz]*r*g[iz]*r*pac[3][i1]  + g[iz]*r*dt*pan[3][i1])/(1+g[iz]*r*dt);
                t2n[k1] = (gprime[iz]*r*r*dt2*mod[1][j]*pu_z[0][j]/mod[2][j] - dt2*g[iz]*r*g[iz]*r*t2c[k1] + 2*t2c[k1]-t2n[k1]+ g[iz]*r*dt*t2n[k1])/(1+g[iz]*r*dt);
                s5[k1] += dt*(t2c[k1] - g[iz]*r*s5[k1]);
                pout[0][j] = pan[0][i1] + s2[i1] + pan[2][i1] +pan[3][i1] + s5[k1];

                pan[3][i2] = (dt2*pout[1][j]/mod[2][j] + 2*pac[3][i2]-pan[3][i2] - dt2*g[iz]*r*g[iz]*r*pac[3][i2]  + g[iz]*r*dt*pan[3][i2])/(1+g[iz]*r*dt);
                t2n[k2] = (gprime[iz]*r*r*dt2*(mod[0][j]+2*mod[1][j])*pu_z[1][j]/mod[2][j] - dt2*g[iz]*r*g[iz]*r*t2c[k2] + 2*t2c[k2]-t2n[k2]+ g[iz]*r*dt*t2n[k2])/(1+g[iz]*r*dt);
                s5[k2] += dt*(t2c[k2] - g[iz]*r*s5[k2]);
                pout[1][j] = pan[0][i2] + s2[i2] + pan[2][i2] +pan[3][i2] + s5[k2];
            }
        }
    }

    if (pml_B){
        #pragma omp parallel for
        for (int ix=0; ix<l; ix++){
            for (int iz=0; iz<l; iz++){

                int i1 = ix*nz + nz-1-iz;
                int k1 = ix*l + l-1-iz;
                int j = i1;
                int i2 = nz*(ix+l) + nz-1-iz;
                int k2 = l*(ix+l) + l-1-iz;

                pan[3][i1] = (dt2*pout[0][j]/mod[2][j] + 2*pac[3][i1]-pan[3][i1] - dt2*g[iz]*r*g[iz]*r*pac[3][i1]  + g[iz]*r*dt*pan[3][i1])/(1+g[iz]*r*dt);
                t3n[k1] = (-gprime[iz]*r*r*dt2*mod[1][j]*pu_z[0][j]/mod[2][j] - dt2*g[iz]*r*g[iz]*r*t3c[k1] + 2*t3c[k1]-t3n[k1]+ g[iz]*r*dt*t3n[k1])/(1+g[iz]*r*dt);
                s6[k1] += dt*(t3c[k1] - g[iz]*r*s6[k1]);
                pout[0][j] = pan[0][i1] + s2[i1] + pan[2][i1] +pan[3][i1] + s6[k1];

                pan[3][i2] = (dt2*pout[1][j]/mod[2][j] + 2*pac[3][i2]-pan[3][i2] - dt2*g[iz]*r*g[iz]*r*pac[3][i2]  + g[iz]*r*dt*pan[3][i2])/(1+g[iz]*r*dt);
                t3n[k2] = (-gprime[iz]*r*r*dt2*(mod[0][j]+2*mod[1][j])*pu_z[1][j]/mod[2][j] - dt2*g[iz]*r*g[iz]*r*t3c[k2] + 2*t3c[k2]-t3n[k2]+ g[iz]*r*dt*t3n[k2])/(1+g[iz]*r*dt);
                s6[k2] += dt*(t3c[k2] - g[iz]*r*s6[k2]);
                pout[1][j] = pan[0][i2] + s2[i2] + pan[2][i2] +pan[3][i2] + s6[k2];
            }
        }
    }

    bucket = ac;
    ac = an;
    an = bucket;
    bucket=t2c; t2c=t2n; t2n=bucket; 
    bucket=t3c; t3c=t3n; t3n=bucket;
}

void pml_right(const data_t * in, data_t * out, const data_t * u_x, const data_t * u_z,
            data_t * &ac, data_t * &an, data_t * s2,
            data_t * &t2c, data_t * &t2n, data_t * &t3c, data_t * &t3n, data_t * &s5, data_t * &s6, // additional auxiliary arrays for the top right and bottom right corners
            bool pml_T, bool pml_B, // flags to apply or not the 2D PML in the corners. If true izmin must = 0 and izmax must = nz
            const data_t * model, const data_t * g, const data_t * gprime,
            int nx, int nz, int l, int izmin, int izmax, data_t dx, data_t dz, data_t dt, int version)
{
    data_t (*pin) [nx*nz] = (data_t (*)[nx*nz]) in;
    data_t (*pout) [nx*nz] = (data_t (*)[nx*nz]) out;
    data_t (*pu_x) [nx*nz] = (data_t (*)[nx*nz]) u_x;
    data_t (*pu_z) [nx*nz] = (data_t (*)[nx*nz]) u_z;
    data_t (*pac) [2*nx*l] = (data_t (*)[2*nx*l]) ac;
    data_t (*pan) [2*nx*l] = (data_t (*)[2*nx*l]) an;
    data_t * bucket;
    const data_t* mod[3] = {model, model+nx*nz, model+2*nx*nz};

    data_t dt2 = dt*dt, r=dx/dz;
    int izmin1=izmin, izmax1=izmax;
    if (pml_T) izmin1=l;
    if (pml_B) izmax1=nz-l;

    // first and second eq. in eq. (21) along with eq. (20)
    if (version==1) Dxx_var<expr4a>(false, pin[0], pout[0], nx, nz, dx, nx-l, nx, izmin, izmax, mod, 1.0); // d/dx((la+2mu).dux/dx)
    else{
        Dxx_var<expr4b>(false, pin[0], pout[0], nx, nz, dx, nx-l, nx, izmin, izmax, mod, 1.0); // d/dx((2mu).dux/dx)
        mult_Dx(true, pu_x[0], pout[0], nx, nz, dx, nx-l, nx, izmin, izmax, mod[0], 1.0); // d/dx(lambda.dux/dx)
    }
    Dxx_var<expr2>(false, pin[1], pout[1], nx, nz, dx, nx-l, nx, izmin, izmax, mod, 1.0); // d/dz(mu.duz/dz)

    #pragma omp parallel for
    for (int ix=0; ix<l; ix++){
        for (int iz=izmin; iz<izmax; iz++){

            int i1 = ix*nz + iz;
            int j = (nx-1-ix)*nz + iz;
            int i2 = nz*(ix+l) + iz;

            pan[0][i1] = (dt2*pout[0][j]/mod[2][j] - dt2*g[ix]*g[ix]*pac[0][i1] + 2*pac[0][i1]-pan[0][i1] + g[ix]*dt*pan[0][i1])/(1+g[ix]*dt);
            pan[1][i1] = (-gprime[ix]*dt2*(mod[0][j]+2*mod[1][j])*pu_x[0][j]/mod[2][j] - dt2*g[ix]*g[ix]*pac[1][i1] + 2*pac[1][i1]-pan[1][i1]+ g[ix]*dt*pan[1][i1])/(1+g[ix]*dt);
            s2[i1] += dt*(pac[1][i1] - g[ix]*s2[i1]);

            pan[0][i2] = (dt2*pout[1][j]/mod[2][j] - dt2*g[ix]*g[ix]*pac[0][i2] + 2*pac[0][i2]-pan[0][i2] + g[ix]*dt*pan[0][i2])/(1+g[ix]*dt);
            pan[1][i2] = (-gprime[ix]*dt2*mod[1][j]*pu_x[1][j]/mod[2][j] - dt2*g[ix]*g[ix]*pac[1][i2] + 2*pac[1][i2]-pan[1][i2]+ g[ix]*dt*pan[1][i2])/(1+g[ix]*dt);
            s2[i2] += dt*(pac[1][i2] - g[ix]*s2[i2]);
        }
    }

    // third eq. in eq. (21)
    mult_Dx(false, pu_z[1], pout[0], nx, nz, dx, nx-l, nx, izmin, izmax, mod[0], 1.0);
    mult_Dz(true, pu_x[1], pout[0], nx, nz, dz, nx-l, nx, izmin, izmax, mod[1], 1.0);
    mult_Dz(false, pu_x[0], pout[1], nx, nz, dz, nx-l, nx, izmin, izmax, mod[0], 1.0);
    mult_Dx(true, pu_z[0], pout[1], nx, nz, dx, nx-l, nx, izmin, izmax, mod[1], 1.0);
    #pragma omp parallel for
    for (int ix=0; ix<l; ix++){
        for (int iz=izmin1; iz<izmax1; iz++){

            int i1 = ix*nz + iz;
            int j = (nx-1-ix)*nz + iz;
            int i2 = nz*(ix+l) + iz;

            pan[2][i1] = (dt2*pout[0][j]/mod[2][j] + 2*pac[2][i1]-pan[2][i1] + 0.5*g[ix]*dt*pan[2][i1])/(1+0.5*g[ix]*dt);
            pan[2][i2] = (dt2*pout[1][j]/mod[2][j] + 2*pac[2][i2]-pan[2][i2] + 0.5*g[ix]*dt*pan[2][i2])/(1+0.5*g[ix]*dt);
        }
    }

    // apply a 2D PML at the top right corner
    if (pml_T){
        #pragma omp parallel for
        for (int ix=0; ix<l; ix++){
            for (int iz=0; iz<l; iz++){

                int i1 = ix*nz + iz;
                int j = (nx-1-ix)*nz + iz;
                int i2 = nz*(ix+l) + iz;

                pan[2][i1] = (dt2*pout[0][j]/mod[2][j] + 2*pac[2][i1]-pan[2][i1] + 0.5*(g[ix]+g[iz]*r)*dt*pan[2][i1] - dt2*g[ix]*g[iz]*r*pac[2][i1])/(1+0.5*(g[ix]+g[iz]*r)*dt);
                pan[2][i2] = (dt2*pout[1][j]/mod[2][j] + 2*pac[2][i2]-pan[2][i2] + 0.5*(g[ix]+g[iz]*r)*dt*pan[2][i2] - dt2*g[ix]*g[iz]*r*pac[2][i2])/(1+0.5*(g[ix]+g[iz]*r)*dt);
            }
        }
    }

    // apply a 2D PML at the bottom right corner
    if (pml_B){
        #pragma omp parallel for
        for (int ix=0; ix<l; ix++){
            for (int iz=0; iz<l; iz++){

                int i1 = ix*nz + nz-1-iz;
                int j = (nx-1-ix)*nz + nz-1-iz;
                int i2 = nz*(ix+l) + nz-1-iz;

                pan[2][i1] = (dt2*pout[0][j]/mod[2][j] + 2*pac[2][i1]-pan[2][i1] + 0.5*(g[ix]+g[iz]*r)*dt*pan[2][i1] - dt2*g[ix]*g[iz]*r*pac[2][i1])/(1+0.5*(g[ix]+g[iz]*r)*dt);
                pan[2][i2] = (dt2*pout[1][j]/mod[2][j] + 2*pac[2][i2]-pan[2][i2] + 0.5*(g[ix]+g[iz]*r)*dt*pan[2][i2] - dt2*g[ix]*g[iz]*r*pac[2][i2])/(1+0.5*(g[ix]+g[iz]*r)*dt);
            }
        }
    }

    // fourth eq. in eq. (21) and sum the split fields
    Dzz_var<expr2>(false, pin[0], pout[0], nx, nz, dz, nx-l, nx, izmin, izmax, mod, 1.0); // d/dz(mu.dux/dz)
    if (version==1) Dzz_var<expr4a>(false, pin[1], pout[1], nx, nz, dz, nx-l, nx, izmin, izmax, mod, 1.0); // d/dz((la+2mu).duz/dz)
    else{
       Dzz_var<expr4b>(false, pin[1], pout[1], nx, nz, dz, nx-l, nx, izmin, izmax, mod, 1.0); // d/dz((2mu).duz/dz)
        mult_Dz(true, pu_z[1], pout[1], nx, nz, dz, nx-l, nx, izmin, izmax, mod[0], 1.0); // d/dz(lambda.duz/dz)
    }
    #pragma omp parallel for
    for (int ix=0; ix<l; ix++){
        for (int iz=izmin1; iz<izmax1; iz++){

            int i1 = ix*nz + iz;
            int j = (nx-1-ix)*nz + iz;
            int i2 = nz*(ix+l) + iz;

            pan[3][i1] = dt2*pout[0][j]/mod[2][j] + 2*pac[3][i1] - pan[3][i1];
            pout[0][j] = pan[0][i1] + s2[i1] + pan[2][i1] + pan[3][i1];

            pan[3][i2] = dt2*pout[1][j]/mod[2][j] + 2*pac[3][i2] - pan[3][i2];
            pout[1][j] = pan[0][i2] + s2[i2] + pan[2][i2] + pan[3][i2];
        }
    }

    // Corner PML will split the fourth equation in eq. (21) into two equations
    if (pml_T){
        #pragma omp parallel for
        for (int ix=0; ix<l; ix++){
            for (int iz=0; iz<l; iz++){

                int i1 = ix*nz + iz;
                int k1 = ix*l + iz;
                int j = (nx-1-ix)*nz + iz;
                int i2 = nz*(ix+l) + iz;
                int k2 = l*(ix+l) + iz;

                pan[3][i1] = (dt2*pout[0][j]/mod[2][j] + 2*pac[3][i1]-pan[3][i1] - dt2*g[iz]*r*g[iz]*r*pac[3][i1]  + g[iz]*r*dt*pan[3][i1])/(1+g[iz]*r*dt);
                t2n[k1] = (gprime[iz]*r*r*dt2*mod[1][j]*pu_z[0][j]/mod[2][j] - dt2*g[iz]*r*g[iz]*r*t2c[k1] + 2*t2c[k1]-t2n[k1]+ g[iz]*r*dt*t2n[k1])/(1+g[iz]*r*dt);
                s5[k1] += dt*(t2c[k1] - g[iz]*r*s5[k1]);
                pout[0][j] = pan[0][i1] + s2[i1] + pan[2][i1] +pan[3][i1] + s5[k1];

                pan[3][i2] = (dt2*pout[1][j]/mod[2][j] + 2*pac[3][i2]-pan[3][i2] - dt2*g[iz]*r*g[iz]*r*pac[3][i2]  + g[iz]*r*dt*pan[3][i2])/(1+g[iz]*r*dt);
                t2n[k2] = (gprime[iz]*r*r*dt2*(mod[0][j]+2*mod[1][j])*pu_z[1][j]/mod[2][j] - dt2*g[iz]*r*g[iz]*r*t2c[k2] + 2*t2c[k2]-t2n[k2]+ g[iz]*r*dt*t2n[k2])/(1+g[iz]*r*dt);
                s5[k2] += dt*(t2c[k2] - g[iz]*r*s5[k2]);
                pout[1][j] = pan[0][i2] + s2[i2] + pan[2][i2] +pan[3][i2] + s5[k2];
            }
        }
    }

    if (pml_B){
        #pragma omp parallel for
        for (int ix=0; ix<l; ix++){
            for (int iz=0; iz<l; iz++){

                int i1 = ix*nz + nz-1-iz;
                int k1 = ix*l + l-1-iz;
                int j = (nx-1-ix)*nz + nz-1-iz;
                int i2 = nz*(ix+l) + nz-1-iz;
                int k2 = l*(ix+l) + l-1-iz;

                pan[3][i1] = (dt2*pout[0][j]/mod[2][j] + 2*pac[3][i1]-pan[3][i1] - dt2*g[iz]*r*g[iz]*r*pac[3][i1]  + g[iz]*r*dt*pan[3][i1])/(1+g[iz]*r*dt);
                t3n[k1] = (-gprime[iz]*r*r*dt2*mod[1][j]*pu_z[0][j]/mod[2][j] - dt2*g[iz]*r*g[iz]*r*t3c[k1] + 2*t3c[k1]-t3n[k1]+ g[iz]*r*dt*t3n[k1])/(1+g[iz]*r*dt);
                s6[k1] += dt*(t3c[k1] - g[iz]*r*s6[k1]);
                pout[0][j] = pan[0][i1] + s2[i1] + pan[2][i1] +pan[3][i1] + s6[k1];

                pan[3][i2] = (dt2*pout[1][j]/mod[2][j] + 2*pac[3][i2]-pan[3][i2] - dt2*g[iz]*r*g[iz]*r*pac[3][i2]  + g[iz]*r*dt*pan[3][i2])/(1+g[iz]*r*dt);
                t3n[k2] = (-gprime[iz]*r*r*dt2*(mod[0][j]+2*mod[1][j])*pu_z[1][j]/mod[2][j] - dt2*g[iz]*r*g[iz]*r*t3c[k2] + 2*t3c[k2]-t3n[k2]+ g[iz]*r*dt*t3n[k2])/(1+g[iz]*r*dt);
                s6[k2] += dt*(t3c[k2] - g[iz]*r*s6[k2]);
                pout[1][j] = pan[0][i2] + s2[i2] + pan[2][i2] +pan[3][i2] + s6[k2];
            }
        }
    }

    bucket = ac;
    ac = an;
    an = bucket;
    bucket=t2c; t2c=t2n; t2n=bucket; 
    bucket=t3c; t3c=t3n; t3n=bucket;
}
    

void pml_allocate(data_t * &top_ac, data_t * &top_an, data_t * &top_s2, data_t * &bottom_ac, data_t * &bottom_an, data_t * &bottom_s2,
            data_t * &left_ac, data_t * &left_an, data_t * &left_s2, data_t * &left_t2c, data_t * &left_t2n, data_t * &left_t3c, data_t * &left_t3n, data_t * &left_s5, data_t * &left_s6,
            data_t * &right_ac, data_t * &right_an, data_t * &right_s2, data_t * &right_t2c, data_t * &right_t2n, data_t * &right_t3c, data_t * &right_t3n, data_t * &right_s5, data_t * &right_s6,
            bool pml_T, bool pml_B, bool pml_L, bool pml_R, int nx, int nz, int l)
{
    if (pml_T){
        top_ac = new data_t[4*2*nx*l];
        top_an = new data_t[4*2*nx*l];
        top_s2 = new data_t[2*nx*l];
        memset(top_ac,0,sizeof(data_t)*4*2*nx*l);
        memset(top_an,0,sizeof(data_t)*4*2*nx*l);
        memset(top_s2 ,0,sizeof(data_t)*2*nx*l);
    }
    if (pml_B){
        bottom_ac = new data_t[4*2*nx*l];
        bottom_an = new data_t[4*2*nx*l];
        bottom_s2 = new data_t[2*nx*l];
        memset(bottom_ac,0,sizeof(data_t)*4*2*nx*l);
        memset(bottom_an,0,sizeof(data_t)*4*2*nx*l);
        memset(bottom_s2 ,0,sizeof(data_t)*2*nx*l);
    }
    if (pml_L){
        left_ac = new data_t[4*2*nx*l];
        left_an = new data_t[4*2*nx*l];
        left_s2 = new data_t[2*nx*l];      
        memset(left_ac,0,sizeof(data_t)*4*2*nx*l);
        memset(left_an,0,sizeof(data_t)*4*2*nx*l);
        memset(left_s2 ,0,sizeof(data_t)*2*nx*l);
        left_t2c = new data_t[2*l*l];
        left_t2n = new data_t[2*l*l];
        left_s5 = new data_t[2*l*l];
        memset(left_t2c,0,sizeof(data_t)*2*l*l);
        memset(left_t2n,0,sizeof(data_t)*2*l*l);
        memset(left_s5,0,sizeof(data_t)*2*l*l);
        left_t3c = new data_t[2*l*l];
        left_t3n = new data_t[2*l*l];
        left_s6 = new data_t[2*l*l];
        memset(left_t3c,0,sizeof(data_t)*2*l*l);
        memset(left_t3n,0,sizeof(data_t)*2*l*l);
        memset(left_s6,0,sizeof(data_t)*2*l*l);
    }
    if (pml_R){
        right_ac = new data_t[4*2*nx*l];
        right_an = new data_t[4*2*nx*l];
        right_s2 = new data_t[2*nx*l];      
        memset(right_ac,0,sizeof(data_t)*4*2*nx*l);
        memset(right_an,0,sizeof(data_t)*4*2*nx*l);
        memset(right_s2 ,0,sizeof(data_t)*2*nx*l);
        right_t2c = new data_t[2*l*l];
        right_t2n = new data_t[2*l*l];
        right_s5 = new data_t[2*l*l];
        memset(right_t2c,0,sizeof(data_t)*2*l*l);
        memset(right_t2n,0,sizeof(data_t)*2*l*l);
        memset(right_s5,0,sizeof(data_t)*2*l*l);
        right_t3c = new data_t[2*l*l];
        right_t3n = new data_t[2*l*l];
        right_s6 = new data_t[2*l*l];
        memset(right_t3c,0,sizeof(data_t)*2*l*l);
        memset(right_t3n,0,sizeof(data_t)*2*l*l);
        memset(right_s6,0,sizeof(data_t)*2*l*l);

    }
}

void pml_deallocate(data_t * &top_ac, data_t * &top_an, data_t * &top_s2, data_t * &bottom_ac, data_t * &bottom_an, data_t * &bottom_s2,
            data_t * &left_ac, data_t * &left_an, data_t * &left_s2, data_t * &left_t2c, data_t * &left_t2n, data_t * &left_t3c, data_t * &left_t3n, data_t * &left_s5, data_t * &left_s6,
            data_t * &right_ac, data_t * &right_an, data_t * &right_s2, data_t * &right_t2c, data_t * &right_t2n, data_t * &right_t3c, data_t * &right_t3n, data_t * &right_s5, data_t * &right_s6,
            bool pml_T, bool pml_B, bool pml_L, bool pml_R, int nx, int nz, int l)
{
    if (pml_T){
        delete [] top_ac ;
        delete [] top_an ;
        delete [] top_s2 ;
    }
    if (pml_B){
        delete [] bottom_ac ;
        delete [] bottom_an ;
        delete [] bottom_s2 ;
    }
    if (pml_L){
        delete [] left_ac ;
        delete [] left_an ;
        delete [] left_s2 ;
        if (pml_T){
            delete [] left_t2c ;
            delete [] left_t2n ;
            delete [] left_s5  ;
        }
        if (pml_B){
            delete [] left_t3c ;
            delete [] left_t3n ;
            delete [] left_s6  ;
        }
    }
    if (pml_R){
        delete [] right_ac ;
        delete [] right_an ;
        delete [] right_s2 ;
        if (pml_T){
            delete [] right_t2c ;
            delete [] right_t2n ;
            delete [] right_s5  ;
        }
        if (pml_B){
            delete [] right_t3c ;
            delete [] right_t3n ;
            delete [] right_s6  ;
        }
    }
}