#pragma once

#include "vecReg.hpp"

typedef data_t (*expr)(const data_t ** par, int i);
typedef data_t (*taper)(data_t x);

// spatial quadrature operator H (or its inverse)
void applyHz(bool inv, bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax);
void applyHx(bool inv, bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax);

// first derivative operators
void Dz(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax);
extern "C" void Dz_ispc(int add, const data_t* pmod, data_t* out, data_t scale, int i1, int izmin, int izmax);

void Dx(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax);
extern "C" void Dx_ispc(int add, const data_t* pmod, data_t* out, data_t scale, int i1, int nz, int izmin, int izmax);

// first derivative operators pre-multiply with variable parameter and a constant
void mult_Dz(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax, const data_t * par, data_t a);
extern "C" void mult_Dz_ispc(int add, const data_t* pmod, data_t* out, data_t scale, int i1, int izmin, int izmax, const data_t * par);

void mult_Dx(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax, const data_t * par, data_t a);
extern "C" void mult_Dx_ispc(int add, const data_t* pmod, data_t* out, data_t scale, int i1, int nz, int izmin, int izmax, const data_t * par);

// second derivative operators
void Dzz(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax);
void Dxx(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax);

// second derivative operators with variable parameters, defined as template function to accomodate different expressions of parameters
template<expr f>
void Dzz_var(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax, const data_t ** par, data_t a){

    data_t coef[10] = {-3.0/4, -5.0/6, -1.0/24, 1.0/6, 1.0/2, 1.0/2, 1.0/6, -1.0/8, 1.0/6, -1.0/8};
    data_t bnd_coef[384] = {920.0/289,-59.0/68,-7549318.0/34157643,-17440994.0/184112825,0,0,0,0,-1740.0/289,0,295314153.0/366719282,262545878.0/1218454145,0,0,0,0,1128.0/289,59.0/68,-19250923.0/26254840,-12192537.0/324892213,0,0,0,0,-308.0/289,0,42283069.0/254173229,-43013531.0/427521546,0,0,0,0,0,0,-18700293.0/1355757959,18700293.0/1355757959,0,0,0,0,0,0,-3.0/833,3.0/833,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                12.0/17,0,89562243.0/385991318,16914691.0/272440171,0,0,0,0,-59.0/68,0,-47979680.0/48852831,-18034913.0/120051851,0,0,0,0,2.0/17,0,299262497.0/373256703,16156647.0/200473586,0,0,0,0,3.0/68,0,-14723969.0/177744748,22633571.0/584543661,0,0,0,0,0,0,46802031.0/1628311862,-46802031.0/1628311862,0,0,0,0,0,0,441.0/181507,-441.0/181507,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                -96.0/731,59.0/172,-47632751.0/164317197,-5570723.0/375470930,0,0,0,0,118.0/731,0,54282598.0/49343777,23802793.0/215253532,0,0,0,0,-16.0/731,-59.0/172,-39119273.0/25083370,-35971870.0/61324629,-26254.0/557679,0,0,0,-6.0/731,0,360454121.0/368940022,17254963.0/80047776,1500708.0/7993399,0,0,0,0,0,-18024731.0/79673021,24178273.0/88099647,-26254.0/185893,0,0,0,0,0,-870707.0/620833782,960119.0/1147305747,13564.0/23980197,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                -36.0/833,0,54741803.0/948483020,-13602043.0/389676498,0,0,0,0,177.0/3332,0,-35820026.0/359121865,24921773.0/534548210,0,0,0,0,-6.0/833,0,813284372.0/948584067,30057666.0/158897885,1500708.0/9108757,0,0,0,-9.0/3332,0,-95056924.0/90903639,-23417695.0/47008088,-7476412.0/9108757,-2.0/49,0,0,0,0,23159719.0/99948527,110687545.0/265515812,4502124.0/9108757,8.0/49,0,0,0,0,-3671038.0/2687426923,-1063649.0/8893843,1473580.0/9108757,-6.0/49,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,-9437957.0/1931986386,9437957.0/1931986386,0,0,0,0,0,0,17289851.0/489388053,-17289851.0/489388053,0,0,0,0,0,0,-66355919.0/327412264,24178273.0/98343792,-564461.0/4461432,0,0,0,0,0,17638343.0/74566894,103749401.0/243793650,375177.0/743572,1.0/6,0,0,0,0,-19321801.0/295845927,-50677283.0/62943042,-280535.0/371786,-5.0/6,-1.0/24,0,0,0,5130905.0/5183662662,35039615.0/213452232,1118749.0/2230716,1.0/2,1.0/6,0,0,0,0,0,-1.0/8,1.0/6,-1.0/8,0,0,0,0,0,0,0,0,0,
                0,0,-1.0/784,1.0/784,0,0,0,0,0,0,8673.0/2904112,-8673.0/2904112,0,0,0,0,0,0,-403062.0/320810033,960119.0/1280713392,3391.0/6692148,0,0,0,0,0,-1920494.0/1377228165,-1063649.0/8712336,368395.0/2230716,-1.0/8,0,0,0,0,5130905.0/5183662662,35039615.0/213452232,1118749.0/2230716,1.0/2,1.0/6,0,0,0,-117380.0/2351569839,-3290636.0/80044587,-5580181.0/6692148,-3.0/4,-5.0/6,-1.0/24,0,0,0,0,1.0/6,1.0/2,1.0/2,1.0/6,0,0,0,0,0,-1.0/8,1.0/6,-1.0/8
    };
    
    int nc1=6, nc2=8, nc3=8;
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);
    data_t  val=0;
    data_t d2 = d*d;
    
    #pragma omp parallel for private(val)
    for (int ix = ixmin; ix < ixmax; ix++){
        int i1=ix*nz;

        // apply the operator near the top boundary if included
        for (int iz=izmin; iz<nc1; iz++){
            out[i1+iz] = add*out[i1+iz];
            int i2=iz*nc2*nc3;
            for (int j = 0; j<nc2; j++){
                val=0;
                for (int k = 0; k<nc3; k++){
                     val += bnd_coef[i2+j*nc3+k] * f(par,i1+k);
                }
                out[i1+iz] += a * val * in[i1+j] / d2;
            }
        }

        // apply the operator near the bottom boundary if included
        for (int iz=nz-izmax; iz<nc1; iz++){
            out[i1+nz-1-iz] = add*out[i1+nz-1-iz];
            int i2=iz*nc2*nc3;
            for (int j = 0; j<nc2; j++){
                val=0;
                for (int k = 0; k<nc3; k++){
                     val += bnd_coef[i2+j*nc3+k] * f(par,i1+nz-1-k);
                }
                out[i1+nz-1-iz] +=  a * val * in[i1+nz-1-j]/d2;
            }
        }

        // apply the operator in the interior
        for (int iz = izminb; iz<izmaxb; iz++){
            out[i1+iz] = add*out[i1+iz]
                                        + a/d2 * ( (coef[0]*f(par,i1+iz)+coef[1]*(f(par,i1+iz-1)+f(par,i1+iz+1))+coef[2]*(f(par,i1+iz-2)+f(par,i1+iz+2)))*in[i1+iz] 
                                        + (coef[3]*f(par,i1+iz-2)+coef[4]*f(par,i1+iz-1)+coef[5]*f(par,i1+iz)+coef[6]*f(par,i1+iz+1))*in[i1+iz-1]
                                        + (coef[3]*f(par,i1+iz+2)+coef[4]*f(par,i1+iz+1)+coef[5]*f(par,i1+iz)+coef[6]*f(par,i1+iz-1))*in[i1+iz+1]
                                        + (coef[7]*f(par,i1+iz-2)+coef[8]*f(par,i1+iz-1)+coef[9]*f(par,i1+iz))*in[i1+iz-2]
                                        + (coef[7]*f(par,i1+iz+2)+coef[8]*f(par,i1+iz+1)+coef[9]*f(par,i1+iz))*in[i1+iz+2]
                                        );
        }
    }
}

template<expr f>
void Dxx_var(bool add, const data_t* in, data_t* out, int nx, int nz, data_t d, int ixmin, int ixmax, int izmin, int izmax, const data_t ** par, data_t a){

    data_t coef[10] = {-3.0/4, -5.0/6, -1.0/24, 1.0/6, 1.0/2, 1.0/2, 1.0/6, -1.0/8, 1.0/6, -1.0/8};
    data_t bnd_coef[384] = {920.0/289,-59.0/68,-7549318.0/34157643,-17440994.0/184112825,0,0,0,0,-1740.0/289,0,295314153.0/366719282,262545878.0/1218454145,0,0,0,0,1128.0/289,59.0/68,-19250923.0/26254840,-12192537.0/324892213,0,0,0,0,-308.0/289,0,42283069.0/254173229,-43013531.0/427521546,0,0,0,0,0,0,-18700293.0/1355757959,18700293.0/1355757959,0,0,0,0,0,0,-3.0/833,3.0/833,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                12.0/17,0,89562243.0/385991318,16914691.0/272440171,0,0,0,0,-59.0/68,0,-47979680.0/48852831,-18034913.0/120051851,0,0,0,0,2.0/17,0,299262497.0/373256703,16156647.0/200473586,0,0,0,0,3.0/68,0,-14723969.0/177744748,22633571.0/584543661,0,0,0,0,0,0,46802031.0/1628311862,-46802031.0/1628311862,0,0,0,0,0,0,441.0/181507,-441.0/181507,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                -96.0/731,59.0/172,-47632751.0/164317197,-5570723.0/375470930,0,0,0,0,118.0/731,0,54282598.0/49343777,23802793.0/215253532,0,0,0,0,-16.0/731,-59.0/172,-39119273.0/25083370,-35971870.0/61324629,-26254.0/557679,0,0,0,-6.0/731,0,360454121.0/368940022,17254963.0/80047776,1500708.0/7993399,0,0,0,0,0,-18024731.0/79673021,24178273.0/88099647,-26254.0/185893,0,0,0,0,0,-870707.0/620833782,960119.0/1147305747,13564.0/23980197,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                -36.0/833,0,54741803.0/948483020,-13602043.0/389676498,0,0,0,0,177.0/3332,0,-35820026.0/359121865,24921773.0/534548210,0,0,0,0,-6.0/833,0,813284372.0/948584067,30057666.0/158897885,1500708.0/9108757,0,0,0,-9.0/3332,0,-95056924.0/90903639,-23417695.0/47008088,-7476412.0/9108757,-2.0/49,0,0,0,0,23159719.0/99948527,110687545.0/265515812,4502124.0/9108757,8.0/49,0,0,0,0,-3671038.0/2687426923,-1063649.0/8893843,1473580.0/9108757,-6.0/49,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,-9437957.0/1931986386,9437957.0/1931986386,0,0,0,0,0,0,17289851.0/489388053,-17289851.0/489388053,0,0,0,0,0,0,-66355919.0/327412264,24178273.0/98343792,-564461.0/4461432,0,0,0,0,0,17638343.0/74566894,103749401.0/243793650,375177.0/743572,1.0/6,0,0,0,0,-19321801.0/295845927,-50677283.0/62943042,-280535.0/371786,-5.0/6,-1.0/24,0,0,0,5130905.0/5183662662,35039615.0/213452232,1118749.0/2230716,1.0/2,1.0/6,0,0,0,0,0,-1.0/8,1.0/6,-1.0/8,0,0,0,0,0,0,0,0,0,
                0,0,-1.0/784,1.0/784,0,0,0,0,0,0,8673.0/2904112,-8673.0/2904112,0,0,0,0,0,0,-403062.0/320810033,960119.0/1280713392,3391.0/6692148,0,0,0,0,0,-1920494.0/1377228165,-1063649.0/8712336,368395.0/2230716,-1.0/8,0,0,0,0,5130905.0/5183662662,35039615.0/213452232,1118749.0/2230716,1.0/2,1.0/6,0,0,0,-117380.0/2351569839,-3290636.0/80044587,-5580181.0/6692148,-3.0/4,-5.0/6,-1.0/24,0,0,0,0,1.0/6,1.0/2,1.0/2,1.0/6,0,0,0,0,0,-1.0/8,1.0/6,-1.0/8
    };

    int nc1=6, nc2=8, nc3=8;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);
    data_t  val=0;
    data_t d2 = d*d;
    
    // apply the operator near the left boundary if included
    #pragma omp parallel for private( val)
    for (int ix = ixmin; ix < nc1; ix++){
        int i1=ix*nz;
        int i2=ix*nc2*nc3;
        for (int iz=izmin; iz<izmax; iz++){
            out[i1+iz] = add*out[i1+iz];
            for (int j = 0; j<nc2; j++){
                val=0;
                for (int k = 0; k<nc3; k++){
                    val += bnd_coef[i2+j*nc3+k] * f(par,k*nz+iz);
                }
                out[i1+iz] += a * val * in[j*nz+iz] / d2;
            }
        }
    }

    // apply the operator near the right boundary if included
    #pragma omp parallel for private( val)
    for (int ix = nx-ixmax; ix < nc1; ix++){
        int i1=(nx-1-ix)*nz;
        int i2=ix*nc2*nc3;
        for (int iz=izmin; iz<izmax; iz++){
            out[i1+iz] = add*out[i1+iz];
            for (int j = 0; j<nc2; j++){
                val=0;
                for (int k = 0; k<nc3; k++){
                    val += bnd_coef[i2+j*nc3+k] * f(par,(nx-1-k)*nz+iz);
                }
                out[i1+iz] += a * val * in[(nx-1-j)*nz+iz] / d2;
            }
        }
    }

    // apply the operator in the interior
    #pragma omp parallel for private( val)
    for (int ix=ixminb; ix<ixmaxb; ix++){
        for (int iz = izmin; iz<izmax; iz++){
            out[ ix* nz+ iz] = add*out[ ix* nz+ iz] + a/d2 * ( (coef[0]* f(par, ix* nz+ iz)+coef[1]*( f(par,( ix-1)* nz+ iz)+ f(par,( ix+1)* nz+ iz))+coef[2]*( f(par,( ix-2)* nz+ iz)+ f(par,( ix+2)* nz+ iz)))*in[ ix* nz+ iz] 
                                        + (coef[3]* f(par,( ix-2)* nz+ iz)+coef[4]* f(par,( ix-1)* nz+ iz)+coef[5]* f(par, ix* nz+ iz)+coef[6]* f(par,( ix+1)* nz+ iz))*in[( ix-1)* nz+ iz]
                                        + (coef[3]* f(par,( ix+2)* nz+ iz)+coef[4]* f(par,( ix+1)* nz+ iz)+coef[5]* f(par, ix* nz+ iz)+coef[6]* f(par,( ix-1)* nz+ iz))*in[( ix+1)* nz+ iz]
                                        + (coef[7]* f(par,( ix-2)* nz+ iz)+coef[8]* f(par,( ix-1)* nz+ iz)+coef[9]* f(par, ix* nz+ iz))*in[( ix-2)* nz+ iz]
                                        + (coef[7]* f(par,( ix+2)* nz+ iz)+coef[8]* f(par,( ix+1)* nz+ iz)+coef[9]* f(par, ix* nz+ iz))*in[( ix+2)* nz+ iz]
                                        );
        }
    }
}

// Neumann SAT to impose free surface BC for elastic WE
// defined as template function to accomodate different expressions of parameters
template<expr f1, expr f2>
void esat_neumann_top(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, const data_t ** par, data_t a){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);

    data_t sumx, sumz;

    // SAT = - H-1 (-f1.Dx.in1 + f2.Sz.in2)_0
    // Sz is the boundary derivative operator pointing outwards

    // top left
    for (int ix=ixmin; ix<nc1; ix++){

        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx += bnd_coef[ix*nc2+j] * in[0][j*nz];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in[1][ix*nz+iz];
        }
        out[ix*nz] = add*out[ix*nz] - a /(dz * h0) * (-f1(par,ix*nz)*sumx/dx + f2(par,ix*nz)*sumz/dz);
    }

    // top right
    for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){

        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx -= bnd_coef[ix*nc2+j] * in[0][(nx-1-j)*nz];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in[1][(nx-1-ix)*nz+iz];
        }
        out[(nx-1-ix)*nz] = add*out[(nx-1-ix)*nz] - a /(dz * h0) * (-f1(par,(nx-1-ix)*nz)*sumx/dx + f2(par,(nx-1-ix)*nz)*sumz/dz);
    }

    // top middle
    #pragma omp parallel for private(sumx, sumz)
    for (int ix=ixminb; ix<ixmaxb; ix++){

       // (Dx.in1)_0
        sumx=0;
        for (int j=1; j<=2; j++){
            sumx += coef[j-1] * (in[0][(ix+j)*nz]-in[0][(ix-j)*nz]);
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in[1][ix*nz+iz];
        }

        out[ix*nz] = add*out[ix*nz] - a /(dz * h0) * (-f1(par,ix*nz)*sumx/dx + f2(par,ix*nz)*sumz/dz);
    }
}

template<expr f1, expr f2>
void esat_neumann_bottom(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, const data_t ** par, data_t a){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);

    data_t sumx, sumz;

    // SAT = - H-1 (f1.Dx.in1 + f2.Sz.in2)_0

    // bottom left
    for (int ix=ixmin; ix<nc1; ix++){

        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx += bnd_coef[ix*nc2+j] * in[0][j*nz+nz-1];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in[1][ix*nz+nz-1-iz];
        }
        out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - a /(dz * h0) * (f1(par,ix*nz+nz-1)*sumx/dx + f2(par,ix*nz+nz-1)*sumz/dz);
    }

    // bottom right
    for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){

        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx -= bnd_coef[ix*nc2+j] * in[0][(nx-1-j)*nz+nz-1];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in[1][(nx-1-ix)*nz+nz-1-iz];
        }
        out[(nx-1-ix)*nz+nz-1] = add*out[(nx-1-ix)*nz+nz-1] - a /(dz * h0) * (f1(par,(nx-1-ix)*nz+nz-1)*sumx/dx + f2(par,(nx-1-ix)*nz+nz-1)*sumz/dz);
    }

    // bottom middle
    #pragma omp parallel for private(sumx, sumz)
    for (int ix=ixminb; ix<ixmaxb; ix++){

       // (Dx.in1)_0
        sumx=0;
        for (int j=1; j<=2; j++){
            sumx += coef[j-1] * (in[0][(ix+j)*nz+nz-1]-in[0][(ix-j)*nz+nz-1]);
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in[1][ix*nz+nz-1-iz];
        }

        out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - a /(dz * h0) * (f1(par,ix*nz+nz-1)*sumx/dx + f2(par,ix*nz+nz-1)*sumz/dz);
    }
}

template<expr f1, expr f2>
void esat_neumann_left(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int izmin, int izmax, const data_t ** par, data_t a){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);

    data_t sumx, sumz;

    // SAT = - H-1 (-f1.Dz.in1 + f2.Sx.in2)_0
    // Sx is the boundary derivative operator pointing outwards

    // left top
    for (int iz=izmin; iz<nc1; iz++){

        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz += bnd_coef[iz*nc2+j] * in[0][j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in[1][ix*nz+iz];
        }
        out[iz] = add*out[iz] - a /(dx * h0) * (-f1(par,iz)*sumz/dz + f2(par,iz)*sumx/dx);
    }

    // left bottom
    for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz -= bnd_coef[iz*nc2+j] * in[0][nz-1-j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in[1][ix*nz+nz-1-iz];
        }
        out[nz-1-iz] = add*out[nz-1-iz] - a /(dx * h0) * (-f1(par,nz-1-iz)*sumz/dz + f2(par,nz-1-iz)*sumx/dx);
    }

    // left middle
    #pragma omp parallel for private(sumx, sumz)
    for (int iz=izminb; iz<izmaxb; iz++){

       // (Dz.in1)_0
        sumz=0;
        for (int j=1; j<=2; j++){
            sumz += coef[j-1] * (in[0][iz+j]-in[0][iz-j]);
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in[1][ix*nz+iz];
        }

        out[iz] = add*out[iz] - a /(dx * h0) * (-f1(par,iz)*sumz/dz + f2(par,iz)*sumx/dx);
    }
}

template<expr f1, expr f2>
void esat_neumann_right(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int izmin, int izmax, const data_t ** par, data_t a){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);

    data_t sumx, sumz;

    // SAT = - H-1 (f1.Dz.in1 + f2.Sx.in2)_0

    // right top
    for (int iz=izmin; iz<nc1; iz++){

        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz += bnd_coef[iz*nc2+j] * in[0][(nx-1)*nz+j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in[1][(nx-1-ix)*nz+iz];
        }
        out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - a /(dx * h0) * (f1(par,(nx-1)*nz+iz)*sumz/dz + f2(par,(nx-1)*nz+iz)*sumx/dx);
    }

    // right bottom
    for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz -= bnd_coef[iz*nc2+j] * in[0][(nx-1)*nz+nz-1-j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in[1][(nx-1-ix)*nz+nz-1-iz];
        }
        out[(nx-1)*nz+nz-1-iz] = add*out[(nx-1)*nz+nz-1-iz] - a /(dx * h0) * (f1(par,(nx-1)*nz+nz-1-iz)*sumz/dz + f2(par,(nx-1)*nz+nz-1-iz)*sumx/dx);
    }

    // right middle
    #pragma omp parallel for private(sumx, sumz)
    for (int iz=izminb; iz<izmaxb; iz++){

       // (Dz.in1)_0
        sumz=0;
        for (int j=1; j<=2; j++){
            sumz += coef[j-1] * (in[0][(nx-1)*nz+iz+j]-in[0][(nx-1)*nz+iz-j]);
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in[1][(nx-1-ix)*nz+iz];
        }

        out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - a /(dx * h0) * (f1(par,(nx-1)*nz+iz)*sumz/dz + f2(par,(nx-1)*nz+iz)*sumx/dx);
    }
}

// Mixed SAT to impose (locally) absorbing BC for elastic WE
// defined as template function to accomodate different expressions of parameters
template<expr f1, expr f2, expr f3>
void esat_absorbing_top(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int ixmin, int ixmax, const data_t ** par, data_t a){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);

    data_t sumx, sumz;

    // SAT = - H-1 (-f1.Dx.in1 + f2.Sz.in2 -f3.in3/dt)_0
    // f3 is often P or S impedance (rho*Vp or rho*Vs) divided by 2

    // top left
    for (int ix=ixmin; ix<nc1; ix++){

        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx += bnd_coef[ix*nc2+j] * in[0][j*nz];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in[1][ix*nz+iz];
        }
        out[ix*nz] = add*out[ix*nz] - a /(dz * h0) * (-f1(par,ix*nz)*sumx/dx + f2(par,ix*nz)*sumz/dz - f3(par,ix*nz)*in[2][ix*nz]/dt);
        //out[ix*nz] = add*out[ix*nz] - a /(dz * h0) * (-f1(par,ix*nz)*sumx/dx + f2(par,ix*nz)*sumz/dz - f3(par,ix*nz)*2*(in[2][ix*nz]-in[1][ix*nz])/dt);
    }

    // top right
    for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){

        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx -= bnd_coef[ix*nc2+j] * in[0][(nx-1-j)*nz];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in[1][(nx-1-ix)*nz+iz];
        }
        out[(nx-1-ix)*nz] = add*out[(nx-1-ix)*nz] - a /(dz * h0) * (-f1(par,(nx-1-ix)*nz)*sumx/dx + f2(par,(nx-1-ix)*nz)*sumz/dz - f3(par,(nx-1-ix)*nz)*in[2][(nx-1-ix)*nz]/dt);
        //out[(nx-1-ix)*nz] = add*out[(nx-1-ix)*nz] - a /(dz * h0) * (-f1(par,(nx-1-ix)*nz)*sumx/dx + f2(par,(nx-1-ix)*nz)*sumz/dz - f3(par,(nx-1-ix)*nz)*2*(in[2][(nx-1-ix)*nz]-in[2][(nx-1-ix)*nz])/dt);
    }

    // top middle
    #pragma omp parallel for private(sumx, sumz)
    for (int ix=ixminb; ix<ixmaxb; ix++){

       // (Dx.in1)_0
        sumx=0;
        for (int j=1; j<=2; j++){
            sumx += coef[j-1] * (in[0][(ix+j)*nz]-in[0][(ix-j)*nz]);
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in[1][ix*nz+iz];
        }

        out[ix*nz] = add*out[ix*nz] - a /(dz * h0) * (-f1(par,ix*nz)*sumx/dx + f2(par,ix*nz)*sumz/dz - f3(par,ix*nz)*in[2][ix*nz]/dt);
        //out[ix*nz] = add*out[ix*nz] - a /(dz * h0) * (-f1(par,ix*nz)*sumx/dx + f2(par,ix*nz)*sumz/dz - f3(par,ix*nz)*2*(in[2][ix*nz]-in[1][ix*nz])/dt);
    }
}

template<expr f1, expr f2, expr f3>
void esat_absorbing_bottom(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int ixmin, int ixmax, const data_t ** par, data_t a){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int ixminb=std::max(nc1,ixmin);
    int ixmaxb=std::min(nx - nc1, ixmax);

    data_t sumx, sumz;

    // SAT = - H-1 (f1.Dx.in1 + f2.Sz.in2 - f3.in3/dt)_0

    // bottom left
    for (int ix=ixmin; ix<nc1; ix++){

        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx += bnd_coef[ix*nc2+j] * in[0][j*nz+nz-1];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in[1][ix*nz+nz-1-iz];
        }
        out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - a /(dz * h0) * (f1(par,ix*nz+nz-1)*sumx/dx + f2(par,ix*nz+nz-1)*sumz/dz - f3(par,ix*nz+nz-1)*in[2][ix*nz+nz-1]/dt);
        //out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - a /(dz * h0) * (f1(par,ix*nz+nz-1)*sumx/dx + f2(par,ix*nz+nz-1)*sumz/dz - f3(par,ix*nz+nz-1)*2*(in[2][ix*nz+nz-1]-in[1][ix*nz+nz-1])/dt);
    }

    // bottom right
    for (int ix=nx-ixmax; ix<std::min(nc1,nx-ixmin); ix++){

        // (Dx.in1)_0
        sumx=0;
        for (int j=0; j<nc2; j++){
            sumx -= bnd_coef[ix*nc2+j] * in[0][(nx-1-j)*nz+nz-1];
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in[1][(nx-1-ix)*nz+nz-1-iz];
        }
        out[(nx-1-ix)*nz+nz-1] = add*out[(nx-1-ix)*nz+nz-1] - a /(dz * h0) * (f1(par,(nx-1-ix)*nz+nz-1)*sumx/dx + f2(par,(nx-1-ix)*nz+nz-1)*sumz/dz - f3(par,(nx-1-ix)*nz+nz-1)*in[2][(nx-1-ix)*nz+nz-1]/dt);
        //out[(nx-1-ix)*nz+nz-1] = add*out[(nx-1-ix)*nz+nz-1] - a /(dz * h0) * (f1(par,(nx-1-ix)*nz+nz-1)*sumx/dx + f2(par,(nx-1-ix)*nz+nz-1)*sumz/dz - f3(par,(nx-1-ix)*nz+nz-1)*2*(in[2][(nx-1-ix)*nz+nz-1]-in[1][(nx-1-ix)*nz+nz-1])/dt);
    }

    // bottom middle
    #pragma omp parallel for private(sumx, sumz)
    for (int ix=ixminb; ix<ixmaxb; ix++){

       // (Dx.in1)_0
        sumx=0;
        for (int j=1; j<=2; j++){
            sumx += coef[j-1] * (in[0][(ix+j)*nz+nz-1]-in[0][(ix-j)*nz+nz-1]);
        }

        // (Sz.in2)_0
        sumz = 0;
        for (int iz = 0; iz < 4; iz++){
            sumz += scoef[iz] * in[1][ix*nz+nz-1-iz];
        }

        out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - a /(dz * h0) * (f1(par,ix*nz+nz-1)*sumx/dx + f2(par,ix*nz+nz-1)*sumz/dz - f3(par,ix*nz+nz-1)*in[2][ix*nz+nz-1]/dt);
        //out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - a /(dz * h0) * (f1(par,ix*nz+nz-1)*sumx/dx + f2(par,ix*nz+nz-1)*sumz/dz - f3(par,ix*nz+nz-1)*2*(in[2][ix*nz+nz-1]-in[1][ix*nz+nz-1])/dt);
    }
}

template<expr f1, expr f2, expr f3>
void esat_absorbing_left(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int izmin, int izmax, const data_t ** par, data_t a){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);

    data_t sumx, sumz;

    // SAT = - H-1 (-f1.Dz.in1 + f2.Sx.in2 - f3.in3/dt)_0

    // left top
    for (int iz=izmin; iz<nc1; iz++){

        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz += bnd_coef[iz*nc2+j] * in[0][j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in[1][ix*nz+iz];
        }
        out[iz] = add*out[iz] - a /(dx * h0) * (-f1(par,iz)*sumz/dz + f2(par,iz)*sumx/dx - f3(par,iz)*in[2][iz]/dt);
        //out[iz] = add*out[iz] - a /(dx * h0) * (-f1(par,iz)*sumz/dz + f2(par,iz)*sumx/dx - f3(par,iz)*2*(in[2][iz]-in[1][iz])/dt);
    }

    // left bottom
    for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz -= bnd_coef[iz*nc2+j] * in[0][nz-1-j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in[1][ix*nz+nz-1-iz];
        }
        out[nz-1-iz] = add*out[nz-1-iz] - a /(dx * h0) * (-f1(par,nz-1-iz)*sumz/dz + f2(par,nz-1-iz)*sumx/dx - f3(par,nz-1-iz)*in[2][nz-1-iz]/dt);
        //out[nz-1-iz] = add*out[nz-1-iz] - a /(dx * h0) * (-f1(par,nz-1-iz)*sumz/dz + f2(par,nz-1-iz)*sumx/dx - f3(par,nz-1-iz)*2*(in[2][nz-1-iz]-in[1][nz-1-iz])/dt);
    }

    // left middle
    #pragma omp parallel for private(sumx, sumz)
    for (int iz=izminb; iz<izmaxb; iz++){

       // (Dz.in1)_0
        sumz=0;
        for (int j=1; j<=2; j++){
            sumz += coef[j-1] * (in[0][iz+j]-in[0][iz-j]);
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in[1][ix*nz+iz];
        }

        out[iz] = add*out[iz] - a /(dx * h0) * (-f1(par,iz)*sumz/dz + f2(par,iz)*sumx/dx - f3(par,iz)*in[2][iz]/dt);
        //out[iz] = add*out[iz] - a /(dx * h0) * (-f1(par,iz)*sumz/dz + f2(par,iz)*sumx/dx - f3(par,iz)*2*(in[2][iz]-in[1][iz])/dt);
    }
}

template<expr f1, expr f2, expr f3>
void esat_absorbing_right(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int izmin, int izmax, const data_t ** par, data_t a){

    data_t coef[2] = {2.0/3,-1.0/12};
    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t scoef[4] = {11.0/6, -3, 1.5, -1.0/3}; 
    data_t h0 = 17.0/48;
    int nc1=4, nc2=6;
    int izminb=std::max(nc1,izmin);
    int izmaxb=std::min(nz - nc1, izmax);

    data_t sumx, sumz;

    // SAT = - H-1 (f1.Dz.in1 + f2.Sx.in2 - f3.in3/dt)_0

    // right top
    for (int iz=izmin; iz<nc1; iz++){

        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz += bnd_coef[iz*nc2+j] * in[0][(nx-1)*nz+j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in[1][(nx-1-ix)*nz+iz];
        }
        out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - a /(dx * h0) * (f1(par,(nx-1)*nz+iz)*sumz/dz + f2(par,(nx-1)*nz+iz)*sumx/dx - f3(par,(nx-1)*nz+iz)*in[2][(nx-1)*nz+iz]/dt);
        //out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - a /(dx * h0) * (f1(par,(nx-1)*nz+iz)*sumz/dz + f2(par,(nx-1)*nz+iz)*sumx/dx - f3(par,(nx-1)*nz+iz)*2*(in[2][(nx-1)*nz+iz]-in[1][(nx-1)*nz+iz])/dt);
    }

    // right bottom
    for (int iz=nz-izmax; iz<std::min(nc1,nz-izmin); iz++){

        // (Dz.in1)_0
        sumz=0;
        for (int j=0; j<nc2; j++){
            sumz -= bnd_coef[iz*nc2+j] * in[0][(nx-1)*nz+nz-1-j];
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in[1][(nx-1-ix)*nz+nz-1-iz];
        }
        out[(nx-1)*nz+nz-1-iz] = add*out[(nx-1)*nz+nz-1-iz] - a /(dx * h0) * (f1(par,(nx-1)*nz+nz-1-iz)*sumz/dz + f2(par,(nx-1)*nz+nz-1-iz)*sumx/dx - f3(par,(nx-1)*nz+nz-1-iz)*in[2][(nx-1)*nz+nz-1-iz]/dt);
        //out[(nx-1)*nz+nz-1-iz] = add*out[(nx-1)*nz+nz-1-iz] - a /(dx * h0) * (f1(par,(nx-1)*nz+nz-1-iz)*sumz/dz + f2(par,(nx-1)*nz+nz-1-iz)*sumx/dx - f3(par,(nx-1)*nz+nz-1-iz)*2*(in[2][(nx-1)*nz+nz-1-iz]-in[1][(nx-1)*nz+nz-1-iz])/dt);
    }

    // right middle
    #pragma omp parallel for private(sumx, sumz)
    for (int iz=izminb; iz<izmaxb; iz++){

       // (Dz.in1)_0
        sumz=0;
        for (int j=1; j<=2; j++){
            sumz += coef[j-1] * (in[0][(nx-1)*nz+iz+j]-in[0][(nx-1)*nz+iz-j]);
        }

        // (Sx.in2)_0
        sumx = 0;
        for (int ix = 0; ix < 4; ix++){
            sumx += scoef[ix] * in[1][(nx-1-ix)*nz+iz];
        }

        out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - a /(dx * h0) * (f1(par,(nx-1)*nz+iz)*sumz/dz + f2(par,(nx-1)*nz+iz)*sumx/dx - f3(par,(nx-1)*nz+iz)*in[2][(nx-1)*nz+iz]/dt);
        //out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - a /(dx * h0) * (f1(par,(nx-1)*nz+iz)*sumz/dz + f2(par,(nx-1)*nz+iz)*sumx/dx - f3(par,(nx-1)*nz+iz)*2*(in[2][(nx-1)*nz+iz]-in[1][(nx-1)*nz+iz])/dt);
    }
}


// Apply derivative operator orthogonal to the boundary, to complement the SAT terms
// defined as template function to accomodate different expressions of parameters
template<expr f>
void esat_Dz_top(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, int ixmin, int ixmax, const data_t ** par, data_t a){

    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t h0 = 17.0/48;
    int nc2=6;
    data_t sc = a /(d*d*h0);

    data_t sum;

    // additional SAT = - H-1 (-f.Dz.in)_0
    #pragma omp parallel for private(sum)
    for (int ix=ixmin; ix<ixmax; ix++){

        // -(Dz.in)_0
        sum=0;
        for (int j=0; j<nc2; j++){
            sum -= bnd_coef[j] * in[ix*nz+j];
        }
        out[ix*nz] = add*out[ix*nz] - sc * f(par,ix*nz)*sum;
    }
}

template<expr f>
void esat_Dz_bottom(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, int ixmin, int ixmax, const data_t ** par, data_t a){

    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t h0 = 17.0/48;
    int nc2=6;
    data_t sc = a /(d*d*h0);

    data_t sum;

    // additional SAT = - H-1 (f.Dz.in)_0
    #pragma omp parallel for private(sum)
    for (int ix=ixmin; ix<ixmax; ix++){

        // (Dz.in)_0
        sum=0;
        for (int j=0; j<nc2; j++){
            sum -= bnd_coef[j] * in[ix*nz+nz-1-j];
        }
        out[ix*nz+nz-1] = add*out[ix*nz+nz-1] - sc * f(par,ix*nz+nz-1)*sum;
    }
}

template<expr f>
void esat_Dx_left(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, int izmin, int izmax, const data_t ** par, data_t a){

    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t h0 = 17.0/48;
    int nc2=6;
    data_t sc = a /(d*d*h0);

    data_t sum;

    // additional SAT = - H-1 (-f.Dx.in)_0
    #pragma omp parallel for private(sum)
    for (int iz=izmin; iz<izmax; iz++){

        // -(Dx.in)_0
        sum=0;
        for (int j=0; j<nc2; j++){
            sum -= bnd_coef[j] * in[j*nz+iz];
        }
        out[iz] = add*out[iz] - sc * f(par,iz)*sum;
    }
}

template<expr f>
void esat_Dx_right(bool add, const data_t * in, data_t * out, int nx, int nz, data_t d, int izmin, int izmax, const data_t ** par, data_t a){

    data_t bnd_coef[24] = {-24.0/17,59.0/34,-4.0/17,-3.0/34,0,0,-0.5,0,0.5,0,0,0,4.0/43,-59.0/86,0,59.0/86,-4.0/43,0,3.0/98,0,-59.0/98,0,32.0/49,-4.0/49};
    data_t h0 = 17.0/48;
    int nc2=6;
    data_t sc = a /(d*d*h0);

    data_t sum;

    // additional SAT = - H-1 (f.Dx.in)_0
    #pragma omp parallel for private(sum)
    for (int iz=izmin; iz<izmax; iz++){

        // (Dx.in)_0
        sum=0;
        for (int j=0; j<nc2; j++){
            sum -= bnd_coef[j] * in[(nx-1-j)*nz+iz];
        }
        out[(nx-1)*nz+iz] = add*out[(nx-1)*nz+iz] - sc * f(par,(nx-1)*nz+iz)*sum;
    }
}

// Scale boundaries when absorbing SAT is used. This is needed for the time recursion of wavefields
void esat_scale_boundaries(data_t** in, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, int izmin, int izmax, const data_t** par, data_t dt, bool top, bool bottom, bool left, bool right);
void vtisat_scale_boundaries(data_t** in, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, int izmin, int izmax, const data_t** par, data_t dt, bool top, bool bottom, bool left, bool right);

// Apply cosine^2 taper to damp the wavefield (to use in conjunction with absorbing SAT)
// taper = cos[a * pi/2 * (i - istart)/(iend-istart)]^2 ; 0 <= a <= 1
void taperz(data_t* in, int nx, int nz, int ixmin, int ixmax, int istart, int iend, data_t a);
void taperx(data_t* in, int nx, int nz, int izmin, int izmax, int istart, int iend, data_t a);
extern "C" void taperz_ispc(data_t* in, int i1, int istart, int iend, data_t a);
extern "C" void taperx_ispc(data_t* in, int i1, int n, int izmin, int izmax, data_t a);

// SAT functions for the acoustic wave equation (variable density)
void asat_dirichlet_top(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, const data_t ** par, data_t a);
void asat_dirichlet_bottom(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, const data_t ** par, data_t a);
void asat_dirichlet_left(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int izmin, int izmax, const data_t ** par, data_t a);
void asat_dirichlet_right(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int izmin, int izmax, const data_t ** par, data_t a);
void asat_neumann_top(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, const data_t ** par, data_t a);
void asat_neumann_bottom(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, const data_t ** par, data_t a);
void asat_neumann_left(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int izmin, int izmax, const data_t ** par, data_t a);
void asat_neumann_right(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int izmin, int izmax, const data_t ** par, data_t a);
void asat_absorbing_top(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int ixmin, int ixmax, const data_t ** par, data_t a);
void asat_absorbing_bottom(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int ixmin, int ixmax, const data_t ** par, data_t a);
void asat_absorbing_left(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int izmin, int izmax, const data_t ** par, data_t a);
void asat_absorbing_right(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int izmin, int izmax, const data_t ** par, data_t a);
void asat_scale_boundaries(data_t** in, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, int izmin, int izmax, const data_t** par, data_t dt, bool top, bool bottom, bool left, bool right);

// BONUS SAT functions to be used in waveguide locally hydraulically connected to a half space
// Mixed Neumann and locally absorbing SAT
void asat_neumann_absorbing_top(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, data_t dt, int ixmin, int ixmax, const data_t ** par, const data_t * a);
// Mixed Neumann and Dirichlet SAT
void asat_neumann_dirichlet_top(bool add, const data_t** in, data_t* out, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, const data_t ** par, data_t a, const data_t * m);
void asat_scale_boundaries_bis(data_t** in, int nx, int nz, data_t dx, data_t dz, int ixmin, int ixmax, int izmin, int izmax, const data_t** par, const data_t * a, data_t dt, bool top, bool bottom, bool left, bool right);



/*
Perfectly matched layer following Komatitsch and Tromp 2003, eq. (20) and (21)
for the corners, the equations are modified to account for stretching along both directions

the wavefields have full size nx*nz
the auxiliary arrays have a size n*l where n=nx or nz and l is the PML thickness in grid points
the additional auxiliary arrays have a size l*l

the stretching function is g = (p+1)/(2L)*vmax*log(1/R)(d/L)^p where L is the PML thickness, d=x or z, p is the polynomial power (e.g. =2)
*/

void pml_top(const data_t * in, data_t * out, // input and output wavefields
            const data_t * u_x, const data_t * u_z, // first derivative wavefields
            data_t * &ac, data_t * &an, data_t * s2, // auxilary arrays to store the split fields s1, t, s3, s4, and s2
            const data_t * model, // elastic parameters lamda, mu, rho
            const data_t * g, const data_t * gprime, // stretching function and its derivative
            int nx, int nz, int l, int ixmin, int ixmax, data_t dx, data_t dz, data_t dt, // sizes and sampling
            int version // version for spatial operators (1 or 2)
            );

void pml_bottom(const data_t * in, data_t * out, const data_t * u_x, const data_t * u_z,
            data_t * &ac, data_t * &an, data_t * s2, const data_t * model, const data_t * g, const data_t * gprime,
            int nx, int nz, int l, int ixmin, int ixmax, data_t dx, data_t dz, data_t dt, int version);

void pml_left(const data_t * in, data_t * out, const data_t * u_x, const data_t * u_z,
            data_t * &ac, data_t * &an, data_t * s2,
            data_t * &t2c, data_t * &t2n, data_t * &t3c, data_t * &t3n, data_t * &s5, data_t * &s6, // additional auxiliary arrays for the top left and bottom left corners
            bool pml_T, bool pml_B, // flags to apply or not the 2D PML in the corners. If true izmin must = 0 and izmax must = nz
            const data_t * model, const data_t * g, const data_t * gprime,
            int nx, int nz, int l, int izmin, int izmax, data_t dx, data_t dz, data_t dt, int version);

void pml_right(const data_t * in, data_t * out, const data_t * u_x, const data_t * u_z,
            data_t * &ac, data_t * &an, data_t * s2,
            data_t * &t2c, data_t * &t2n, data_t * &t3c, data_t * &t3n, data_t * &s5, data_t * &s6, // additional auxiliary arrays for the top right and bottom right corners
            bool pml_T, bool pml_B, // flags to apply or not the 2D PML in the corners. If true izmin must = 0 and izmax must = nz
            const data_t * model, const data_t * g, const data_t * gprime,
            int nx, int nz, int l, int izmin, int izmax, data_t dx, data_t dz, data_t dt, int version);

void pml_allocate(data_t * &top_ac, data_t * &top_an, data_t * &top_s2, data_t * &bottom_ac, data_t * &bottom_an, data_t * &bottom_s2,
            data_t * &left_ac, data_t * &left_an, data_t * &left_s2, data_t * &left_t2c, data_t * &left_t2n, data_t * &left_t3c, data_t * &left_t3n, data_t * &left_s5, data_t * &left_s6,
            data_t * &right_ac, data_t * &right_an, data_t * &right_s2, data_t * &right_t2c, data_t * &right_t2n, data_t * &right_t3c, data_t * &right_t3n, data_t * &right_s5, data_t * &right_s6,
            bool pml_T, bool pml_B, bool pml_L, bool pml_R, int nx, int nz, int l);

void pml_deallocate(data_t * &top_ac, data_t * &top_an, data_t * &top_s2, data_t * &bottom_ac, data_t * &bottom_an, data_t * &bottom_s2,
            data_t * &left_ac, data_t * &left_an, data_t * &left_s2, data_t * &left_t2c, data_t * &left_t2n, data_t * &left_t3c, data_t * &left_t3n, data_t * &left_s5, data_t * &left_s6,
            data_t * &right_ac, data_t * &right_an, data_t * &right_s2, data_t * &right_t2c, data_t * &right_t2n, data_t * &right_t3c, data_t * &right_t3n, data_t * &right_s5, data_t * &right_s6,
            bool pml_T, bool pml_B, bool pml_L, bool pml_R, int nx, int nz, int l);