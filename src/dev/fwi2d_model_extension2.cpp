#include <string.h>
#include "we_op.hpp"
#include "bsplines.hpp"
#include "nlsolver.hpp"
#include "IO.hpp"
#include "seplib.h"

typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;

// operator that maps from several models (one for each source) to the physical space with damping (cosine squared) around the sources location 
// m = W.ms where ms is the extended model across sources, 'W' is the row operator of normalized cosine squared weighting operator for each source 's'
class sWeighting : public loper {
protected:
    data_t _dpower; // damping strength for cosine^2
    data_t _dwidthx; // damping width (radius)
    data_t _dwidthz;
    std::shared_ptr<vec> _nW; // normalization for the weighting (sum of all weighting operators)
    std::vector<std::vector<data_t> > _sxz; // sources location
    std::vector<std::vector<int> > _box; // box delimiting the size of the zone of influence of each shot point

public:
    sWeighting(){}
    ~sWeighting(){}
    sWeighting(const hypercube<data_t> &domain,const std::vector<std::vector<data_t> > &sxz, const std::vector<std::vector<std::vector<data_t> > > &rxz, const data_t dwidthx, const data_t dwidthz, const data_t dpower, const data_t xextension, const data_t zextension){
        std::vector<ax > axes = domain.getAxes();
        successCheck(axes.size()>=3,__FILE__,__LINE__,"The domain must contain at least 3 axes\n");
        int ns=domain.getN123()/(axes[0].n*axes[1].n*axes[2].n);
        ax S(ns,0,1);
        successCheck(ns==sxz.size(),__FILE__,__LINE__,"The domain size must be consistent with the number of sources\n");
        
        _domain = hyper(axes[0],axes[1],axes[2],S);
        _range = hyper(axes[0],axes[1],axes[2]);
        _dwidthx = dwidthx;
        _dwidthz = dwidthz;
        _dpower = dpower;
        _sxz=sxz;

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        _nW = std::make_shared<vec> (vec(hyper(Z,X)));
        _nW->zero();
        data_t (*pw) [Z.n] = (data_t (*)[Z.n]) _nW->getVals();

        data_t val0 = pow(cos(_dpower*0.5*M_PI),2);

        for (int s=0; s<ns; s++){

            // find the bounding box for each shot
            int ixmin = ceil((sxz[s][0]-dwidthx-X.o)/X.d);
            int ixmax = ceil((sxz[s][0]+dwidthx-X.o)/X.d);
            int izmin = ceil((sxz[s][1]-dwidthz-Z.o)/Z.d);
            int izmax = ceil((sxz[s][1]+dwidthz-Z.o)/Z.d);
            ixmin=std::max(0,ixmin);
            ixmax=std::min(X.n,ixmax);
            izmin=std::max(0,izmin);
            izmax=std::min(Z.n,izmax);

            // find the bounding box for each shot+receivers spread (shot zone of influence)
            std::vector<int> v (4,0);
            data_t rxmin=rxz[s][0][0];
            data_t rxmax=rxz[s][0][0];
            data_t rzmin=rxz[s][0][1];
            data_t rzmax=rxz[s][0][1];
            for (int r=0; r<rxz[s].size(); r++){
                if (rxmin>rxz[s][r][0]) rxmin=rxz[s][r][0];
                if (rxmax<rxz[s][r][0]) rxmax=rxz[s][r][0];
                if (rzmin>rxz[s][r][1]) rzmin=rxz[s][r][1];
                if (rzmax<rxz[s][r][1]) rzmax=rxz[s][r][1];
            }
            rxmin = std::min(rxmin,_sxz[s][0]-dwidthx);
            rxmax = std::max(rxmax,_sxz[s][0]+dwidthx);
            rzmin = std::min(rzmin,_sxz[s][1]-dwidthz);
            rzmax = std::max(rzmax,_sxz[s][1]+dwidthz);
            rxmin -= xextension; 
            rxmax += xextension;
            rzmin -= zextension; 
            rzmax += zextension; 
            v[0] = floor((rxmin-X.o)/X.d);
            v[1] = ceil((rxmax-X.o)/X.d);
            v[2] = floor((rzmin-Z.o)/Z.d);
            v[3] = ceil((rzmax-Z.o)/Z.d);
            v[0]=std::max(0,v[0]);
            v[1]=std::min(X.n,v[1]);
            v[2]=std::max(0,v[2]);
            v[3]=std::min(Z.n,v[3]);
            
            if (xextension<0) {v[0]=0; v[1]=X.n;}
            if (zextension<0) {v[2]=0; v[3]=Z.n;}

            _box.push_back(v);

            int ixw=ceil(dwidthx/X.d);
            int izw=ceil(dwidthz/Z.d);

            // normalization weight
            for (int ix=0; ix<v[0]-ixw; ix++){
                data_t val0x=val0;
                for (int iz=0; iz<v[2]-izw; iz++) pw[ix][iz] += val0x*val0;
                for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pw[ix][iz] += val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2);
                for (int iz=v[2]; iz<v[3]; iz++) pw[ix][iz] += val0x;
                for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pw[ix][iz] += val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2);
                for (int iz=v[3]+izw; iz<Z.n; iz++) pw[ix][iz] += val0x*val0;
            }
            for (int ix=std::max(v[0]-ixw,0); ix<v[0]; ix++){
                data_t val0x=pow(cos(_dpower*0.5*M_PI*(v[0]-ix)/ixw),2);
                for (int iz=0; iz<v[2]-izw; iz++) pw[ix][iz] += val0x*val0;
                for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pw[ix][iz] += val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2);
                for (int iz=v[2]; iz<v[3]; iz++) pw[ix][iz] += val0x;
                for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pw[ix][iz] += val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2);
                for (int iz=v[3]+izw; iz<Z.n; iz++) pw[ix][iz] += val0x*val0;
            }
            for (int ix=v[0]; ix<v[1]; ix++){
                data_t val0x=1;
                for (int iz=0; iz<v[2]-izw; iz++) pw[ix][iz] += val0x*val0;
                for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pw[ix][iz] += val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2);
                for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pw[ix][iz] += val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2);
                for (int iz=v[3]+izw; iz<Z.n; iz++) pw[ix][iz] += val0x*val0;
            }
            for (int ix=v[1]; ix<std::min(v[1]+ixw,X.n); ix++){
                data_t val0x=pow(cos(_dpower*0.5*M_PI*(v[1]-ix)/ixw),2);
                for (int iz=0; iz<v[2]-izw; iz++) pw[ix][iz] += val0x*val0;
                for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pw[ix][iz] += val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2);
                for (int iz=v[2]; iz<v[3]; iz++) pw[ix][iz] += val0x;
                for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pw[ix][iz] += val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2);
                for (int iz=v[3]+izw; iz<Z.n; iz++) pw[ix][iz] += val0x*val0;
            }
            for (int ix=v[1]+ixw; ix<X.n; ix++){
                data_t val0x=val0;
                for (int iz=0; iz<v[2]-izw; iz++) pw[ix][iz] += val0x*val0;
                for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pw[ix][iz] += val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2);
                for (int iz=v[2]; iz<v[3]; iz++) pw[ix][iz] += val0x;
                for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pw[ix][iz] += val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2);
                for (int iz=v[3]+izw; iz<Z.n; iz++) pw[ix][iz] += val0x*val0;
            }

            for (int ix=ixmin; ix<ixmax; ix++){
                data_t x = X.o + ix*X.d;
                data_t rx = (x-_sxz[s][0])/_dwidthx;
                rx*=rx;
                for (int iz=izmin; iz<izmax; iz++){
                    data_t z = Z.o + iz*Z.d;
                    data_t rz = (z-_sxz[s][1])/_dwidthz;
                    rz*=rz;
                    data_t d = sqrt(rx+rz);
                    data_t val = 1;
                    if (d<1) val = cos(_dpower*0.5*M_PI*(1-d));
                    pw[ix][iz] += val*val;
                }
            }

            for (int ix=v[0]; ix<ixmin; ix++){
                for (int iz=v[2]; iz<v[3]; iz++) pw[ix][iz] += 1;
            }
            for (int ix=ixmax; ix<v[1]; ix++){
                for (int iz=v[2]; iz<v[3]; iz++) pw[ix][iz] += 1;
            }
            for (int ix=ixmin; ix<ixmax; ix++){
                for (int iz=v[2]; iz<izmin; iz++) pw[ix][iz] += 1;
                for (int iz=izmax; iz<v[3]; iz++) pw[ix][iz] += 1;
            }
        }
    }

    sWeighting * clone() const {
        //sWeighting * op = new sWeighting(_domain, _sxz, _dwidthx, _dwidthz, _dpower, _xextension, _zextension);
        sWeighting * op = new sWeighting();
        op->_domain=_domain;
        op->_range=_range;
        op->_dpower=_dpower;
        op->_dwidthx=_dwidthx;
        op->_dwidthz=_dwidthz;
        op->_sxz=_sxz;
        op->_box=_box;
        op->_nW = _nW->clone();
        return op;
    }    
    int getNs(){return _sxz.size();}
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){

        if (!add) memset(pdat,0,_range.getN123()*sizeof(data_t));

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        int ns = _sxz.size();
        const data_t (*pm) [C.n][X.n][Z.n] = (const data_t (*)[C.n][X.n][Z.n]) pmod;
        data_t (*pd) [X.n][Z.n] = (data_t (*)[X.n][Z.n]) pdat;
        const data_t (*pw) [Z.n] = (const data_t (*)[Z.n]) _nW->getVals();

        data_t val0 = pow(cos(_dpower*0.5*M_PI),2);

        for (int s=0; s<ns; s++){
            
            int ixmin = ceil((_sxz[s][0]-_dwidthx-X.o)/X.d);
            int ixmax = ceil((_sxz[s][0]+_dwidthx-X.o)/X.d);
            int izmin = ceil((_sxz[s][1]-_dwidthz-Z.o)/Z.d);
            int izmax = ceil((_sxz[s][1]+_dwidthz-Z.o)/Z.d);
            ixmin=std::max(0,ixmin);
            ixmax=std::min(X.n,ixmax);
            izmin=std::max(0,izmin);
            izmax=std::min(Z.n,izmax);

            std::vector<int> v = _box[s];
            int ixw=ceil(_dwidthx/X.d);
            int izw=ceil(_dwidthz/Z.d);
            
            for (int c=0; c<C.n; c++){

                for (int ix=0; ix<v[0]-ixw; ix++){
                    data_t val0x=val0;
                    for (int iz=0; iz<v[2]-izw; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*val0/pw[ix][iz];
                    for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[2]; iz<v[3]; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x/pw[ix][iz];
                    for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[3]+izw; iz<Z.n; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*val0/pw[ix][iz];
                }
                for (int ix=std::max(v[0]-ixw,0); ix<v[0]; ix++){
                    data_t val0x=pow(cos(_dpower*0.5*M_PI*(v[0]-ix)/ixw),2);
                    for (int iz=0; iz<v[2]-izw; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*val0/pw[ix][iz];
                    for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[2]; iz<v[3]; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x/pw[ix][iz];
                    for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[3]+izw; iz<Z.n; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*val0/pw[ix][iz];
                }
                for (int ix=v[0]; ix<v[1]; ix++){
                    data_t val0x=1;
                    for (int iz=0; iz<v[2]-izw; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*val0/pw[ix][iz];
                    for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[3]+izw; iz<Z.n; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*val0/pw[ix][iz];
                }
                for (int ix=v[1]; ix<std::min(v[1]+ixw,X.n); ix++){
                    data_t val0x=pow(cos(_dpower*0.5*M_PI*(v[1]-ix)/ixw),2);
                    for (int iz=0; iz<v[2]-izw; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*val0/pw[ix][iz];
                    for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[2]; iz<v[3]; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x/pw[ix][iz];
                    for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[3]+izw; iz<Z.n; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*val0/pw[ix][iz];
                }
                for (int ix=v[1]+ixw; ix<X.n; ix++){
                    data_t val0x=val0;
                    for (int iz=0; iz<v[2]-izw; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*val0/pw[ix][iz];
                    for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[2]; iz<v[3]; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x/pw[ix][iz];
                    for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[3]+izw; iz<Z.n; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]*val0x*val0/pw[ix][iz];
                }

                #pragma omp parallel for
                for (int ix=ixmin; ix<ixmax; ix++){
                    data_t x = X.o + ix*X.d;
                    data_t rx = (x-_sxz[s][0])/_dwidthx;
                    rx*=rx;
                    for (int iz=izmin; iz<izmax; iz++){
                        data_t z = Z.o + iz*Z.d;
                        data_t rz = (z-_sxz[s][1])/_dwidthz;
                        rz*=rz;
                        data_t d = sqrt(rx+rz);
                        data_t val = 1;
                        if (d<1) val = cos(_dpower*0.5*M_PI*(1-d));
                        pd[c][ix][iz] += pm[s][c][ix][iz]*val*val/pw[ix][iz];
                    }
                }

                #pragma omp parallel for
                for (int ix=v[0]; ix<ixmin; ix++){
                    for (int iz=v[2]; iz<v[3]; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]/pw[ix][iz];
                }
                #pragma omp parallel for
                for (int ix=ixmax; ix<v[1]; ix++){
                    for (int iz=v[2]; iz<v[3]; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]/pw[ix][iz];
                }
                #pragma omp parallel for
                for (int ix=ixmin; ix<ixmax; ix++){
                    for (int iz=v[2]; iz<izmin; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]/pw[ix][iz];
                    for (int iz=izmax; iz<v[3]; iz++) pd[c][ix][iz] += pm[s][c][ix][iz]/pw[ix][iz];
                }
            }            
        }
    }

    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        int ns = _sxz.size();
        data_t (*pm) [C.n][X.n][Z.n] = (data_t (*)[C.n][X.n][Z.n]) pmod;
        const data_t (*pd) [X.n][Z.n] = (const data_t (*)[X.n][Z.n]) pdat;
        const data_t (*pw) [Z.n] = (const data_t (*)[Z.n]) _nW->getVals();

        data_t val0 = pow(cos(_dpower*0.5*M_PI),2);

        for (int s=0; s<ns; s++){
            int ixmin = ceil((_sxz[s][0]-_dwidthx-X.o)/X.d);
            int ixmax = ceil((_sxz[s][0]+_dwidthx-X.o)/X.d);
            int izmin = ceil((_sxz[s][1]-_dwidthz-Z.o)/Z.d);
            int izmax = ceil((_sxz[s][1]+_dwidthz-Z.o)/Z.d);
            ixmin=std::max(0,ixmin);
            ixmax=std::min(X.n,ixmax);
            izmin=std::max(0,izmin);
            izmax=std::min(Z.n,izmax);

            std::vector<int> v = _box[s];
            int ixw=ceil(_dwidthx/X.d);
            int izw=ceil(_dwidthz/Z.d);

            for (int c=0; c<C.n; c++){

                for (int ix=0; ix<v[0]-ixw; ix++){
                    data_t val0x=val0;
                    for (int iz=0; iz<v[2]-izw; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*val0/pw[ix][iz];
                    for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[2]; iz<v[3]; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x/pw[ix][iz];
                    for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[3]+izw; iz<Z.n; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*val0/pw[ix][iz];
                }
                for (int ix=std::max(v[0]-ixw,0); ix<v[0]; ix++){
                    data_t val0x=pow(cos(_dpower*0.5*M_PI*(v[0]-ix)/ixw),2);
                    for (int iz=0; iz<v[2]-izw; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*val0/pw[ix][iz];
                    for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[2]; iz<v[3]; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x/pw[ix][iz];
                    for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[3]+izw; iz<Z.n; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*val0/pw[ix][iz];
                }
                for (int ix=v[0]; ix<v[1]; ix++){
                    data_t val0x=1;
                    for (int iz=0; iz<v[2]-izw; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*val0/pw[ix][iz];
                    for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[3]+izw; iz<Z.n; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*val0/pw[ix][iz];
                }
                for (int ix=v[1]; ix<std::min(v[1]+ixw,X.n); ix++){
                    data_t val0x=pow(cos(_dpower*0.5*M_PI*(v[1]-ix)/ixw),2);
                    for (int iz=0; iz<v[2]-izw; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*val0/pw[ix][iz];
                    for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[2]; iz<v[3]; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x/pw[ix][iz];
                    for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[3]+izw; iz<Z.n; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*val0/pw[ix][iz];
                }
                for (int ix=v[1]+ixw; ix<X.n; ix++){
                    data_t val0x=val0;
                    for (int iz=0; iz<v[2]-izw; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*val0/pw[ix][iz];
                    for (int iz=std::max(v[2]-izw,0); iz<v[2]; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[2]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[2]; iz<v[3]; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x/pw[ix][iz];
                    for (int iz=v[3]; iz<std::min(v[3]+izw,Z.n); iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*pow(cos(_dpower*0.5*M_PI*(v[3]-iz)/izw),2)/pw[ix][iz];
                    for (int iz=v[3]+izw; iz<Z.n; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val0x*val0/pw[ix][iz];
                }

                for (int ix=ixmin; ix<ixmax; ix++){
                    data_t x = X.o + ix*X.d;
                    data_t rx = (x-_sxz[s][0])/_dwidthx;
                    rx*=rx;
                    for (int iz=izmin; iz<izmax; iz++){
                        data_t z = Z.o + iz*Z.d;
                        data_t rz = (z-_sxz[s][1])/_dwidthz;
                        rz*=rz;
                        data_t d = sqrt(rx+rz);
                        data_t val = 1;
                        if (d<1) val = cos(_dpower*0.5*M_PI*(1-d));
                        pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]*val*val/pw[ix][iz];
                    }
                }
                for (int ix=v[0]; ix<ixmin; ix++){
                    for (int iz=v[2]; iz<v[3]; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]/pw[ix][iz];
                }
                for (int ix=ixmax; ix<v[1]; ix++){
                    for (int iz=v[2]; iz<v[3]; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]/pw[ix][iz];
                }
                for (int ix=ixmin; ix<ixmax; ix++){
                    for (int iz=v[2]; iz<izmin; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]/pw[ix][iz];
                    for (int iz=izmax; iz<v[3]; iz++) pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[c][ix][iz]/pw[ix][iz];
                }
            }
        }
    }
};

// operator that maps from several models (one for each source) to the physical space with damping (cosine squared) around the sources location, and then spray the difference (initial - collapsed) across sources 
// Op = (I-S.W).ms where ms is the extended model across sources, 'W' is the row operator of normalized cosine squared weighting operator for each source 's', 'S' is a spraying operator
class model_extension : public loper {
protected:
    sWeighting * _sW;

public:
    model_extension(){}
    ~model_extension(){delete _sW;}
    model_extension(sWeighting * sW){
        _sW = sW->clone();
        _domain = *sW->getDomain();
        _range = _domain;
    }
    model_extension * clone() const {
        model_extension * op = new model_extension(_sW);
        return op;
    }
    sWeighting * getsW(){return _sW;}
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        int ns = _sW->getNs();
        const data_t (*pm) [C.n][X.n][Z.n] = (const data_t (*)[C.n][X.n][Z.n]) pmod;
        data_t (*pd) [C.n][X.n][Z.n] = (data_t (*)[C.n][X.n][Z.n]) pdat;

        data_t * vec = new data_t[C.n*X.n*Z.n];
        memset(vec, 0, C.n*X.n*Z.n*sizeof(data_t));
        data_t (*pv) [X.n][Z.n] = (data_t (*)[X.n][Z.n]) vec;

        // W.m
        _sW->apply_forward(false,pmod,vec);

        // (I-S.W).m
        for (int s=0; s<ns; s++){
            for (int c=0; c<C.n; c++){
                #pragma omp parallel for
                for (int ix=0; ix<X.n;ix++){
                    for (int iz=0; iz<Z.n;iz++) pd[s][c][ix][iz] = add*pd[s][c][ix][iz] + pm[s][c][ix][iz] - pv[c][ix][iz];
                }
            }
        }
        delete [] vec;
    }

    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        int ns = _sW->getNs();
        data_t (*pm) [C.n][X.n][Z.n] = (data_t (*)[C.n][X.n][Z.n]) pmod;
        const data_t (*pd) [C.n][X.n][Z.n] = (const data_t (*)[C.n][X.n][Z.n]) pdat;

        data_t * vec = new data_t[C.n*X.n*Z.n];
        memset(vec, 0, C.n*X.n*Z.n*sizeof(data_t));
        data_t (*pv) [X.n][Z.n] = (data_t (*)[X.n][Z.n]) vec;

        // -S'.d and I.d
        for (int s=0; s<ns; s++){
            for (int c=0; c<C.n; c++){
                #pragma omp parallel for
                for (int ix=0; ix<X.n;ix++){
                    for (int iz=0; iz<Z.n;iz++) {
                        pv[c][ix][iz] -= pd[s][c][ix][iz];
                        pm[s][c][ix][iz] = add*pm[s][c][ix][iz] + pd[s][c][ix][iz];
                    }
                }
            }
        }

        // (I-W'.S').d
        _sW->apply_adjoint(true,pmod,vec);
        delete [] vec;
    }
};

// soft clip of the elastic model and the Vs/Vp ratio for the extended model
class emodelSoftClipExt : public nloper {
    emodelSoftClip * _S;
public:
    emodelSoftClipExt(){}
    ~emodelSoftClipExt(){delete _S;}
    emodelSoftClipExt(const hypercube<data_t> &domain, emodelSoftClip * S){
        successCheck((domain.getNdim()>=3) && (domain.getAxis(3).n>=3),__FILE__,__LINE__,"The domain must be at least 3D with 3rd dimension containing at least 3 fields\n");
        _domain = domain;
        _range = domain;
        hyper hyp(domain.getAxis(1),domain.getAxis(2),domain.getAxis(3));
        _S = S->clone();
    }
    emodelSoftClipExt * clone() const {
        emodelSoftClipExt * op = new emodelSoftClipExt(_domain, _S);
        return op;
    }   
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){
        int n=_domain.getAxis(1).n*_domain.getAxis(2).n*_domain.getAxis(3).n;
        int ns=_domain.getN123()/n;
        for (int s=0; s<ns; s++) _S->apply_forward(add,pmod+s*n,pdat+s*n);
    }
    void apply_jacobianT(bool add, data_t * pmod, const data_t * pmod0, const data_t * pdat){
        int n=_domain.getAxis(1).n*_domain.getAxis(2).n*_domain.getAxis(3).n;
        int ns=_domain.getN123()/n;
        for (int s=0; s<ns; s++) _S->apply_jacobianT(add,pmod+s*n,pmod0+s*n,pdat+s*n);
    }
};

// Executable to run 2D FWI

int main(int argc, char **argv){

int rank=0, size=0;
#ifdef ENABLE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    fprintf (stderr,"\n====================\nSize of MPI communicator = %d ; current rank = %d\n====================\n",size,rank);
#endif
    
    initpar(argc,argv);

// Read parameters for wave propagation and inversion
    param par;
    readParameters(argc, argv, par);
    int verbose=par.verbose;
    if (rank>0) par.verbose=0;
    par.device+=rank;

// Set the maximum number of threads
    if (par.nthreads>0) omp_set_num_threads(par.nthreads);

// Read inputs/outputs files
    std::string source_file="none", model_file="none", data_file="none", output_file="none", ioutput_file="none", obj_func_file="none";
    readParam<std::string>(argc, argv, "source", source_file);
    readParam<std::string>(argc, argv, "model", model_file);
    readParam<std::string>(argc, argv, "data", data_file);
    readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "ioutput", ioutput_file);
    readParam<std::string>(argc, argv, "obj_func", obj_func_file);

// Read some parameters
    data_t dwidthx=0, dwidthz=0, dpower=0, xextension=-1, zextension=-1;
    readParam<data_t>(argc, argv, "damping_widthx", dwidthx);
    readParam<data_t>(argc, argv, "damping_widthz", dwidthz);
    readParam<data_t>(argc, argv, "damping_power", dpower);
    readParam<data_t>(argc, argv, "xextension", xextension);
    readParam<data_t>(argc, argv, "zextension", zextension);

    successCheck(source_file!="none",__FILE__,__LINE__,"Source wavelet is not provided\n");
    successCheck(model_file!="none",__FILE__,__LINE__,"Earth model is not provided\n");
    successCheck(data_file!="none",__FILE__,__LINE__,"Data to be inverted is not provided\n");

    std::shared_ptr<vec> src = read<data_t>(source_file, par.format);
    std::shared_ptr<vec> data = read<data_t>(data_file, par.format);
    std::shared_ptr<vec> model = read<data_t>(model_file, par.format);
    hyper hyp0 = *model->getHyper();

    std::shared_ptr<vec> gmask = nullptr;
    std::shared_ptr<vec> w = nullptr;
    std::shared_ptr<vec> filter = nullptr;
    std::shared_ptr<vec> invDiagH = nullptr;
    std::shared_ptr<vec> prior = nullptr;
    if (par.mask_file!="none") {gmask = read<data_t>(par.mask_file, par.format); successCheck(gmask->getN123()==model->getN123(),__FILE__,__LINE__,"Gradient mask must have the same number of samples as the model\n");}
    if (par.weights_file!="none") {w = read<data_t>(par.weights_file, par.format); successCheck(w->getN123()==data->getN123(),__FILE__,__LINE__,"Data weights must have the same number of samples as the data\n");}
    if (par.filter_file!="none") {filter = read<data_t>(par.filter_file, par.format); successCheck(filter->getHyper()->getAxis(1).d==data->getHyper()->getAxis(1).d,__FILE__,__LINE__,"Filter and data must have the same sampling rate\n");}
    if (par.inverse_diagonal_hessian_file!="none") {invDiagH = read<data_t>(par.inverse_diagonal_hessian_file, par.format); successCheck(invDiagH->getN123()==model->getN123(),__FILE__,__LINE__,"Inverse diagonal Hessian must have the same number of samples as the model\n");}
    if (par.prior_file!="none") {prior = read<data_t>(par.prior_file, par.format); successCheck(prior->getN123()==model->getN123(),__FILE__,__LINE__,"Prior model must have the same number of samples as the model\n");}

// Analyze the inputs and parameters and modify if necessary
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0);
    analyzeBsplines(*model->getHyper(),par);
    analyzeNLInversion(par);
    par.sextension=true;

// ----------------------------------------------------------------------------------------//
// Extend model along sources, extend prior, mask and inverse Hessian if necessary
// ----------------------------------------------------------------------------------------//
    std::vector<ax > axes = model->getHyper()->getAxes();
    int n=model->getN123();
    int nxz=n/par.nmodels;
    data_t vs0 = model->sum(nxz, 2*nxz);
    data_t rho0 = model->sum(2*nxz, 3*nxz);
    vs0/=nxz;
    rho0/=nxz;
{
    ax S(par.ns,0,1);
    axes.push_back(S);
    std::shared_ptr<vec> tmp = std::make_shared<vec>(vec(hyper(axes)));
    for (int s=0; s<par.ns; s++) memcpy(tmp->getVals()+s*n,model->getVals(),model->getN123()*sizeof(data_t));
    model=tmp;
}
    if (par.mask_file != "none"){
        if (gmask->getN123() == n){
            std::shared_ptr<vec> tmp = std::make_shared<vec>(vec(hyper(axes)));
            for (int s=0; s<par.ns; s++) memcpy(tmp->getVals()+s*n,gmask->getVals(),gmask->getN123()*sizeof(data_t));
            gmask=tmp;
        }
        else successCheck(gmask->getN123()==n*par.ns,__FILE__,__LINE__,"The extended gradient mask has an incorrect size\n");
    }
    if (par.inverse_diagonal_hessian_file != "none"){
        if (invDiagH->getN123() == n){
            std::shared_ptr<vec> tmp = std::make_shared<vec>(vec(hyper(axes)));
            for (int s=0; s<par.ns; s++) memcpy(tmp->getVals()+s*n,invDiagH->getVals(),invDiagH->getN123()*sizeof(data_t));
            invDiagH=tmp;
        }
        else successCheck(gmask->getN123()==n*par.ns,__FILE__,__LINE__,"The extended inverse Hessian has an incorrect size\n");
    }
    if (par.prior_file != "none"){
        if (prior->getN123() == n){
            std::shared_ptr<vec> tmp = std::make_shared<vec>(vec(hyper(axes)));
            for (int s=0; s<par.ns; s++) memcpy(tmp->getVals()+s*n,prior->getVals(),prior->getN123()*sizeof(data_t));
            prior=tmp;
        }
        else successCheck(prior->getN123()==n*par.ns,__FILE__,__LINE__,"The prior model has an incorrect size\n");
    }

    sWeighting * sW = new sWeighting(*model->getHyper(), par.sxz, par.rxz, dwidthx, dwidthz, dpower, xextension, zextension);
    model_extension * E = new model_extension(sW);
    delete sW;
// ----------------------------------------------------------------------------------------//

// Build model parameterization precon
// ----------------------------------------------------------------------------------------//
    std::shared_ptr<vec> model_temp = model;
    std::shared_ptr<vec> prior_temp = prior;

    nloper * P = nullptr;
    if (par.model_parameterization==0) P = new lam_mu_rho(*model->getHyper());
    else if (par.model_parameterization==2) P = new ip_is_rho(*model->getHyper());
    else if (par.model_parameterization==3) {
        P = new vs_vpvs_rho(*model->getHyper(), vs0, rho0);
        if (par.verbose>0) fprintf(stderr,"Average Vs and rho used in model parameterization are %.3f and %.3f\n",vs0/nxz, rho0/nxz);
    }
    if (P != nullptr){
        model_temp = std::make_shared<vec>(*model->getHyper());
        model_temp->zero();
        P->inverse(false,model_temp,model);
        if (par.prior_file != "none") {
            prior_temp = std::make_shared<vec>(*model->getHyper());
            prior_temp->zero();
            P->inverse(false,prior_temp,prior);
        }
    }
    model = model_temp;
    prior = prior_temp;
// ----------------------------------------------------------------------------------------//
  
// ----------------------------------------------------------------------------------------//
// Build model precon for model extension along sources with B-splines included
// ----------------------------------------------------------------------------------------//
    std::shared_ptr<vec> bsmodel = model;
    std::shared_ptr<vec> bsmask = gmask;
    std::shared_ptr<vec> bsinvDiagH = invDiagH;
    std::shared_ptr<vec> bsprior = prior;
    loper * BD;

if (par.bsplines)
{
    std::vector<data_t> kx;
    std::vector<data_t> kz;
    setKnot(kx,par.bs_controlx,par.bs_mx);
    setKnot(kz,par.bs_controlz,par.bs_mz);

    bsfillin F(*model->getHyper(),par.bs_controlx,par.bs_controlz);
    bsmodel = std::make_shared<vec>(vec(*F.getRange()));
    bsmodel->zero();
    F.apply_forward(false, model->getVals(), bsmodel->getVals());
    if (par.mask_file != "none") {
        bsmask = std::make_shared<vec>(vec(*F.getRange()));
        bsmask->zero();
        F.apply_forward(false, gmask->getVals(), bsmask->getVals());
    }
    if (par.inverse_diagonal_hessian_file != "none") {
        bsinvDiagH = std::make_shared<vec>(vec(*F.getRange()));
        bsinvDiagH->zero();
        F.apply_forward(false, invDiagH->getVals(), bsinvDiagH->getVals());
    }
    if (prior != nullptr) {
        bsprior = std::make_shared<vec>(vec(*F.getRange()));
        bsprior->zero();
        F.apply_forward(false, prior->getVals(), bsprior->getVals());
    }
   
    duplicate D(*F.getRange(),par.bs_mx,par.bs_mz);
    bsplines3 B(*D.getRange(),*model->getHyper(),kx,kz);
    BD = new chainLOper(&B,&D);
}
 
// ----------------------------------------------------------------------------------------//
// Build model precon if soft clipping is activated
// ----------------------------------------------------------------------------------------//
    emodelSoftClipExt * S;
    if (par.soft_clip) 
    {
        if (par.model_parameterization==0) {
            data_t mu_min = par.rhomin*par.vsmin*par.vsmin;
            data_t mu_max = par.rhomax*par.vsmax*par.vsmax;
            data_t lam_min = par.rhomin*(par.vpmin*par.vpmin - 2*par.vsmax*par.vsmax);
            data_t lam_max = par.rhomax*(par.vpmax*par.vpmax - 2*par.vsmin*par.vsmin);
            lam_min = std::max((data_t)0.0 , lam_min);
            emodelSoftClip S0(hyp0, lam_min, lam_max, mu_min, mu_max, par.rhomin, par.rhomax, 1, 9, 9);
            S = new emodelSoftClipExt(*model->getHyper(), &S0);
        }
        else if (par.model_parameterization==2){
            data_t ip_min = par.rhomin*par.vpmin;
            data_t ip_max = par.rhomax*par.vpmax;
            data_t is_min = par.rhomin*par.vsmin;
            data_t is_max = par.rhomax*par.vsmax;
            emodelSoftClip S0(hyp0, ip_min, ip_max, is_min, is_max, par.rhomin, par.rhomax, 1/sqrt(2.00001), 9, 9);
            S = new emodelSoftClipExt(*model->getHyper(), &S0);
        }
        else if (par.model_parameterization==3){
            emodelSoftClip S0(hyp0, log(par.vsmin/vs0), log(par.vsmax/vs0), -10.0, log(par.vpmax/par.vsmin - sqrt(2)), log(par.rhomin/rho0), log(par.rhomax/rho0), 1, 9, 9);
            S = new emodelSoftClipExt(*model->getHyper(), &S0);
        }
        else {
            emodelSoftClip S0(hyp0, par.vpmin, par.vpmax, par.vsmin, par.vsmax, par.rhomin, par.rhomax, 1/sqrt(2.00001), 9, 9);
            S = new emodelSoftClipExt(*model->getHyper(), &S0);
        }
    }
// ----------------------------------------------------------------------------------------//

    if (rank>0) par.verbose=verbose;
    nloper * op = nullptr;
    nl_we_op * L;
    if (par.nmodels==2) L=new nl_we_op_a(*model->getHyper(),allsrc,par);
    else if (par.nmodels==3 && !par.acoustic_elastic) L=new nl_we_op_e(*model->getHyper(),allsrc,par);
    else if (par.nmodels==3 && par.acoustic_elastic) L=new nl_we_op_ae(*model->getHyper(),allsrc,par);
    else if (par.nmodels==5) L=new nl_we_op_vti(*model->getHyper(),allsrc,par);

    if (rank>0) par.verbose=0;

    if (par.bsplines)
    {
        if (P != nullptr)
        {
            if (par.soft_clip) {
                chainNLOper SBD(S,BD);
                op  = new chainNLOper(P,&SBD);
            } 
            else op = new chainNLOper(P,BD);
        }
        else
        {
            if (par.soft_clip) op  = new chainNLOper(S,BD);
            else op = BD->clone();
        }
    }
    else
    {
        if (P != nullptr){
            if (par.soft_clip) op = new chainNLOper(P,S);
            else op = P->clone();
        }
        else{
            if (par.soft_clip) op = S->clone();
        }
    }
    if (par.bsplines) delete BD;
    if (P != nullptr) delete P;
    if (par.soft_clip) delete S;

    nloper * D = nullptr;
    
    if (op==nullptr) D = E->clone();
    else D = new chainNLOper(E, op);

    nlls_fwi_reg * prob = new nlls_fwi_reg(L, D, bsmodel, data, par.lambda, bsprior, op, bsmask, w, filter);

    lsearch * ls;
    if (par.lsearch=="weak_wolfe") ls = new weak_wolfe(par.ls_c1, par.ls_a0, par.ls_a1, par.ls_version);
    else if(par.lsearch=="strong_wolfe") ls = new strong_wolfe(par.ls_c1, par.ls_c2, par.ls_a0, par.ls_a1, par.ls_max_step, par.ls_version);
    else ls = new regular_wolfe(par.ls_c1, par.ls_c2, par.ls_a0, par.ls_a1, par.ls_max_step, par.ls_version);

    nlsolver * solver;
    if (par.nlsolver=="nlsd") solver = new nlsd(par.niter, par.max_trial, par.threshold, ls); 
    else if (par.nlsolver=="nlcg") solver = new nlcg(par.niter, par.max_trial, par.threshold, ls); 
    else if (par.nlsolver=="bfgs") solver = new bfgs(par.niter, par.max_trial, par.threshold, ls); 
    else solver = new lbfgs(par.niter, par.max_trial, par.threshold, ls, bsinvDiagH, par.lbfgs_m); 

    solver->run(prob, par.verbose>0, ioutput_file, par.isave, par.format, par.datapath);
    model = bsmodel;
    
    if (D != nullptr) delete D;
    if (op != nullptr)
    {
        std::shared_ptr<vec> tmp = std::make_shared<vec>(*op->getRange());
        tmp->zero();
        op->forward(false, model, tmp);
        model = tmp;
        delete op;
    }
    
    std::shared_ptr<vec> fmodel = std::make_shared<vec>(*E->getsW()->getRange());
    fmodel->zero();
    E->getsW()->forward(false,model,fmodel);

    if (rank==0 && output_file!="none") {
        write<data_t>(fmodel, output_file, par.format, par.datapath);
        write<data_t>(model, output_file+".ext", par.format, par.datapath);
    }
    if (rank==0 && obj_func_file!="none") {
        std::shared_ptr<vec > func = std::make_shared<vec > (hyper(solver->_func.size()));
        memcpy(func->getVals(), solver->_func.data(), solver->_func.size()*sizeof(data_t));
        write(func,obj_func_file, par.format, par.datapath);
        if (D != nullptr){
            std::shared_ptr<vec > dfunc = std::make_shared<vec > (hyper(prob->_dfunc.size()));
            std::shared_ptr<vec > mfunc = std::make_shared<vec > (hyper(prob->_mfunc.size()));
            memcpy(dfunc->getVals(), prob->_dfunc.data(), prob->_dfunc.size()*sizeof(data_t));
            memcpy(mfunc->getVals(), prob->_mfunc.data(), prob->_mfunc.size()*sizeof(data_t));
            write(dfunc,obj_func_file+".d", par.format, par.datapath);
            write(mfunc,obj_func_file+".m", par.format, par.datapath);
        }
    }

    delete E;
    delete L;
    delete prob;

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif

return 0;
}