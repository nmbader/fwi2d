#include "operator.hpp"

// B-splines functions of order 0, 1, 2, 3
// N0
data_t N0(int i, data_t u, const std::vector<data_t> &uk){
    if ((u<uk[i]) || (u>=uk[i+1])) return 0;
    else return 1;
}
// N1
data_t N1(int i, data_t u, const std::vector<data_t> &uk){
    if ((uk[i+1]-uk[i]==0) && (uk[i+2]-uk[i+1]==0)) return 0;
    else if  ((uk[i+1]-uk[i]==0) && (uk[i+2]-uk[i+1]!=0)) return (uk[i+2]-u)/(uk[i+2]-uk[i+1])*N0(i+1,u,uk);
    else if ((uk[i+1]-uk[i]!=0) && (uk[i+2]-uk[i+1]==0)) return (u-uk[i])/(uk[i+1]-uk[i])*N0(i,u,uk);
    else return (u-uk[i])/(uk[i+1]-uk[i])*N0(i,u,uk) + (uk[i+2]-u)/(uk[i+2]-uk[i+1])*N0(i+1,u,uk);
}

// N2
data_t N2(int i, data_t u, const std::vector<data_t> &uk){
    if ((uk[i+2]-uk[i]==0) && (uk[i+3]-uk[i+1]==0)) return 0;
    else if  ((uk[i+2]-uk[i]==0) && (uk[i+3]-uk[i+1]!=0))return (uk[i+3]-u)/(uk[i+3]-uk[i+1])*N1(i+1,u,uk);
    else if ((uk[i+2]-uk[i]!=0) && (uk[i+3]-uk[i+1]==0)) return (u-uk[i])/(uk[i+2]-uk[i])*N1(i,u,uk);
    else return (u-uk[i])/(uk[i+2]-uk[i])*N1(i,u,uk) + (uk[i+3]-u)/(uk[i+3]-uk[i+1])*N1(i+1,u,uk);
}

// N3
data_t N3(int i, data_t u, const std::vector<data_t> &uk){
    if ((uk[i+3]-uk[i]==0) && (uk[i+4]-uk[i+1]==0)) return 0;
    else if  ((uk[i+3]-uk[i]==0) && (uk[i+4]-uk[i+1]!=0)) return (uk[i+4]-u)/(uk[i+4]-uk[i+1])*N2(i+1,u,uk);
    else if ((uk[i+3]-uk[i]!=0) && (uk[i+4]-uk[i+1]==0)) return (u-uk[i])/(uk[i+3]-uk[i])*N2(i,u,uk);
    else return (u-uk[i])/(uk[i+3]-uk[i])*N2(i,u,uk) + (uk[i+4]-u)/(uk[i+4]-uk[i+1])*N2(i+1,u,uk);
}

// duplicate z and x slices according to a multiplicity vector for each direction
class duplicate : public loper {
protected: 
    std::vector<int> _mx, _mz; // multiplicity of slices
public:
    duplicate(){}
    ~duplicate(){}
    duplicate(const hypercube<data_t> &domain, const std::vector<int> &mx, const std::vector<int> &mz){
        successCheck(domain.getNdim()>=2,__FILE__,__LINE__,"Domain for duplicate operator must have at least 2 dimensions\n");
        std::vector<axis<data_t> > axes = domain.getAxes();
        successCheck(axes[0].n == mz.size(),__FILE__,__LINE__,"Domain first axis must match the size of multiplicity vector mz\n");
        successCheck(axes[1].n == mx.size(),__FILE__,__LINE__,"Domain second axis must match the size of multiplicity vector mx\n");

        int nz=0, nx=0;
        for (int i=0; i<mz.size(); i++) {successCheck(mz[i]>=1,__FILE__,__LINE__,"Multiplicity must be >=1\n") ; nz += mz[i];}
        for (int i=0; i<mx.size(); i++) {successCheck(mx[i]>=1,__FILE__,__LINE__,"Multiplicity must be >=1\n") ; nx += mx[i];}

        axes[0].n = nz;
        axes[1].n = nx;

        _domain = domain;
        _range = hypercube<data_t>(axes);
        _mx = mx;
        _mz = mz;
    }
    duplicate * clone() const {
        duplicate * op = new duplicate(_domain, _mx,_mz);
        return op;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        return checkCompatible(mod, dat);
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain.isCompatible(_domain) && range.isCompatible(_range),__FILE__,__LINE__,"Domain or range and incompatible with the operator\n");
        _domain = domain;
        _range = range;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat) {
        
        int nz = _domain.getAxis(1).n;
        int nz2 = _range.getAxis(1).n;
        int nx = _domain.getAxis(2).n;
        int nx2 = _range.getAxis(2).n;
        int ny = _domain.getN123()/(nx*nz);

        data_t * p = new data_t[ny*nx*nz2];

        int count;
        #pragma omp parallel for private(count)
        for (int i=0; i<ny*nx; i++){
            int count = 0;
            for (int iz=0; iz<nz; iz++){
                for (int j=count; j<count+_mz[iz]; j++) p[i*nz2+j] = pmod[i*nz+iz];
                count += _mz[iz];
            }
        }
        #pragma omp parallel for private(count)
        for (int iy=0; iy<ny; iy++){
            count = 0;
            for (int ix=0; ix<nx; ix++){
                for (int j=count; j<count+_mx[ix]; j++){
                    for (int iz=0; iz<nz2; iz++) pdat[iy*nx2*nz2+j*nz2+iz] = add*pdat[iy*nx2*nz2+j*nz2+iz] + p[iy*nx*nz2+ix*nz2+iz];
                }
                count += _mx[ix];
            }
        }
        delete [] p;
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

        int nz = _domain.getAxis(1).n;
        int nz2 = _range.getAxis(1).n;
        int nx = _domain.getAxis(2).n;
        int nx2 = _range.getAxis(2).n;
        int ny = _domain.getN123()/(nx*nz);

        if (add == false) memset(pmod, 0, _domain.getN123()*sizeof(data_t));

        data_t * p = new data_t[ny*nx2*nz];
        int count;
        #pragma omp parallel for private(count)
        for (int i=0; i<ny*nx2; i++){
            count = 0;
            for (int iz=0; iz<nz; iz++){
                p[i*nz+iz] = 0;
                for (int j=count; j<count+_mz[iz]; j++) p[i*nz+iz] += pdat[i*nz2+j];
                count += _mz[iz];
            }
        }
        #pragma omp parallel for private(count)
        for (int iy=0; iy<ny; iy++){
            count = 0;
            for (int ix=0; ix<nx; ix++){
                for (int j=count; j<count+_mx[ix]; j++){
                    for (int iz=0; iz<nz; iz++) pmod[iy*nx*nz+ix*nz+iz] += p[iy*nx2*nz+j*nz+iz];
                }
                count += _mx[ix];
            }
        }
        delete [] p;
    }
};

// Cubic B-splines model preconditioner
// The 1st and last knots are always repeated 4 times
// The control points are assumed to be the same as the knots except the 1st and last
// control points that are repeated only twice
class bsplines3 : public loper {
protected:
    std::vector<data_t> _kx, _kz; // knot vectors
    std::vector<int> _kxmin, _kzmin; // first useful index of the knot vectors
public:
    bsplines3(){}
    ~bsplines3(){}
    bsplines3(const hypercube<data_t> &domain, const hypercube<data_t> &range, const std::vector<data_t> &kx, const std::vector<data_t> &kz){
        successCheck(domain.getAxis(1).n == kz.size()-4,__FILE__,__LINE__,"Domain first axis must have 4 samples less than the z knot vector\n");
        successCheck(domain.getAxis(2).n == kx.size()-4,__FILE__,__LINE__,"Domain second axis must have 4 samples less than the x knot vector\n");
        successCheck(domain.getNdim() == range.getNdim(),__FILE__,__LINE__,"Domain and range must have the same number of dimensions\n");
        if (domain.getNdim() == 3) successCheck(domain.getAxis(3).n == range.getAxis(3).n,__FILE__,__LINE__,"Domain and range must have the same size in the 3rd dimension\n");
        for (int i=0; i<kx.size()-1; i++){
            successCheck(kx[i]<=kx[i+1],__FILE__,__LINE__,"Ux knot vector entries must be in ascending order.\n");
        }
        for (int i=0; i<kz.size()-1; i++){
            successCheck(kz[i]<=kz[i+1],__FILE__,__LINE__,"Uz knot vector entries must be in ascending order.\n");
        }
        for (int i=1; i<4; i++){
            successCheck((kx[i] == kx[0]) && (kz[i]==kz[0]) && (kx[kx.size()-1-i]==kx[kx.size()-1]) && (kz[kz.size()-1-i]==kz[kz.size()-1]),__FILE__,__LINE__,"The first and last 4 knots in the knot vectors must be the same.\n");
        }
        successCheck(kz[0] == range.getAxis(1).o,__FILE__,__LINE__,"The origin of z knot and range z axis must be the same\n");
        successCheck(kx[0] == range.getAxis(2).o,__FILE__,__LINE__,"The origin of x knot and range x axis must be the same\n");

        _domain = domain;
        _range = range;
        _kx = kx;
        _kz = kz;

        _kx[kx.size()-1] += 1e-06;
        _kx[kx.size()-2] += 1e-06;
        _kx[kx.size()-3] += 1e-06;
        _kx[kx.size()-4] += 1e-06;
        _kz[kz.size()-1] += 1e-06;
        _kz[kz.size()-2] += 1e-06;
        _kz[kz.size()-3] += 1e-06;
        _kz[kz.size()-4] += 1e-06;

        axis<data_t> Z = _range.getAxis(1);
        axis<data_t> X = _range.getAxis(2);
        data_t x,z;
        _kzmin.resize(Z.n,0);
        _kxmin.resize(X.n,0);
        int count=0;
        for (int i=1; i<Z.n; i++){
            z = Z.o + i*Z.d;
            while (z > _kz[count]) count++;
            _kzmin[i] = count-4;
        }
        count = 0;
        for (int i=1; i<X.n; i++){
            x = X.o + i*X.d;
            while (x > _kx[count]) count++;
            _kxmin[i] = count-4;
        }
    }
    bsplines3 * clone() const {
        bsplines3 * op = new bsplines3(_domain, _range, _kx,_kz);
        return op;
    }
    bool checkDomainRange(const std::shared_ptr<vecReg<data_t> > mod, const std::shared_ptr<vecReg<data_t> > dat) const {
        return checkSame(mod, dat);
    }
    void setDomainRange(const hypercube<data_t> &domain, const hypercube<data_t> &range){
        successCheck(domain==_domain && range==_range,__FILE__,__LINE__,"Domain or range and different than the operator\n");
        _domain = domain;
        _range = range;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){

        axis<data_t> Z = _range.getAxis(1);
        axis<data_t> X = _range.getAxis(2);
        int ny = _range.getN123()/(X.n*Z.n);
        int kx = _kx.size()-4;
        int kz = _kz.size()-4;
        data_t x, z;

        if (add == false) memset(pdat, 0, _range.getN123()*sizeof(data_t));

        for (int iy=0; iy<ny; iy++){
            #pragma omp parallel for private(x,z)
            for (int ix=0; ix<X.n; ix++){
                x = ix*X.d+X.o;
                for (int iz=0; iz<Z.n; iz++){
                    z = iz*Z.d+Z.o;
                    //pdat[iy*X.n*Z.n+ix*Z.n+iz] = add*pdat[iy*X.n*Z.n+ix*Z.n+iz];
                    for (int i=_kxmin[ix]; i<_kxmin[ix]+4; i++){
                        for (int j=_kzmin[iz]; j<_kzmin[iz]+4; j++){
                            pdat[iy*X.n*Z.n+ix*Z.n+iz] += N3(i,x,_kx)*N3(j,z,_kz)*pmod[iy*kx*kz+i*kz+j];
                        }
                    }
                }
            }
        }
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        
        axis<data_t> Z = _range.getAxis(1);
        axis<data_t> X = _range.getAxis(2);
        int ny = _range.getN123()/(X.n*Z.n);
        int kx = _kx.size()-4;
        int kz = _kz.size()-4;
        data_t x, z;

        if (add == false) memset(pmod, 0, _domain.getN123()*sizeof(data_t));

        for (int iy=0; iy<ny; iy++){
            //#pragma omp parallel for private(x,z)
            for (int ix=0; ix<X.n; ix++){
                x = ix*X.d+X.o;
                for (int iz=0; iz<Z.n; iz++){
                    z = iz*Z.d+Z.o;
                    for (int i=_kxmin[ix]; i<_kxmin[ix]+4; i++){
                        for (int j=_kzmin[iz]; j<_kzmin[iz]+4; j++){
                            pmod[iy*kx*kz+i*kz+j] += N3(i,x,_kx)*N3(j,z,_kz)*pdat[iy*X.n*Z.n+ix*Z.n+iz];
                        }
                    }
                }
            }
        }
    }
};

// set the knot vector from control points and multiplicity vectors
void setKnot(std::vector<data_t> &u, std::vector<data_t> &c, std::vector<int> &m){
    successCheck(c.size() == m.size(),__FILE__,__LINE__,"Control and multiplicity vectors must have the same size\n");
    u.clear();
    u.push_back(c[0]);
    u.push_back(c[0]);
    for (int i=0; i<m.size(); i++){
        for (int j=0; j<m[i]; j++) u.push_back(c[i]);
    }
    u.push_back(c[c.size()-1]);
    u.push_back(c[c.size()-1]);
}

// populate the coarse spline vector from a dense regular vector using nearest neighbor
void fillin(std::shared_ptr<vecReg<data_t> > c, const std::shared_ptr<vecReg<data_t> > v, const std::vector<data_t> &cx, const std::vector<data_t> &cz){
    axis<data_t> Z = c->getHyper()->getAxis(1);
    axis<data_t> Z2 = v->getHyper()->getAxis(1);
    axis<data_t> X = c->getHyper()->getAxis(2);
    axis<data_t> X2 = v->getHyper()->getAxis(2);
    successCheck(Z.n == cz.size(),__FILE__,__LINE__,"The coarse spline vector must have the same z size as the control vector\n");
    successCheck(X.n == cx.size(),__FILE__,__LINE__,"The coarse spline vector must have the same x size as the control vector\n");
    int ny = c->getN123()/(X.n*Z.n);
    data_t * pc = c->getVals();
    const data_t * pv = v->getCVals();
    int ix2, iz2;
    for (int iy=0; iy<ny; iy++){
        for (int ix=0; ix<X.n; ix++){
            ix2 = round((cx[ix]-X2.o)/X2.d);
            for (int iz=0; iz<Z.n; iz++){
                iz2 = round((cz[iz]-Z2.o)/Z2.d);
                pc[iy*X.n*Z.n+ix*Z.n+iz] = pv[iy*X2.n*Z2.n+ix2*Z2.n+iz2];
            }
        }
    }
}
