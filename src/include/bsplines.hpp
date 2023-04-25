#include "operator.hpp"

#define ZERO 1e-16

// B-splines functions of order 0, 1, 2, 3
// N0
data_t N0(int i, data_t u, const std::vector<data_t> &uk){
    if ((u<uk[i]) || (u>=uk[i+1])) return 0;
    else return 1;
}
// N1
data_t N1(int i, data_t u, const std::vector<data_t> &uk){
    if (std::abs(uk[i+1]-uk[i])<ZERO && std::abs(uk[i+2]-uk[i+1])<ZERO) return 0;
    else if  (std::abs(uk[i+1]-uk[i])<ZERO && std::abs(uk[i+2]-uk[i+1])>=ZERO) return (uk[i+2]-u)/(uk[i+2]-uk[i+1])*N0(i+1,u,uk);
    else if (std::abs(uk[i+1]-uk[i])>=ZERO && std::abs(uk[i+2]-uk[i+1])<ZERO) return (u-uk[i])/(uk[i+1]-uk[i])*N0(i,u,uk);
    else return (u-uk[i])/(uk[i+1]-uk[i])*N0(i,u,uk) + (uk[i+2]-u)/(uk[i+2]-uk[i+1])*N0(i+1,u,uk);
}

// N2
data_t N2(int i, data_t u, const std::vector<data_t> &uk){
    if (std::abs(uk[i+2]-uk[i])<ZERO && std::abs(uk[i+3]-uk[i+1])<ZERO) return 0;
    else if  (std::abs(uk[i+2]-uk[i])<ZERO && std::abs(uk[i+3]-uk[i+1])>=ZERO) return (uk[i+3]-u)/(uk[i+3]-uk[i+1])*N1(i+1,u,uk);
    else if (std::abs(uk[i+2]-uk[i])>=ZERO && std::abs(uk[i+3]-uk[i+1])<ZERO) return (u-uk[i])/(uk[i+2]-uk[i])*N1(i,u,uk);
    else return (u-uk[i])/(uk[i+2]-uk[i])*N1(i,u,uk) + (uk[i+3]-u)/(uk[i+3]-uk[i+1])*N1(i+1,u,uk);
}

// N3
data_t N3(int i, data_t u, const std::vector<data_t> &uk){
    if (std::abs(uk[i+3]-uk[i])<ZERO && std::abs(uk[i+4]-uk[i+1])<ZERO) return 0;
    else if  (std::abs(uk[i+3]-uk[i])<ZERO && std::abs(uk[i+4]-uk[i+1])>=ZERO) return (uk[i+4]-u)/(uk[i+4]-uk[i+1])*N2(i+1,u,uk);
    else if (std::abs(uk[i+3]-uk[i])>=ZERO && std::abs(uk[i+4]-uk[i+1])<ZERO) return (u-uk[i])/(uk[i+3]-uk[i])*N2(i,u,uk);
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

// same as duplicate but in one dimension
class duplicate1d : public loper {
protected: 
    std::vector<int> _mz; // multiplicity vector
public:
    duplicate1d(){}
    ~duplicate1d(){}
    duplicate1d(const hypercube<data_t> &domain, const std::vector<int> &mz){
        std::vector<axis<data_t> > axes = domain.getAxes();
        successCheck(axes[0].n == mz.size(),__FILE__,__LINE__,"Domain first axis must match the size of multiplicity vector mz\n");
        int nz=0;
        for (int i=0; i<mz.size(); i++) {successCheck(mz[i]>=1,__FILE__,__LINE__,"Multiplicity must be >=1\n") ; nz += mz[i];}
        axes[0].n = nz;
        _domain = domain;
        _range = hypercube<data_t>(axes);
        _mz = mz;
    }
    duplicate1d * clone() const {
        duplicate1d * op = new duplicate1d(_domain, _mz);
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
        int ny = _domain.getN123()/(nz);

        int count;
        #pragma omp parallel for private(count)
        for (int i=0; i<ny; i++){
            int count = 0;
            for (int iz=0; iz<nz; iz++){
                for (int j=count; j<count+_mz[iz]; j++) pdat[i*nz2+j] = add*pdat[i*nz2+j] + pmod[i*nz+iz];
                count += _mz[iz];
            }
        }
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

        int nz = _domain.getAxis(1).n;
        int nz2 = _range.getAxis(1).n;
        int ny = _domain.getN123()/(nz);

        if (add == false) memset(pmod, 0, _domain.getN123()*sizeof(data_t));

        int count;
        #pragma omp parallel for private(count)
        for (int i=0; i<ny; i++){
            count = 0;
            for (int iz=0; iz<nz; iz++){
                for (int j=count; j<count+_mz[iz]; j++) pmod[i*nz+iz] += pdat[i*nz2+j];
                count += _mz[iz];
            }
        }
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
        successCheck(domain.getNdim()>=2,__FILE__,__LINE__,"Domain must have at least 2 dimensions\n");
        successCheck(domain.getAxis(1).n == kz.size()-4,__FILE__,__LINE__,"Domain first axis must have 4 samples less than the z knot vector\n");
        successCheck(domain.getAxis(2).n == kx.size()-4,__FILE__,__LINE__,"Domain second axis must have 4 samples less than the x knot vector\n");
        successCheck(domain.getNdim() == range.getNdim(),__FILE__,__LINE__,"Domain and range must have the same number of dimensions\n");
        for (int i=3; i<=domain.getNdim(); i++) successCheck(domain.getAxis(i).n == range.getAxis(i).n,__FILE__,__LINE__,"Domain and range must have the same size in the other dimensions\n");
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

        _kx[kx.size()-1] += range.getAxis(2).d*1e-03;
        _kx[kx.size()-2] += range.getAxis(2).d*1e-03;
        _kx[kx.size()-3] += range.getAxis(2).d*1e-03;
        _kx[kx.size()-4] += range.getAxis(2).d*1e-03;
        _kz[kz.size()-1] += range.getAxis(1).d*1e-03;
        _kz[kz.size()-2] += range.getAxis(1).d*1e-03;
        _kz[kz.size()-3] += range.getAxis(1).d*1e-03;
        _kz[kz.size()-4] += range.getAxis(1).d*1e-03;

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

        #pragma omp parallel for private(x,z)
        for (int iy=0; iy<ny; iy++){
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

// same as bsplines3 but in one dimension 
class bsplines31d : public loper {
protected:
    std::vector<data_t> _kz; // knot vector
    std::vector<int> _kzmin; // first useful index of the knot vector
public:
    bsplines31d(){}
    ~bsplines31d(){}
    bsplines31d(const hypercube<data_t> &domain, const hypercube<data_t> &range, const std::vector<data_t> &kz){
        successCheck(domain.getAxis(1).n == kz.size()-4,__FILE__,__LINE__,"Domain first axis must have 4 samples less than the z knot vector\n");
        successCheck(domain.getNdim() == range.getNdim(),__FILE__,__LINE__,"Domain and range must have the same number of dimensions\n");
        for (int i=2; i<=domain.getNdim(); i++) successCheck(domain.getAxis(i).n == range.getAxis(i).n,__FILE__,__LINE__,"Domain and range must have the same size in the other dimensions\n");
        for (int i=0; i<kz.size()-1; i++){
            successCheck(kz[i]<=kz[i+1],__FILE__,__LINE__,"Uz knot vector entries must be in ascending order.\n");
        }
        for (int i=1; i<4; i++){
            successCheck( (kz[i]==kz[0]) && (kz[kz.size()-1-i]==kz[kz.size()-1]),__FILE__,__LINE__,"The first and last 4 knots in the knot vector must be the same.\n");
        }
        successCheck(kz[0] == range.getAxis(1).o,__FILE__,__LINE__,"The origin of z knot and range z axis must be the same\n");

        _domain = domain;
        _range = range;
        _kz = kz;

        _kz[kz.size()-1] += range.getAxis(1).d*1e-03;
        _kz[kz.size()-2] += range.getAxis(1).d*1e-03;
        _kz[kz.size()-3] += range.getAxis(1).d*1e-03;
        _kz[kz.size()-4] += range.getAxis(1).d*1e-03;

        axis<data_t> Z = _range.getAxis(1);
        data_t z;
        _kzmin.resize(Z.n,0);
        int count=0;
        for (int i=1; i<Z.n; i++){
            z = Z.o + i*Z.d;
            while (z > _kz[count]) count++;
            _kzmin[i] = count-4;
        }
    }
    bsplines31d * clone() const {
        bsplines31d * op = new bsplines31d(_domain, _range, _kz);
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
        int ny = _range.getN123()/(Z.n);
        int kz = _kz.size()-4;
        data_t z;

        if (add == false) memset(pdat, 0, _range.getN123()*sizeof(data_t));
        #pragma omp parallel for private(z)
        for (int iy=0; iy<ny; iy++){
            for (int iz=0; iz<Z.n; iz++){
                z = iz*Z.d+Z.o;
                for (int j=_kzmin[iz]; j<_kzmin[iz]+4; j++){
                    pdat[iy*Z.n+iz] += N3(j,z,_kz)*pmod[iy*kz+j];
                }
            }
        }
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){
        
        axis<data_t> Z = _range.getAxis(1);
        int ny = _range.getN123()/(Z.n);
        int kz = _kz.size()-4;
        data_t z;

        if (add == false) memset(pmod, 0, _domain.getN123()*sizeof(data_t));

        for (int iy=0; iy<ny; iy++){
            for (int iz=0; iz<Z.n; iz++){
                z = iz*Z.d+Z.o;
                for (int j=_kzmin[iz]; j<_kzmin[iz]+4; j++){
                    pmod[iy*kz+j] += N3(j,z,_kz)*pdat[iy*Z.n+iz];
                }
            }
        }
    }
};

// operator doing the opposite (not the inverse) of bsplines3; it populates the B-spline model using a gridded model with bilinear interpolation
class bsfillin : public loper {
protected:
    std::vector<data_t> _controlx; // control points
    std::vector<data_t> _controlz; // control points
public:
    bsfillin(){}
    ~bsfillin(){}
    bsfillin(const hypercube<data_t> &domain, const std::vector<data_t> &controlx, const std::vector<data_t> &controlz){

        _domain = domain;
        std::vector<axis<data_t> > axes = domain.getAxes();
        axes[0]=controlz.size();
        axes[1]=controlx.size();
        _range = hypercube<data_t>(axes);
        _controlx = controlx;
        _controlz = controlz;
    }
    bsfillin * clone() const {
        bsfillin * op = new bsfillin(_domain, _controlx, _controlz);
        return op;
    }
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){

        axis<data_t> Z = _range.getAxis(1);
        axis<data_t> Z2 = _domain.getAxis(1);
        axis<data_t> X = _range.getAxis(2);
        axis<data_t> X2 = _domain.getAxis(2);
        int ny = _range.getN123()/(X.n*Z.n);
        data_t (*pc) [X.n][Z.n] = (data_t (*) [X.n][Z.n]) pdat;
        const data_t (*pv) [X2.n][Z2.n] = (const data_t (*) [X2.n][Z2.n]) pmod;
        for (int iy=0; iy<ny; iy++){
            #pragma omp parallel for
            for (int ix=0; ix<X.n; ix++){
                int ix2 = floor((_controlx[ix]-X2.o)/X2.d);
                int ix3 = std::min(X2.n-1,ix2+1);
                data_t wx = (_controlx[ix] - X2.o - ix2*X2.d)/X2.d;
                for (int iz=0; iz<Z.n; iz++){
                    int iz2 = floor((_controlz[iz]-Z2.o)/Z2.d);
                    int iz3 = std::min(Z2.n-1,iz2+1);
                    data_t wz = (_controlz[iz] - Z2.o - iz2*Z2.d)/Z2.d;
                    pc[iy][ix][iz] = add*pc[iy][ix][iz] + (1-wx)*(1-wz)*pv[iy][ix2][iz2]  + wx*(1-wz)*pv[iy][ix3][iz2] + (1-wx)*wz*pv[iy][ix2][iz3] + wx*wz*pv[iy][ix3][iz3];
                }
            }
        }
    }
    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

        if (!add) memset(pmod, 0, _domain.getN123()*sizeof(data_t));
        
        axis<data_t> Z = _range.getAxis(1);
        axis<data_t> Z2 = _domain.getAxis(1);
        axis<data_t> X = _range.getAxis(2);
        axis<data_t> X2 = _domain.getAxis(2);
        int ny = _range.getN123()/(X.n*Z.n);
        const data_t (*pc) [X.n][Z.n] = (const data_t (*) [X.n][Z.n]) pdat;
        data_t (*pv) [X2.n][Z2.n] = (data_t (*) [X2.n][Z2.n]) pmod;
        #pragma omp parallel for
        for (int iy=0; iy<ny; iy++){
            for (int ix=0; ix<X.n; ix++){
                int ix2 = floor((_controlx[ix]-X2.o)/X2.d);
                int ix3 = std::min(X2.n-1,ix2+1);
                data_t wx = (_controlx[ix] - X2.o - ix2*X2.d)/X2.d;
                for (int iz=0; iz<Z.n; iz++){
                    int iz2 = floor((_controlz[iz]-Z2.o)/Z2.d);
                    int iz3 = std::min(Z2.n-1,iz2+1);
                    data_t wz = (_controlz[iz] - Z2.o - iz2*Z2.d)/Z2.d;                    
                    pv[iy][ix2][iz2] += (1-wx)*(1-wz)*pc[iy][ix][iz];
                    pv[iy][ix3][iz2] += wx*(1-wz)*pc[iy][ix][iz];
                    pv[iy][ix2][iz3] += (1-wx)*wz*pc[iy][ix][iz];
                    pv[iy][ix3][iz3] += wx*wz*pc[iy][ix][iz];
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

// populate the coarse spline vector from a dense regular vector using linear interpolation
void fillin(std::shared_ptr<vecReg<data_t> > c, const std::shared_ptr<vecReg<data_t> > v, const std::vector<data_t> &cx, const std::vector<data_t> &cz){
    axis<data_t> Z = c->getHyper()->getAxis(1);
    axis<data_t> Z2 = v->getHyper()->getAxis(1);
    axis<data_t> X = c->getHyper()->getAxis(2);
    axis<data_t> X2 = v->getHyper()->getAxis(2);
    successCheck(Z.n == cz.size(),__FILE__,__LINE__,"The coarse spline vector must have the same z size as the control vector\n");
    successCheck(X.n == cx.size(),__FILE__,__LINE__,"The coarse spline vector must have the same x size as the control vector\n");
    int ny = c->getN123()/(X.n*Z.n);
    data_t (*pc) [X.n][Z.n] = (data_t (*) [X.n][Z.n]) c->getVals();
    const data_t (*pv) [X2.n][Z2.n] = (const data_t (*) [X2.n][Z2.n]) v->getCVals();
    int ix2, iz2, ix3, iz3;
    data_t wx, wz;
    for (int iy=0; iy<ny; iy++){
        for (int ix=0; ix<X.n; ix++){
            ix2 = floor((cx[ix]-X2.o)/X2.d);
            ix3 = std::min(X2.n-1,ix2+1);
            wx = (cx[ix] - X2.o - ix2*X2.d)/X2.d;
            for (int iz=0; iz<Z.n; iz++){
                iz2 = floor((cz[iz]-Z2.o)/Z2.d);
                iz3 = std::min(Z2.n-1,iz2+1);
                wz = (cz[iz] - Z2.o - iz2*Z2.d)/Z2.d;
                pc[iy][ix][iz] = (1-wx)*(1-wz)*pv[iy][ix2][iz2]  + wx*(1-wz)*pv[iy][ix3][iz2] + (1-wx)*wz*pv[iy][ix2][iz3] + wx*wz*pv[iy][ix3][iz3];
            }
        }
    }
}

// same as fillin but in one dimension
void fillin1d(std::shared_ptr<vecReg<data_t> > c, const std::shared_ptr<vecReg<data_t> > v, const std::vector<data_t> &cz){
    axis<data_t> Z = c->getHyper()->getAxis(1);
    axis<data_t> Z2 = v->getHyper()->getAxis(1);
    successCheck(Z.n == cz.size(),__FILE__,__LINE__,"The coarse spline vector must have the same z size as the control vector\n");
    int ny = c->getN123()/(Z.n);
    data_t * pc = c->getVals();
    const data_t * pv = v->getCVals();
    int iz2, iz3;
    data_t wz;
    for (int iy=0; iy<ny; iy++){
        for (int iz=0; iz<Z.n; iz++){
            iz2 = floor((cz[iz]-Z2.o)/Z2.d);
            iz3 = std::min(Z2.n-1,iz2+1);
            wz = (cz[iz] - Z2.o - iz2*Z2.d)/Z2.d;
            pc[iy*Z.n+iz] = (1-wz)*pv[iy*Z2.n+iz2] + wz*pv[iy*Z2.n+iz3];
        }
    }
}

#undef ZERO