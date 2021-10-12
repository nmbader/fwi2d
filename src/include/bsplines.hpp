#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>

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

namespace OP {

// duplicate z and x slices according to a multiplicity vector for each direction
class duplicate : public oper {
protected: 
    std::vector<int> _mx, _mz; // multiplicity of slices
public:
    duplicate(){}
    ~duplicate(){}
    duplicate(const std::shared_ptr<VEC::hypercube<data_t> > domain, const std::vector<int> &mx, const std::vector<int> &mz){
        assert(domain->getNdim()>=2);
        std::vector<VEC::axis<data_t> > axes = domain->getAxes();
        assert(axes[0].n == mz.size());
        assert(axes[1].n == mx.size());

        size_t nz=0, nx=0;
        for (int i=0; i<mz.size(); i++) {assert(mz[i]>=1) ; nz += mz[i];}
        for (int i=0; i<mx.size(); i++) {assert(mx[i]>=1) ; nx += mx[i];}

        axes[0].n = nz;
        axes[1].n = nx;

        _domain = domain->clone();
        _range = std::make_shared<VEC::hypercube<data_t> >(axes);
        _mx = mx;
        _mz = mz;
    }
    duplicate * clone() const {
        duplicate * op = new duplicate(_domain, _mx,_mz);
        return op;
    }
    void forward(bool add, const std::shared_ptr<VEC::vecReg<data_t> > mod, std::shared_ptr<VEC::vecReg<data_t> > dat) const{
        checkDomainRange(mod,dat);
        forward(add,mod->getVals(),dat->getVals());
    }
    void forward(bool add, data_t * pmod, data_t * pdat) const{
        size_t nz = _domain->getAxis(1).n;
        size_t nz2 = _range->getAxis(1).n;
        size_t nx = _domain->getAxis(2).n;
        size_t nx2 = _range->getAxis(2).n;
        size_t ny = _domain->getN123()/(nx*nz);

        data_t * p = new data_t[ny*nx*nz2];

        size_t count;
        #pragma omp parallel for private(count)
        for (size_t i=0; i<ny*nx; i++){
            size_t count = 0;
            for (size_t iz=0; iz<nz; iz++){
                for (size_t j=count; j<count+_mz[iz]; j++) p[i*nz2+j] = pmod[i*nz+iz];
                count += _mz[iz];
            }
        }
        #pragma omp parallel for private(count)
        for (size_t iy=0; iy<ny; iy++){
            count = 0;
            for (size_t ix=0; ix<nx; ix++){
                for (size_t j=count; j<count+_mx[ix]; j++){
                    for (size_t iz=0; iz<nz2; iz++) pdat[iy*nx2*nz2+j*nz2+iz] = add*pdat[iy*nx2*nz2+j*nz2+iz] + p[iy*nx*nz2+ix*nz2+iz];
                }
                count += _mx[ix];
            }
        }
        delete [] p;
    }
    void adjoint(bool add, std::shared_ptr<VEC::vecReg<data_t> > mod, const std::shared_ptr<VEC::vecReg<data_t> > dat) const{
        checkDomainRange(mod,dat);
        adjoint(add,mod->getVals(),dat->getVals());
    }
    void adjoint(bool add, data_t * pmod, data_t * pdat) const{

        size_t nz = _domain->getAxis(1).n;
        size_t nz2 = _range->getAxis(1).n;
        size_t nx = _domain->getAxis(2).n;
        size_t nx2 = _range->getAxis(2).n;
        size_t ny = _domain->getN123()/(nx*nz);

        if (add == false){
            #pragma omp parallel for
            for (size_t i=0; i<ny*nx*nz; i++) pmod[i] = 0;
        }

        data_t * p = new data_t[ny*nx2*nz];
        size_t count;
        #pragma omp parallel for private(count)
        for (size_t i=0; i<ny*nx2; i++){
            count = 0;
            for (size_t iz=0; iz<nz; iz++){
                p[i*nz+iz] = 0;
                for (size_t j=count; j<count+_mz[iz]; j++) p[i*nz+iz] += pdat[i*nz2+j];
                count += _mz[iz];
            }
        }
        #pragma omp parallel for private(count)
        for (size_t iy=0; iy<ny; iy++){
            count = 0;
            for (size_t ix=0; ix<nx; ix++){
                for (size_t j=count; j<count+_mx[ix]; j++){
                    for (size_t iz=0; iz<nz; iz++) pmod[iy*nx*nz+ix*nz+iz] += p[iy*nx2*nz+j*nz+iz];
                }
                count += _mx[ix];
            }
        }
        delete [] p;
    }
    void adjointn(bool add, std::shared_ptr<VEC::vecReg<data_t> > mod, const std::shared_ptr<VEC::vecReg<data_t> > dat) const{
        checkDomainRange(mod,dat);
        std::shared_ptr<VEC::vecReg<data_t> > temp = mod->clone();
        adjoint(false,temp->getVals(),dat->getVals());
        for (size_t ix=0; ix<_mx.size(); ix++){
            for (size_t iz=0; iz<_mz.size(); iz++) mod->getVals()[ix*_mz.size()+iz] = add*mod->getVals()[ix*_mz.size()+iz] + temp->getVals()[ix*_mz.size()+iz]/(_mx[ix]*_mz[iz]);
        }
    }
};

// Cubic B-splines model preconditioner
// The 1st and last knots are always repeated 4 times
// The control points are assumed to be the same as the knots except the 1st and last
// control points that are repeated only twice
class bsplines3 : public oper {
protected:
    std::vector<data_t> _kx, _kz; // knot vectors
    std::vector<int> _kxmin, _kzmin; // first useful index of the knot vectors
public:
    bsplines3(){}
    ~bsplines3(){}
    bsplines3(const std::shared_ptr<VEC::hypercube<data_t> > domain, const std::shared_ptr<VEC::hypercube<data_t> > range, const std::vector<data_t> &kx, const std::vector<data_t> &kz){
        assert(domain->getAxis(1).n == kz.size()-4);
        assert(domain->getAxis(2).n == kx.size()-4);
        assert(domain->getNdim() == range->getNdim());
        if (domain->getNdim() == 3) assert(domain->getAxis(3).n == range->getAxis(3).n);
        for (int i=0; i<kx.size()-1; i++){
            if (kx[i]>kx[i+1]){
                fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
                throw std::logic_error("Ux knot vector entries must be in ascending order.\n");
            }
        }
        for (int i=0; i<kz.size()-1; i++){
            if (kz[i]>kz[i+1]){
                fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
                throw std::logic_error("Uz knot vector entries must be in ascending order.\n");
            }
        }
        for (int i=1; i<4; i++){
            if ((kx[i] != kx[0]) || (kz[i]!=kz[0]) || (kx[kx.size()-1-i]!=kx[kx.size()-1]) || (kz[kz.size()-1-i]!=kz[kz.size()-1])) {
                fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
                throw std::logic_error("The first and last 4 knots in the knot vectors must be the same.\n");
            }
        }
        assert(kz[0] == range->getAxis(1).o);
        assert(kx[0] == range->getAxis(2).o);
        //assert(kz[kz.size()-1] == range->getAxis(1).o+range->getAxis(1).d*(range->getAxis(1).n-1));
        //assert(kx[kx.size()-1] == range->getAxis(2).o+range->getAxis(2).d*(range->getAxis(2).n-1));

        _domain = domain->clone();
        _range = range->clone();
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

        VEC::axis<data_t> Z = _range->getAxis(1);
        VEC::axis<data_t> X = _range->getAxis(2);
        data_t x,z;
        _kzmin.resize(Z.n,0);
        _kxmin.resize(X.n,0);
        size_t count=0;
        for (size_t i=1; i<Z.n; i++){
            z = Z.o + i*Z.d;
            while (z > _kz[count]) count++;
            _kzmin[i] = count-4;
        }
        count = 0;
        for (size_t i=1; i<X.n; i++){
            x = X.o + i*X.d;
            while (x > _kx[count]) count++;
            _kxmin[i] = count-4;
        }
    }
    bsplines3 * clone() const {
        bsplines3 * op = new bsplines3(_domain, _range, _kx,_kz);
        return op;
    }

    void forward(bool add, const std::shared_ptr<VEC::vecReg<data_t> > mod, std::shared_ptr<VEC::vecReg<data_t> > dat) const{
        checkDomainRange(mod,dat);
        forward(add,mod->getVals(),dat->getVals());
    }
    void forward(bool add, data_t * pmod, data_t * pdat) const{

        VEC::axis<data_t> Z = _range->getAxis(1);
        VEC::axis<data_t> X = _range->getAxis(2);
        size_t ny = _range->getN123()/(X.n*Z.n);
        size_t kx = _kx.size()-4;
        size_t kz = _kz.size()-4;
        data_t x, z;

        for (size_t iy=0; iy<ny; iy++){
            #pragma omp parallel for private(x,z)
            for (size_t ix=0; ix<X.n; ix++){
                x = ix*X.d+X.o;
                for (int iz=0; iz<Z.n; iz++){
                    z = iz*Z.d+Z.o;
                    pdat[iy*X.n*Z.n+ix*Z.n+iz] = add*pdat[iy*X.n*Z.n+ix*Z.n+iz];
                    for (int i=_kxmin[ix]; i<_kxmin[ix]+4; i++){
                        for (int j=_kzmin[iz]; j<_kzmin[iz]+4; j++){
                            pdat[iy*X.n*Z.n+ix*Z.n+iz] += N3(i,x,_kx)*N3(j,z,_kz)*pmod[iy*kx*kz+i*kz+j];
                        }
                    }
                }
            }
        }
    }
    void adjoint(bool add, std::shared_ptr<VEC::vecReg<data_t> > mod, const std::shared_ptr<VEC::vecReg<data_t> > dat) const{
        checkDomainRange(mod,dat);
        adjoint(add,mod->getVals(),dat->getVals());
    }
    void adjoint(bool add, data_t * pmod, data_t * pdat) const{
        VEC::axis<data_t> Z = _range->getAxis(1);
        VEC::axis<data_t> X = _range->getAxis(2);
        size_t ny = _range->getN123()/(X.n*Z.n);
        size_t kx = _kx.size()-4;
        size_t kz = _kz.size()-4;
        data_t x, z;

        if (add == false){
            #pragma omp parallel for
            for (size_t i=0; i<_domain->getN123(); i++) pmod[i] = 0;
        }

        for (size_t iy=0; iy<ny; iy++){
            #pragma omp parallel for private(x,z)
            for (size_t ix=0; ix<X.n; ix++){
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
}

// set the knot vector from control points and multiplicity vectors
void setKnot(std::vector<data_t> &u, std::vector<data_t> &c, std::vector<int> &m){
    assert(c.size() == m.size());
    u.clear();
    u.push_back(c[0]);
    u.push_back(c[0]);
    for (size_t i=0; i<m.size(); i++){
        for (size_t j=0; j<m[i]; j++) u.push_back(c[i]);
    }
    u.push_back(c[c.size()-1]);
    u.push_back(c[c.size()-1]);
}

// populate the coarse spline vector from a dense regular vector using nearest neighbor
void fillin(std::shared_ptr<VEC::vecReg<data_t> > c, const std::shared_ptr<VEC::vecReg<data_t> > v, const std::vector<data_t> &cx, const std::vector<data_t> &cz){
    VEC::axis<data_t> Z = c->getHyper()->getAxis(1);
    VEC::axis<data_t> Z2 = v->getHyper()->getAxis(1);
    VEC::axis<data_t> X = c->getHyper()->getAxis(2);
    VEC::axis<data_t> X2 = v->getHyper()->getAxis(2);
    assert(Z.n == cz.size());
    assert(X.n == cx.size());
    size_t ny = c->getN123()/(X.n*Z.n);
    data_t * pc = c->getVals();
    const data_t * pv = v->getCVals();
    size_t ix2, iz2;
    for (size_t iy=0; iy<ny; iy++){
        for (size_t ix=0; ix<X.n; ix++){
            ix2 = round((cx[ix]-X2.o)/X2.d);
            for (size_t iz=0; iz<Z.n; iz++){
                iz2 = round((cz[iz]-Z2.o)/Z2.d);
                pc[iy*X.n*Z.n+ix*Z.n+iz] = pv[iy*X2.n*Z2.n+ix2*Z2.n+iz2];
            }
        }
    }
}
