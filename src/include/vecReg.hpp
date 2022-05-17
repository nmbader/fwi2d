#pragma once

#include <string>
#include <vector>
#include <memory>
#include <cstring>
#include <stdlib.h>
#include <ctime>
#include <time.h>
#include <cmath>
#include <complex>
#include <omp.h>

#ifdef DOUBLE_PRECISION
    typedef double data_t;
#else
    typedef float data_t;
#endif
#ifdef DOUBLE_PRECISION
#undef DOUBLE_PRECISION
#endif

#ifdef ENABLE_MPI
    #include "mpi.h"
#endif

template <typename T1>
class axis {
public:
	int n;
	T1 o, d;
	std::string label, unit;

	axis() {}  // Default constructor does nothing
	~axis() {} // Destructor

	axis(const int n1, T1 o1=0.0, T1 d1=1.0, const std::string label1="", const std::string unit1=""){
		n=n1; o=o1; d=d1; label=label1; unit=unit1;
	}
	axis (const axis &ax){
		n=ax.n; o=ax.o; d=ax.d; label=ax.label; unit=ax.unit;
	}
	axis &operator=(const axis &ax) {
		n=ax.n; o=ax.o; d=ax.d; label=ax.label; unit=ax.unit;
		return *this;
	}

	bool operator==(const axis &ax) const{
		bool answer=false;
		if ((n==ax.n) && (o==ax.o) && (d==ax.d)) answer=true;
		return answer;
	}
    bool operator!=(const axis &ax) const{
		bool answer=false;
		if ((n==ax.n) && (o==ax.o) && (d==ax.d)) answer=true;
		return !answer;
	}
	void print() const{
		fprintf(stderr,"n=%d ; o=%0.3f ; d=%0.3f\n",n,o,d);
	}
};

template <typename T1>
class hypercube {
private:
	long _n123;
	std::vector<axis<T1> > _axes;
public:
	hypercube() {}
	~hypercube() {}
	hypercube(const long n1){
		_axes.emplace_back(axis<T1>(n1));
		_n123=n1;
	}
	hypercube(const long n1, const long n2){
		_axes.emplace_back(axis<T1>(n1));
		_axes.emplace_back(axis<T1>(n2));
		_n123=n1*n2;
	}
	hypercube(const long n1, const long n2, const long n3){
		_axes.emplace_back(axis<T1>(n1));
		_axes.emplace_back(axis<T1>(n2));
		_axes.emplace_back(axis<T1>(n3));
		_n123=n1*n2*n3;
	}
	hypercube(const long n1, const long n2, const long n3, const long n4){
		_axes.emplace_back(axis<T1>(n1));
		_axes.emplace_back(axis<T1>(n2));
		_axes.emplace_back(axis<T1>(n3));
		_axes.emplace_back(axis<T1>(n4));
		_n123=n1*n2*n3*n4;
	}
    hypercube(const long n1, const long n2, const long n3, const long n4, const long n5){
		_axes.emplace_back(axis<T1>(n1));
		_axes.emplace_back(axis<T1>(n2));
		_axes.emplace_back(axis<T1>(n3));
		_axes.emplace_back(axis<T1>(n4));
		_axes.emplace_back(axis<T1>(n5));
		_n123=n1*n2*n3*n4*n5;
	}
	hypercube(const axis<T1> &a1) {
		_axes.push_back(a1);
		_n123 = (long)a1.n;
	}
	hypercube(const axis<T1> &a1, const axis<T1> &a2) {
		_axes.push_back(a1);
		_axes.push_back(a2);
		_n123 = (long)a1.n*(long)a2.n;
	}
	hypercube(const axis<T1> &a1, const axis<T1> &a2, const axis<T1> &a3) {
		_axes.push_back(a1);
		_axes.push_back(a2);
		_axes.push_back(a3);
		_n123 = (long)a1.n*(long)a2.n*(long)a3.n;
	}
	hypercube(const axis<T1> &a1, const axis<T1> &a2, const axis<T1> &a3, const axis<T1> &a4) {
		_axes.push_back(a1);
		_axes.push_back(a2);
		_axes.push_back(a3);
		_axes.push_back(a4);
		_n123 = (long)a1.n*(long)a2.n*(long)a3.n*(long)a4.n;
	}
    hypercube(const axis<T1> &a1, const axis<T1> &a2, const axis<T1> &a3, const axis<T1> &a4, const axis<T1> &a5) {
		_axes.push_back(a1);
		_axes.push_back(a2);
		_axes.push_back(a3);
		_axes.push_back(a4);
		_axes.push_back(a5);
		_n123 = (long)a1.n*(long)a2.n*(long)a3.n*(long)a4.n*(long)a5.n;
	}
	hypercube(const std::vector<axis<T1> > &axes){
		_axes = axes;
		_n123=1;
		for (int i=0; i<axes.size(); i++) _n123 *= axes[i].n;
	}
	hypercube(const hypercube<T1> &hyper){
		_axes = hyper._axes;
		_n123 = hyper._n123;
	}
    void operator=(const hypercube<T1> &hyper){
        _axes = hyper._axes;
        _n123 = hyper._n123;
    }
	bool operator==(const hypercube<T1> &hyper) const{
		bool answer=true;
		if (_n123 != hyper._n123) answer = false;
        else if (_axes.size() != hyper._axes.size()) answer = false;
		else{
			for (int i=0; i<_axes.size(); i++){
				if ((_axes[i].n!=hyper._axes[i].n) || (_axes[i].d!=hyper._axes[i].d) || (_axes[i].o!=hyper._axes[i].o)){
					answer=false;
					break;
				}
			}
		}
		return answer;
	}
    bool operator!=(const hypercube<T1> &hyper) const{
		bool answer=true;
		if (_n123 != hyper._n123) answer = false;
        else if (_axes.size() != hyper._axes.size()) answer = false;
		else{
			for (int i=0; i<_axes.size(); i++){
				if ((_axes[i].n!=hyper._axes[i].n) || (_axes[i].d!=hyper._axes[i].d) || (_axes[i].o!=hyper._axes[i].o)){
					answer=false;
					break;
				}
			}
		}
		return !answer;
	}
	bool isCompatible(const hypercube<T1> &hyper) const{
		bool answer=true;
		if (_n123 != hyper._n123) answer = false;
        else if (_axes.size() != hyper._axes.size()) answer = false;
		else{
			for (int i=0; i<_axes.size(); i++){
				if (_axes[i].n!=hyper._axes[i].n){
					answer=false;
					break;
				}
			}
		}
		return answer;
	}
	void setAxes(const std::vector<axis<T1> > &axes){
		_axes=axes;
		_n123=1;
		for (int i=0; i<axes.size(); i++) _n123 *= axes[i].n;
	}
	const std::vector<axis<T1> > getAxes() const {return _axes;}
	const axis<T1> getAxis(int n) const {return _axes[n-1];}
	long getN123() const { return _n123; }
	int getNdim() const { return (int)_axes.size();}
	void print() const{
        fprintf(stderr,"\n");
		for (int i=0; i<_axes.size(); i++){
			fprintf(stderr,"Axis %d: n=%d ; o=%0.3f ; d=%0.3f\n",i+1,_axes[i].n,_axes[i].o,_axes[i].d);
		}
        fprintf(stderr,"\n");
	}	
};

template<typename T1>
static T1 static_pow(T1 a, int n){ return pow(a,n);}

template <typename T1>
class vecReg {
protected:
    T1 * _vals;
    hypercube<T1> _hyper;

public:
	vecReg() {}
	~vecReg() {delete [] _vals; _vals=nullptr;}
    vecReg(const hypercube<T1> &hyper){
        _hyper = hyper;
        _vals = new T1[hyper.getN123()];
    }
    vecReg(const std::shared_ptr<vecReg<T1> > vec ){
		_hyper = *(vec->getHyper());
		_vals = new T1[_hyper.getN123()];
        long num = _hyper.getN123()*sizeof(T1);
        std::memcpy(_vals,vec->getCVals(),num);
	}
    vecReg(const vecReg<T1> &  vec){
		_hyper = *(vec.getHyper());
		_vals = new T1[_hyper.getN123()];
        long num = _hyper.getN123()*sizeof(T1);
        std::memcpy(_vals,vec.getCVals(),num);
	}
    std::shared_ptr<vecReg<T1> > clone() const{
        std::shared_ptr<vecReg<T1> > newvec = std::make_shared<vecReg<T1> > (*this);
        return newvec;
    }
	const hypercube<T1> * getHyper() const {return &_hyper;}
    const T1 * getCVals() const {return _vals;}
    T1 * getVals() {return _vals;}
    long getN123() const {return _hyper.getN123();}
    void setHyper(const hypercube<T1> &hyper){
        if (_hyper.getN123()!= hyper.getN123()){
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::logic_error("Vector hypercube re-assignment is allowed only for same size N123.\n");
        }
        _hyper = hyper;
    }
    void zero(long imin=0, long imax=-1){
        if (imax==-1) imax=_hyper.getN123();
        #pragma omp parallel for
        for (long i=imin; i<imax; i++) _vals[i]=0;
    }
    void set(const T1 val, long imin=0, long imax=-1){
        if (imax==-1) imax=_hyper.getN123();
        #pragma omp parallel for
        for (long i=imin; i<imax; i++) _vals[i]=val;
    }
    T1 sum(long imin=0, long imax=-1) const {
        if (imax==-1) imax=_hyper.getN123();
        T1 val = 0;
        #pragma omp parallel for reduction(+: val)
        for (long i=imin; i<imax; i++) val += _vals[i];
        return val;
    }
    T1 norm(int p=2, long imin=0, long imax=-1) const {
        if (imax==-1) imax=_hyper.getN123();
        T1 val = 0;
        if (p == 2){
            #pragma omp parallel for reduction(+: val)
            for (long i=imin; i<imax; i++) val += _vals[i]*_vals[i];
            val = sqrt(val);
        }
        else {
            #pragma omp parallel for reduction(+: val)
            for (long i=imin; i<imax; i++) val += pow(std::abs(_vals[i]),p);
            val = pow(val,1.0/p);
        }
        return val;
    }
    T1 norm2(long imin=0, long imax=-1) const {
        if (imax==-1) imax=_hyper.getN123();
        T1 val = 0;
        #pragma omp parallel for reduction(+: val)
        for (long i=imin; i<imax; i++) val += _vals[i]*_vals[i];
        return val;
    }
    T1 rms(long imin=0, long imax=-1) const {
        if (imax==-1) imax=_hyper.getN123();
        T1 val = 0;
        #pragma omp parallel for reduction(+: val)
        for (long i=imin; i<imax; i++) val += _vals[i]*_vals[i];
        return sqrt(val/imax);
    }
    T1 min(long imin=0, long imax=-1) const {
        if (imax==-1) imax=_hyper.getN123();
        T1 val=_vals[imin];
        #pragma omp parallel for reduction(min: val)
        for (long i=imin+1; i<imax; i++) val=std::min(val,_vals[i]);
        return val;
    }
    T1 max(long imin=0, long imax=-1) const {
        if (imax==-1) imax=_hyper.getN123();
        T1 val=_vals[imin];
        #pragma omp parallel for reduction(max: val)
        for (long i=imin+1; i<imax; i++) val=std::max(val,_vals[i]);
        return val;
    }
    T1 absMax(long imin=0, long imax=-1) const {
        if (imax==-1) imax=_hyper.getN123();
        T1 val=0;
        #pragma omp parallel for reduction(max: val)
        for (long i=imin; i<imax; i++) val=std::max(val,std::abs(_vals[i]));
        return val;
    }
    T1 dot(const std::shared_ptr<vecReg<T1> > vec) const {
        if (this->getN123() != vec->getN123()){
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::logic_error("Vectors have different sizes.\n");
        }
        T1 val=0;
        long n123=_hyper.getN123();
        const T1 * pvec = vec->getCVals();
        #pragma omp parallel for reduction(+: val)
        for (long i=0; i<n123; i++) val += _vals[i]*pvec[i];
        return val;
    }
    void scale(const T1 sc, long imin=0, long imax=-1){
        if (imax==-1) imax=_hyper.getN123();
        #pragma omp parallel for
        for (long i=imin; i<imax; i++) _vals[i] *= sc;
    }
    void add(const T1 sc, long imin=0, long imax=-1){
        if (imax==-1) imax=_hyper.getN123();
        #pragma omp parallel for
        for (long i=imin; i<imax; i++) _vals[i] += sc;
    }
    void add(const std::shared_ptr<vecReg<T1> > vec){
        if (this->getN123() != vec->getN123()){
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::logic_error("Vectors have different sizes.\n");
        }
        const T1 * pvec = vec->getCVals();
        #pragma omp parallel for
        for (long i=0; i<_hyper.getN123(); i++) _vals[i] += pvec[i];
    }
    void scaleAdd(const std::shared_ptr<vecReg<T1> > vec, const T1 sc1=1.0, const T1 sc2=1.0){
        if (this->getN123() != vec->getN123()){
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::logic_error("Vectors have different sizes.\n");
        }
        const T1 * pvec = vec->getCVals();
        #pragma omp parallel for
        for (long i=0; i<_hyper.getN123(); i++) _vals[i] = sc1*_vals[i] + sc2*pvec[i];
    }
    void mult(const std::shared_ptr<vecReg<T1> > vec){
        if (this->getN123() != vec->getN123()){
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::logic_error("Vectors have different sizes.\n");
        }
        const T1 * pvec = vec->getCVals();
        #pragma omp parallel for
        for (long i=0; i<_hyper.getN123(); i++) _vals[i] *= pvec[i];
    }
    void divide(const std::shared_ptr<vecReg<T1> > vec, T1 eps=static_pow<T1> (10.0,-((int)sizeof(T1))*2)){
        const T1 * pvec = vec->getCVals();
        #pragma omp parallel for
        for (long i=0; i<_hyper.getN123(); i++) _vals[i] /= (pvec[i] + eps);
    }
    void random(T1 min=0, T1 max=1, unsigned int seed = 0){
        if (seed==0) {time_t timer; seed=time(&timer);}
        srand(seed);
        for (long i=0; i<_hyper.getN123(); i++) _vals[i] = min + (max-min)*(rand() % 10000)/10000.0;
    }
    void clip(const T1 min, const T1 max, long imin=0, long imax=-1){
        if (imax==-1) imax=_hyper.getN123();
        #pragma omp parallel for
        for (long i=imin; i<imax; i++) _vals[i] = std::max(min,std::min(max,_vals[i]));
    }
    void revert(const int dim = 1){
        if ((dim<1) || (dim > _hyper.getNdim())){
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::logic_error("The dimension to reverse does not exist.\n");
        }
        T1* vals = new T1[_hyper.getN123()];
        if (dim==1){
            int n1=_hyper.getAxis(1).n;
            int n = _hyper.getN123()/n1;
            #pragma omp parallel for
            for (int i=0; i<n; i++){
                for (int j=0; j<n1; j++){
                    vals[(long)i*n1+j] = _vals[(long)i*n1+n1-1-j];
                }
            }
            delete [] _vals;
            _vals = vals;
        }
        else if (dim==2){
            int n1=_hyper.getAxis(1).n;
            int n2=_hyper.getAxis(2).n;
            int n = _hyper.getN123()/(n1*n2);
            #pragma omp parallel for
            for (int i=0; i<n; i++){
                for (int j=0; j<n2; j++){
                    for (int k=0; k<n1; k++){
                        vals[(long)i*(n1*n2)+j*n1+k] = _vals[(long)i*(n1*n2)+(n2-1-j)*n1+k];
                    }
                }
            }
            delete [] _vals;
            _vals = vals;
        }

    }
    void reciprocal(){
        #pragma omp parallel for
        for (long i=0; i<_hyper.getN123(); i++) _vals[i] = 1.0/_vals[i];
    }
};

template <typename T1>
class cvecReg {
protected:
    std::complex<T1> * _vals;
    hypercube<T1> _hyper;
public:
	cvecReg() {}
	~cvecReg() {delete [] _vals; _vals=nullptr;}
    cvecReg(const hypercube<T1> &hyper){
        _hyper = hyper;
        _vals = new std::complex<T1>[hyper.getN123()];
    }
    cvecReg(const std::shared_ptr<cvecReg<T1> > vec ){
		_hyper = *(vec->getHyper());
		_vals = new std::complex<T1>[_hyper.getN123()];
        long num = _hyper.getN123()*sizeof(std::complex<T1>);
        std::memcpy(_vals,vec->getCVals(),num);
	}
    cvecReg(const cvecReg<T1> &  vec){
		_hyper = *(vec.getHyper());
		_vals = new std::complex<T1>[_hyper.getN123()];
        long num = _hyper.getN123()*sizeof(std::complex<T1>);
        std::memcpy(_vals,vec.getCVals(),num);
	}
    std::shared_ptr<cvecReg<T1> > clone() const{
        std::shared_ptr<cvecReg<T1> > newvec = std::make_shared<cvecReg<T1> > (*this);
        return newvec;
    }
	const hypercube<T1> * getHyper() const {return &_hyper;}
    const std::complex<T1> * getCVals() const {return _vals;}
    std::complex<T1> * getVals() {return _vals;}
    long getN123() const {return _hyper.getN123();}
    void setHyper(const hypercube<T1> &hyper){
        if (_hyper.getN123()!= hyper.getN123()){
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::logic_error("Vector hypercube re-assignment is allowed only for same size N123.\n");
        }
        _hyper = hyper;
    }
    void zero(long imin=0, long imax=-1){
        if (imax==-1) imax=_hyper.getN123();
        #pragma omp parallel for
        for (long i=imin; i<imax; i++) _vals[i]=0;
    }
    void set(const std::complex<T1> val, long imin=0, long imax=-1){
        if (imax==-1) imax=_hyper.getN123();
        #pragma omp parallel for
        for (long i=imin; i<imax; i++) _vals[i]=val;
    }
    std::complex<T1> sum(long imin=0, long imax=-1) const {
        if (imax==-1) imax=_hyper.getN123();
        data_t real=0, imag=0;
        #pragma omp parallel for reduction(+: real,imag)
        for (long i=imin; i<imax; i++) {real += _vals[i].rea(); imag += _vals[i].imag();}
        return std::complex<T1>(real,imag);
    }
    T1 norm(int p=2, long imin=0, long imax=-1) const {
        if (imax==-1) imax=_hyper.getN123();
        T1 val = 0;
        if (p == 2){
            #pragma omp parallel for reduction(+: val)
            for (long i=imin; i<imax; i++) val += std::norm(_vals[i]);
            val = sqrt(val);
        }
        else {
            #pragma omp parallel for reduction(+: val)
            for (long i=imin; i<imax; i++) val += pow(std::abs(_vals[i].real()),p)+pow(std::abs(_vals[i].imag()),p);
            val = pow(val,1.0/p);
        }
        return val;
    }
    T1 norm2(long imin=0, long imax=-1) const {
        if (imax==-1) imax=_hyper.getN123();
        T1 val = 0;
        #pragma omp parallel for reduction(+: val)
        for (int i=imin; i<imax; i++) val += std::norm(_vals[i]);
        return val;
    }
    T1 absMax(long imin=0, long imax=-1) const {
        if (imax==-1) imax=_hyper.getN123();
        T1 val=0;
        #pragma omp parallel for reduction(max: val)
        for (long i=imin; i<imax; i++) val=std::max(val,std::abs(_vals[i]));
        return val;
    }
    std::complex<T1> dot(const std::shared_ptr<cvecReg<T1> > vec) const {
        if (this->getN123() != vec->getN123()){
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::logic_error("Vectors have different sizes.\n");
        }
        data_t real=0, imag=0;
        long n123=_hyper.getN123();
        const std::complex<T1> * pvec = vec->getCVals();
        #pragma omp parallel for reduction(+: real,imag)
        for (long i=0; i<n123; i++) {
            real += _vals[i].real()*pvec[i].real()+_vals[i].imag()*pvec[i].imag();
            imag += _vals[i].real()*pvec[i].imag()-_vals[i].imag()*pvec[i].real();
        }
        return std::complex<T1>(real,imag);
    }
    void scale(const T1 sc, long imin=0, long imax=-1){
        if (imax==-1) imax=_hyper.getN123();
        #pragma omp parallel for
        for (long i=imin; i<imax; i++) _vals[i] *= sc;
    }
    void scale(const std::complex<T1> sc, long imin=0, long imax=-1){
        if (imax==-1) imax=_hyper.getN123();
        #pragma omp parallel for
        for (long i=imin; i<imax; i++) _vals[i] *= sc;
    }
    void add(const T1 sc, long imin=0, long imax=-1){
        if (imax==-1) imax=_hyper.getN123();
        #pragma omp parallel for
        for (long i=imin; i<imax; i++) _vals[i] += sc;
    }
    void add(const std::complex<T1> sc, long imin=0, long imax=-1){
        if (imax==-1) imax=_hyper.getN123();
        #pragma omp parallel for
        for (long i=imin; i<imax; i++) _vals[i] += sc;
    }
    void add(const std::shared_ptr<cvecReg<T1> > vec){
        if (this->getN123() != vec->getN123()){
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::logic_error("Vectors have different sizes.\n");
        }
        const std::complex<T1> * pvec = vec->getCVals();
        #pragma omp parallel for
        for (long i=0; i<_hyper.getN123(); i++) _vals[i] += pvec[i];
    }
    void scaleAdd(const std::shared_ptr<cvecReg<T1> > vec, const T1 sc1=1.0, const T1 sc2=1.0){
        if (this->getN123() != vec->getN123()){
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::logic_error("Vectors have different sizes.\n");
        }
        const std::complex<T1> * pvec = vec->getCVals();
        #pragma omp parallel for
        for (long i=0; i<_hyper.getN123(); i++) _vals[i] = sc1*_vals[i] + sc2*pvec[i];
    }
    void scaleAdd(const std::shared_ptr<cvecReg<T1> > vec, const std::complex<T1> sc1=1.0, const std::complex<T1> sc2=1.0){
        if (this->getN123() != vec->getN123()){
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::logic_error("Vectors have different sizes.\n");
        }
        const std::complex<T1> * pvec = vec->getCVals();
        #pragma omp parallel for
        for (long i=0; i<_hyper.getN123(); i++) _vals[i] = sc1*_vals[i] + sc2*pvec[i];
    }
    void mult(const std::shared_ptr<cvecReg<T1> > vec){
        if (this->getN123() != vec->getN123()){
            fprintf(stderr,"\nERROR in %s line %d\n",__FILE__,__LINE__);
            throw std::logic_error("Vectors have different sizes.\n");
        }
        const std::complex<T1> * pvec = vec->getCVals();
        #pragma omp parallel for
        for (long i=0; i<_hyper.getN123(); i++) _vals[i] *= pvec[i];
    }
    void divide(const std::shared_ptr<cvecReg<T1> > vec, T1 eps=static_pow<T1> (10.0,-((int)sizeof(T1))*2)){
        const std::complex<T1> * pvec = vec->getCVals();
        #pragma omp parallel for
        for (long i=0; i<_hyper.getN123(); i++) _vals[i] /= (pvec[i] + eps);
    }
    void random(T1 min=0, T1 max=1, unsigned int seed = 0){
        if (seed==0) {time_t timer; seed=time(&timer);}
        srand(seed);
        for (long i=0; i<_hyper.getN123(); i++) {
            _vals[i].real(min + (max-min)*(rand() % 10000)/10000.0);
            _vals[i].imag(min + (max-min)*(rand() % 10000)/10000.0);
        }
    }
    void clip(const T1 min, const T1 max, long imin=0, long imax=-1){
        if (imax==-1) imax=_hyper.getN123();
        #pragma omp parallel for
        for (long i=imin; i<imax; i++) {_vals[i].real(std::max(min,std::min(max,_vals[i].real()))); _vals[i].imag(std::max(min,std::min(max,_vals[i].imag())));}
    }
    void reciprocal(){
        #pragma omp parallel for
        for (long i=0; i<_hyper.getN123(); i++) _vals[i] = 1.0/_vals[i];
    }
    std::shared_ptr<cvecReg<T1> > conj() const{
        std::shared_ptr<cvecReg<T1> > newvec;
        newvec = std::make_shared<cvecReg<T1> > (*this);
        long n123=_hyper.getN123();
        std::complex<T1> * pvec = newvec->getVals();
        #pragma omp parallel for
        for (long i=0; i<n123; i++) pvec[i] = std::conj(_vals[i]);
        return newvec;
    }
    std::shared_ptr<vecReg<T1> > real() const{
        std::shared_ptr<vecReg<T1> > newvec;
        newvec = std::make_shared<vecReg<T1> > (_hyper);
        long n123=_hyper.getN123();
        T1 * pvec = newvec->getVals();
        #pragma omp parallel for
        for (long i=0; i<n123; i++) pvec[i] = _vals[i].real();
        return newvec;
    }
    std::shared_ptr<vecReg<T1> > imag() const{
        std::shared_ptr<vecReg<T1> > newvec;
        newvec = std::make_shared<vecReg<T1> > (_hyper);
        long n123=_hyper.getN123();
        T1 * pvec = newvec->getVals();
        #pragma omp parallel for
        for (long i=0; i<n123; i++) pvec[i] = _vals[i].imag();
        return newvec;
    }
    std::shared_ptr<vecReg<T1> > modulus() const{
        std::shared_ptr<vecReg<T1> > newvec;
        newvec = std::make_shared<vecReg<T1> > (_hyper);
        long n123=_hyper.getN123();
        T1 * pvec = newvec->getVals();
        #pragma omp parallel for
        for (long i=0; i<n123; i++) pvec[i] = std::abs(_vals[i]);
        return newvec;
    }
    std::shared_ptr<vecReg<T1> > modulus2() const{
        std::shared_ptr<vecReg<T1> > newvec;
        newvec = std::make_shared<vecReg<T1> > (_hyper);
        long n123=_hyper.getN123();
        T1 * pvec = newvec->getVals();
        #pragma omp parallel for
        for (long i=0; i<n123; i++) pvec[i] = std::norm(_vals[i]);
        return newvec;
    }
    std::shared_ptr<vecReg<T1> > arg() const{
        std::shared_ptr<vecReg<T1> > newvec;
        newvec = std::make_shared<vecReg<T1> > (_hyper);
        long n123=_hyper.getN123();
        T1 * pvec = newvec->getVals();
        #pragma omp parallel for
        for (long i=0; i<n123; i++) pvec[i] = std::arg(_vals[i]);
        return newvec;
    }
};