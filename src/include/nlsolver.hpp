#pragma once

#include "optimization.hpp"
#include "operator.hpp"

// General class for line search procedures to be used by the non-linear solvers
class lsearch{
protected:
    bool _flag_g; // Flag to indicate if the gradient needs to be computed after exiting the line search

public:
    lsearch(){}
    virtual ~lsearch(){}
    virtual lsearch * clone() const = 0;
    virtual void reset() = 0 ;
    bool getFlagG() {return _flag_g;}
    virtual bool lineSearch(optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose) = 0;
};

// Static line search with constant or scaled step size at every iteration
class static_ls : public lsearch{
protected:
    data_t _a0; // multiplier of the step length (yields step=a0/|g|)
    data_t _a1; // desired fraction of model update (yields step=a1*|m|/|g|), may overwrite a0 if a1>0 and |m|>0
    data_t _a2; // multiplier of the initial step length to be used in subsequent iterations (yields step=a2*step0), may overwite a1 if a2>0
    data_t _stp0; // first step length

public:
    static_ls(data_t a0=1, data_t a1=0, data_t a2=0){
        _a0 = a0;
        _a1 = a1;
        _a2 = a2;
        _stp0=0;
        _flag_g=true;
    }
    ~static_ls(){}

    static_ls * clone() const{
        static_ls * ls = new static_ls(_a0,_a1,_a2);
        return ls;
    }

    void reset(){
        _stp0=0;
    }

    bool lineSearch(optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose);
};

// Quadratic line search
class quadratic_ls : public lsearch{
protected:
    data_t _a0; // multiplier of the initial step length @ iter=0 (yields step0=a0/|g|)
    data_t _a1; // desired fraction of model update for the initial step length @ iter=0 (yields step0=a1*|m|/|g|), may overwrite a0 if a1>0 and |m|>0
    data_t _stp0; // most recent trial step
    data_t _stp; // most recent step length
    data_t _f, _df; // function and its derivative (with respect to step length) from current iteration
    data_t _f0; // most recent trial function evaluations
    bool _version; // defines the method for initializing the trial step length @ iter>0

public:
    quadratic_ls(data_t a0=1, data_t a1=0, bool version = 0){
        _a0 = a0;
        _a1 = a1;
        _version = version;
        _stp0 = 0;
        _stp = 0;
        _f=0;
        _df=0;
        _f0=0;
        _flag_g=true;
    }
    ~quadratic_ls(){}

    quadratic_ls * clone() const{
        quadratic_ls * ls = new quadratic_ls(_a0,_a1, _version);
        return ls;
    }

    void reset(){
        _stp0 = 0;
        _stp = 0;
        _f=0;
        _df=0;
        _f0=0;
    }

    bool lineSearch(optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose);
};

// Line search using quadratic/cubic interpolation and satisfying the sufficient descrease condition (Armijo condition)
class weak_wolfe : public lsearch{
protected:
    data_t _c1; // constant in (0,1) determining the sufficient descrease condition
    data_t _a0; // multiplier of the initial step length @ iter=0 (yields step0=a0/|g|)
    data_t _a1; // desired fraction of model update for the initial step length @ iter=0 (yields step0=a1*|m|/|g|), may overwrite a0 if a1>0 and |m|>0
    data_t _stp0, _stp1; // two most recent trial steps
    data_t _stp; // last accepted step length
    data_t _f, _df; // function and its derivative (with respect to step length) from current iteration
    data_t _f0, _f1; // two most recent trial function evaluations
    bool _version; // defines the method for initializing the trial step length @ iter>0

public:
    weak_wolfe(data_t c1=1e-4, data_t a0=1, data_t a1=0, bool version = 0){
        successCheck(((c1>0) && (c1<1)),__FILE__,__LINE__,"The constant c1 must be in (0,1).\n");
        _c1 = c1;
        _a0 = a0;
        _a1 = a1;
        _version = version;
        _stp0=0;
        _stp1=0;
        _stp=0;
        _f=0;
        _df=0;
        _f0=0;
        _f1=0;
        _flag_g=true;
    }
    ~weak_wolfe(){}

    weak_wolfe * clone() const{
        weak_wolfe * ls = new weak_wolfe(_c1,_a0,_a1,_version);
        return ls;
    }

    void reset(){
        _stp0=0;
        _stp1=0;
        _stp=0;
        _f=0;
        _df=0;
        _f0=0;
        _f1=0;
    }

    bool lineSearch(optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose);
};

// Line search using quadratic/cubic interpolation and satisfying the regular Wolfe conditions
class regular_wolfe : public lsearch{
protected:
    data_t _c1; // constant in (0,1) determining the sufficient descrease condition
    data_t _c2; // constant in (0,1) determining the curvature condition 0<c1<c2<1
    data_t _a0; // multiplier of the initial step length @ iter=0 (yields step0=a0/|g|)
    data_t _a1; // desired fraction of model update for the initial step length @ iter=0 (yields step0=a1*|m|/|g|), may overwrite a0 if a1>0 and |m|>0
    data_t _stp0, _stp1; // step length interval containing the desired step
    data_t _stp; // last accepted step length
    data_t _stp_max; // maximum step length allowed at any iteration
    data_t _f, _df; // function and its derivative (with respect to step length) from current iteration
    data_t _f0, _f1; // two most recent trial function evaluations
    bool _version; // defines the method for initializing the trial step length @ iter>0

public:
    regular_wolfe(data_t c1=1e-4, data_t c2=0.9, data_t a0=1, data_t a1=0, data_t step_max=1e6, bool version = 0){
        successCheck(((c1>0) && (c1<1) && (c2>c1) && (c2<1)),__FILE__,__LINE__,"The constants c1, c2 must satisfy 0<c1<c2<1\n");
        _c1 = c1;
        _c2 = c2;
        _a0 = a0;
        _a1 = a1;
        _stp_max = step_max;
        _version = version;
        _stp0=0;
        _stp1=0;
        _stp=0;
        _f=0;
        _df=0;
        _f0=0;
        _f1=0;
        _flag_g=false;
    }
    virtual ~regular_wolfe(){}

    regular_wolfe * clone() const{
        regular_wolfe * ls = new regular_wolfe(_c1,_c2,_a0,_a1,_stp_max,_version);
        return ls;
    }

    void reset(){
        _stp0=0;
        _stp1=0;
        _stp=0;
        _f=0;
        _df=0;
        _f0=0;
        _f1=0;
    }

    bool lineSearch(optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose);
    
    bool zoom(bool reverse,
                        optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int &trial,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose);
};

// Line search using quadratic/cubic interpolation and satisfying the strong Wolfe conditions
class strong_wolfe : public regular_wolfe{
public:
    strong_wolfe(data_t c1=1e-4, data_t c2=0.9, data_t a0=1, data_t a1=0, data_t step_max=1e6, bool version = 0){
        successCheck(((c1>0) && (c1<1) && (c2>c1) && (c2<1)),__FILE__,__LINE__,"The constants c1, c2 must satisfy 0<c1<c2<1\n");
        _c1 = c1;
        _c2 = c2;
        _a0 = a0;
        _a1 = a1;
        _stp_max = step_max;
        _version = version;
        _stp0=0;
        _stp1=0;
        _stp=0;
        _f=0;
        _df=0;
        _f0=0;
        _f1=0;
        _flag_g=false;
    }
    ~strong_wolfe(){};
    strong_wolfe * clone() const{
        strong_wolfe * ls = new strong_wolfe(_c1,_c2,_a0,_a1,_stp_max,_version);
        return ls;
    }
    bool lineSearch(optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose);
    
    bool zoom(bool reverse,
                        optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int &trial,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose);
};

// General class for iterative solver for non-linear optimization problems
class nlsolver{
protected:
    int _niter; // number of iterations
    int _feval; // number of function evaluations
    int _geval; // number of gradient evaluations
    data_t _threshold; // stopper based on rate of convergence
    int _max_trial; // maximum number of trial step lengths at any given iteration
    lsearch * _lsearch; // line search procedure

public:
    std::vector <data_t> _func; // vector of objective function values
    nlsolver(){}
    virtual ~nlsolver(){delete _lsearch;}
    nlsolver(const int niter, const int max_trial, const data_t threshold=0, lsearch * ls=nullptr){
        _niter = niter;
        _max_trial = max_trial;
        _threshold = threshold;
        _feval = 0;
        _geval = 0;
        _func = {};

        if (ls==nullptr) _lsearch = new weak_wolfe(1e-4,1,0,0);
        else _lsearch = ls->clone();
    }
    int getNiter() const {return _niter;}
    data_t getFeval() const {return _feval;}
    data_t getGeval() const {return _geval;}
    void setNiter(int niter) {_niter = niter;}
    void setMaxTrial(int max_trial) {_max_trial = max_trial;}
    void setThreshold(data_t threshold){_threshold = threshold;}

    // reset variables
    void reset(){
        _func.clear();
        _feval = 0;
        _geval = 0;
        _lsearch->reset();
    }

    // run the solver for a given problem
    virtual void run(optimization * prob, const bool verbose=false, std::string output="none", int isave = 10, int format = 0, std::string datapath="none") = 0;

    // test the solver for a full rank system Lm = d with random vectors
    void testLinear(const bool verbose=true);

    // test the solver for Rosenbrock function f(x,y)=(a-x)^2 + b(y-x^2)^2 true minimum = 0 at (a,a^2)
    void testRosenbrock(const bool verbose=true);

    // test the solver for the paraboloid function f(x,y)=1/2(ax)^2 + 1/2(by)^2 true minimum = 0 at (0,0)
    void testParaboloid(const bool verbose=true);
};

// Non-linear steepest descent solver
class nlsd : public nlsolver{
public:

    using nlsolver::nlsolver;
    ~nlsd(){}
    void run(optimization * prob, const bool verbose=false, std::string output="none", int isave = 10, int format = 0, std::string datapath="none");
};

// Non-linear conjugate gradient solver
class nlcg : public nlsolver{
    int _method; // 0 for Fletcher-Reeves FR, 1 for Polak-Ribiere PR, 2 for Polak-Ribiere PR+
                 //PR and PR+ should be avoided when Wolfe or strong Wolfe line search is used
public:

    nlcg(){}
    ~nlcg(){}
    nlcg(const int niter, const int max_trial, const data_t threshold=0, lsearch * ls=nullptr):nlsolver(niter,max_trial,threshold,ls){
        _method = 0;
    }
    void setMethod(int method){
        if (method == 1) _method = 1;
        else if (method == 2) _method = 2;
        else _method = 0;
    }
    void run(optimization * prob, const bool verbose=false, std::string output="none", int isave = 10, int format = 0, std::string datapath="none");
};

// BFGS solver
class bfgs : public nlsolver{
protected:
    matrix * _H; // Approximate inverse Hessian; must be SPD
public:

    bfgs(){}
    ~bfgs(){delete _H;}
    bfgs(const int niter, const int max_trial, const data_t threshold=0, lsearch * ls=nullptr, std::shared_ptr<vecReg<data_t> > H0=nullptr):nlsolver(niter,max_trial,threshold,ls){
        if (H0 != nullptr) {
            successCheck(((H0->getHyper()->getNdim()== 2) && (H0->getHyper()->getAxis(1).n == H0->getHyper()->getAxis(2).n)),__FILE__,__LINE__,"Initial inverse Hessian must be a square matrix.\n");
            _H = new matrix(H0);
            successCheck(testSPD(_H),__FILE__,__LINE__,"The provided initial Hessian matrix must be SPD\n");
        }
        else _H=nullptr;
    }

    void run(optimization * prob, const bool verbose=false, std::string output="none", int isave = 10, int format = 0, std::string datapath="none");
    bool testSPD(matrix * H);
    void updateH(matrix * H ,std::shared_ptr<vecReg<data_t> > s, std::shared_ptr<vecReg<data_t> > y);
};

// l-BFGS solver
class lbfgs : public nlsolver{
protected:
    std::shared_ptr<vecReg<data_t> > _H0; // Initial guess for the diagonal inverse Hessian; must be SPD
    int _m; // number of previous s,y pairs to keep for the l-BFGS update
public:

    lbfgs(){}
    ~lbfgs(){}
    lbfgs(const int niter, const int max_trial, const data_t threshold=0, lsearch * ls=nullptr, std::shared_ptr<vecReg<data_t> > H0=nullptr, int m=5):nlsolver(niter,max_trial,threshold,ls){
        if (H0 != nullptr) {
            successCheck((H0->min() > 0),__FILE__,__LINE__,"Initial diagonal inverse Hessian must be definite positive.\n");
            _H0 = H0->clone();
        }
        else _H0=nullptr; 
        _m = m;
    }

    void run(optimization * prob, const bool verbose=false, std::string output="none", int isave = 10, int format = 0, std::string datapath="none");
    void updateSY(std::vector<std::shared_ptr<vecReg<data_t> > > &s,
                    std::vector<std::shared_ptr<vecReg<data_t> > > &y,
                    std::vector<data_t> &rho,
                    data_t * pm, data_t * pm0,
                    data_t * pg, data_t * pg0,
                    int iter);
    void computeHg(std::shared_ptr<vecReg<data_t> > &p, std::shared_ptr<vecReg<data_t> > g,
                    std::vector<std::shared_ptr<vecReg<data_t> > > s, std::vector<std::shared_ptr<vecReg<data_t> > > y, std::vector<data_t> rho, int iter);
};

// Newton solver
class newton : public nlsolver{

public:

    newton(){}
    ~newton(){}
    newton(const int niter, const int max_trial, const data_t threshold=0, lsearch * ls=nullptr):nlsolver(niter,max_trial,threshold,ls){}
    void run(optimization * prob, const bool verbose=false, std::string output="none", int isave = 10, int format = 0, std::string datapath="none");
};