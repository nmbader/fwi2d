#pragma once

#include "optimization.hpp"

// General class for iterative solver for linear optimization problems
class lsolver{
protected:
    int _niter; // number of iterations
    data_t _threshold; // stopper based on rate of convergence
    bool _successful; // flag to indicate if the solver finished successfully

public:
    std::vector <data_t> _func; // vector of objective function values
    lsolver(){}
    virtual ~lsolver(){}
    lsolver(const int niter, const data_t threshold=0){
        _niter = niter;
        _threshold = threshold;
        _successful = true;
        _func = {};
    }
    int getNiter() const {return _niter;}
    data_t getThreshold() const {return _threshold;}
    bool isSuccessful() const {return _successful;}
    void setNiter(int niter) {_niter = niter;}
    void setThreshold(data_t threshold){_threshold = threshold;}

    // run the solver for a given problem
    virtual void run(llsq * prob, const bool verbose=false) = 0;

    // reset the objective function
    void resetFunc(){
        _func.clear();
        _successful = true;
    }

    // test the solver for a full rank system Lm = d with random vectors
    void test(const bool verbose=true);

    // test the iterative solver for min(1/2.|Lm-d|^2+1/2.lambda.|D(m-m_prior)|^2)
    void testReg(const bool verbose=true);

    // test the solver for a SPD system Lm = d with random vectors
    void testSPD(const bool verbose=true);
};

// Linear steepest descent least squares SDLS solver
class sdls : public lsolver{
public:

    using lsolver::lsolver;
    ~sdls(){}
    void run(llsq * prob, const bool verbose=false);
};

// Linear conjugate gradient CG solver (for SPD system)
// according to Aster (Parameter Estimation and Inverse Problems), chap 6
class cg : public lsolver{
public:
    using lsolver::lsolver;
    ~cg(){}
    void run(llsq * prob, const bool verbose=false);
};

// Linear conjugate gradient least squares CGLS solver
// according to Aster (Parameter Estimation and Inverse Problems), chap 6
class cgls : public lsolver{
public:

    using lsolver::lsolver;
    ~cgls(){}
    void run(llsq * prob, const bool verbose=false);
};