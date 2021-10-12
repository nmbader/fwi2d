#include <iostream>
#include "lsolver.hpp"
#include "operator.hpp"

void lsolver::test(const bool verbose){

    data_t seed = 1234;
    std::shared_ptr<vecReg<data_t> > mat = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100,100));
    mat->random(-1,1,seed);
    matrix M(mat);
    std::shared_ptr<vecReg<data_t> > mod = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100));
    std::shared_ptr<vecReg<data_t> > mod0 = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100));
    std::shared_ptr<vecReg<data_t> > dat = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100));

    mod->random(-1,1,seed);
    mod0->zero();
    M.forward(false,mod,dat);
    llsq prob(&M, mod0,dat);

    this->run(&prob, verbose);

    // L2 squared error of estimated model
    data_t trueNorm, estimatedNorm, dataNorm, error, residual;
    trueNorm = mod->norm2();
    estimatedNorm = mod0->norm2();
    dataNorm = dat->norm2();
    mod->scaleAdd(mod0,1,-1);
    error = mod->norm2();
    dat->scale(-1);
    M.forward(true,mod0,dat);
    residual = dat->norm2();

    std::clog << "data norm squared = " << dataNorm << std::endl;
    std::clog << "true model norm squared = " << trueNorm << std::endl;
    std::clog << "estimated model norm squared = " << estimatedNorm << std::endl;
    std::clog << "residual norm squared = " << residual << std::endl;
    std::clog << "model error squared = " << error << std::endl;

}

void lsolver::testReg(const bool verbose){

    std::shared_ptr<vecReg<data_t> > mat = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100,100));
    mat->random(-1,1);
    matrix M(mat);
    std::shared_ptr<vecReg<data_t> > mod = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100));
    std::shared_ptr<vecReg<data_t> > mod0 = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100));
    std::shared_ptr<vecReg<data_t> > mod1 = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100));
    std::shared_ptr<vecReg<data_t> > mod_prior = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100));
    std::shared_ptr<vecReg<data_t> > dat = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100));

    mod->random(-1,1);
    mod_prior->random(-10,10);
    mod0->zero();
    mod1->zero();
    M.forward(false,mod,dat);

    matrix Id (mod->getN123(),1);
    data_t la0 = 1e-03;
    data_t la1 = 1e+03;

    llsq_reg prob0(&M, &Id, mod0,dat,la0,mod_prior);
    llsq_reg prob1(&M, &Id, mod1,dat,la1,mod_prior);

    this->run(&prob0,verbose);
    this->run(&prob1,verbose);

    // L2 squared error of estimated model
    data_t trueNorm, priorNorm, estimatedNorm0, estimatedNorm1, dataNorm, error0, error1;
    trueNorm = mod->norm2();
    priorNorm = mod_prior->norm2();
    estimatedNorm0 = mod0->norm2();
    estimatedNorm1 = mod1->norm2();
    dataNorm = dat->norm2();
    mod0->scaleAdd(mod,-1,1);
    error0 = mod0->norm2();
    mod1->scaleAdd(mod,-1,1);
    error1 = mod1->norm2();

    std::clog << "data norm squared = " << dataNorm << std::endl;
    std::clog << "true model norm squared = " << trueNorm << std::endl;
    std::clog << "prior model norm squared = " << priorNorm << std::endl;
    std::clog << "estimated model norm squared (weak reg) = " << estimatedNorm0 << std::endl;
    std::clog << "estimated model norm squared (strong reg) = " << estimatedNorm1 << std::endl;
    std::clog << "model error squared (weak reg) = " << error0 << std::endl;
    std::clog << "model error squared (strong reg) = " << error1 << std::endl;
}

void lsolver::testSPD(const bool verbose){

    std::shared_ptr<vecReg<data_t> > mat = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100,100));
    std::shared_ptr<vecReg<data_t> > diag = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100));
    std::shared_ptr<vecReg<data_t> > mod = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100));
    std::shared_ptr<vecReg<data_t> > mod0 = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100));
    std::shared_ptr<vecReg<data_t> > dat = std::make_shared<vecReg<data_t> > (hypercube<data_t>(100));

    // Build a SPD matrix M = LDL'
    diag->random(1,2);
    mat->random(-0.01,0.01);
    for (int i=0; i<100; i++){
        for (int j=i+1; j<100; j++) mat->getVals()[100*i+j] = mat->getVals()[100*j+i];
        mat->getVals()[100*i+i] = diag->getVals()[i];
    }
    
    matrix M(mat);
    mod->random(-1,1);
    mod0->set(1);
    M.forward(false,mod,dat);
    lspd prob(&M, mod0,dat);

    this->run(&prob, verbose);

    // L2 squared error of estimated model
    data_t trueNorm, estimatedNorm, dataNorm, error, residual;
    trueNorm = mod->norm2();
    estimatedNorm = mod0->norm2();
    dataNorm = dat->norm2();
    mod->scaleAdd(mod0,1,-1);
    error = mod->norm2();
    dat->scale(-1);
    M.forward(true,mod0,dat);
    residual = dat->norm2();

    std::clog << "data norm squared = " << dataNorm << std::endl;
    std::clog << "true model norm squared = " << trueNorm << std::endl;
    std::clog << "estimated model norm squared = " << estimatedNorm << std::endl;
    std::clog << "residual norm squared = " << residual << std::endl;
    std::clog << "model error squared = " << error << std::endl;
}

void sdls::run(llsq * prob, const bool verbose){

    prob->initGrad();
    prob->initRes();
    prob->initDRes();

    // bottom limit of the objective function
    data_t ZERO = prob->getZero();
    
    //residuals
    prob->res();
    std::shared_ptr<vecReg<data_t> > r = prob->getRes();

    //dresiduals
    std::shared_ptr<vecReg<data_t> > dr = prob->getDRes();

    //model
    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    
    //gradients
    prob->grad();
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();

    //objective function
    _func.clear();
    _func.push_back(prob->getFunc());

    if (_func[0] <= ZERO){
        std::clog << "WARNING: The initial objective function is close to zero. The solver will stop.\n";
        _successful = false;
        return;
    }
    if (std::isnan(_func[0]) || (_func[0]!=_func[0])) throw std::runtime_error("==============================\nERROR: the objective function is NaN.\n==============================\n");

    //step length
    data_t alpha = 0;

    //iteration number
    int k = 0;

    //convergence rate
    data_t rate = 1;

if (verbose) std::clog << "Iteration = 0; functional = "<<_func[0]<<"; normalized functional = 1\n";

    data_t gnorm2 = g->norm2();
    if (gnorm2 <= ZERO){
        std::clog << "WARNING: The initial gradient is close to zero. The solver will stop.\n";
        _successful = false;
        return;
    }

    //start the SD loop consisting of 4 steps
    while (k<_niter & rate>_threshold){

        //step 1: alpha = grad^2 / (L.grad)^2
        prob->dres(g);
        alpha = gnorm2 / dr->norm2();

        //step 2: update the model
        m->scaleAdd(g,1,-alpha);

        //step 3: update residuals
        r->scaleAdd(dr,1,-alpha);
        
        //store the normalized residuals
        data_t f = prob->getFunc();
if (verbose) std::clog << "Iteration = "<<k+1<<"; functional = "<<f<<"; ";
if (verbose) std::clog <<"normalized functional = "<<f / _func[0]<<"\n";
        _func.push_back(f);

        //compute the rate of convergence
        if (_func[k+1] <= ZERO){
            std::clog << "WARNING: The objective function is close to zero. The solver will stop.\n";
            _successful = false;
            k++;
            break;
        }
        if(k>0){
            rate = (_func[k]-_func[k+1])/(std::abs(_func[k])+std::abs(_func[k+1]));
        }

        //step4: compute the gradient
        prob->grad();
        gnorm2 = g->norm2();
        if (gnorm2 <= ZERO){
            std::clog << "WARNING: The gradient is close to zero. The solver will stop.\n";
            _successful = false;
            k++;
            break;
        }

        //iterate
        k++;
    }

if (verbose) std::clog << "Total number of SD iterations: "<<k<<"\n";
}

void cgls::run(llsq * prob, const bool verbose){

    prob->initGrad();
    prob->initRes();
    prob->initDRes();

    // bottom limit of the objective function
    data_t ZERO = prob->getZero();
    
    //residuals
    prob->res();
    std::shared_ptr<vecReg<data_t> > r = prob->getRes();

    //dresiduals
    std::shared_ptr<vecReg<data_t> > dr = prob->getDRes();

    //model
    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    
    //gradients
    prob->grad();
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();

    //conjugate direction
    std::shared_ptr<vecReg<data_t> > p = std::make_shared<vecReg<data_t> >(*g->getHyper());
    p->zero();

    //directions weights
    data_t alpha = 0;
    data_t beta = 0;

    //objective function
    _func.clear();
    _func.push_back(prob->getFunc());

    if (_func[0] <= ZERO){
        std::clog << "WARNING: The initial objective function is close to zero. The solver will stop.\n";
        _successful = false;
        return;
    }

    //iteration number
    int k = 0;

    //convergence rate
    data_t rate = 1;

if (verbose) std::clog << "Iteration = 0; functional = "<<_func[0]<<"; normalized functional = 1\n";

    data_t gnorm2 = g->norm2();
    if (gnorm2 <= ZERO){
        std::clog << "WARNING: The initial gradient is close to zero. The solver will stop.\n";
        _successful = false;
        return;
    }

    //start the CG loop consisting of 7 steps
    while (k<_niter & rate>_threshold){

        //step 1: p = -g + beta*p
        p->scaleAdd(g,beta,-1);

        //step 2: compute alpha = g^2/(Lp)^2
        prob->dres(p);
        alpha = gnorm2 / dr->norm2();

        //step 3: update model
        m->scaleAdd(p,1,alpha);

        //step 4: update residuals
        r->scaleAdd(dr,1,alpha);
        
        //store the residuals norm2
        data_t f = prob->getFunc();
if (verbose) std::clog << "Iteration = "<<k+1<<"; functional = "<<f<<"; ";
if (verbose) std::clog <<"normalized functional = "<<f / _func[0]<<"\n";
        _func.push_back(f);

        //compute the rate of convergence
        if (_func[k+1] <= ZERO){
            std::clog << "WARNING: The objective function is close to zero. The solver will stop.\n";
            _successful = false;
            k++;
            break;
        }
        if(k>0){
            rate = (_func[k]-_func[k+1])/(std::abs(_func[k])+std::abs(_func[k+1]));
        }

        //step 5: update gradient
        beta = 1.0/gnorm2;
        prob->grad();

        //step 6: compute beta
        gnorm2 = g->norm2();
        if (gnorm2 <= ZERO){
            std::clog << "WARNING: The gradient is close to zero. The solver will stop.\n";
            _successful = false;
            k++;
            break;
        }
        beta *= gnorm2;

        //step 7: iterate
        k++;
    }

if (verbose) std::clog << "Total number of CG iterations: "<<k<<"\n";
}

void cg::run(llsq * prob, const bool verbose){

    prob->initGrad();
    prob->initRes();
    prob->initDRes();

    // bottom limit of the objective function
    data_t ZERO = prob->getZero();

    //residuals
    prob->res();
    std::shared_ptr<vecReg<data_t> > r = prob->getRes();

    //dresiduals
    std::shared_ptr<vecReg<data_t> > dr = prob->getDRes();

    //model
    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    
    //gradients
    prob->grad();
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();

    //conjugate direction
    std::shared_ptr<vecReg<data_t> > p = std::make_shared<vecReg<data_t> >(*g->getHyper());
    p->zero();

    //directions weights
    data_t alpha = 0;
    data_t beta = 0;

    //objective function
    _func.clear();
    _func.push_back(prob->getFunc());

    if (_func[0] <= ZERO){
        std::clog << "WARNING: The initial objective function is too small. The solver will stop.\n";
        _successful = false;
        return;
    }

    //iteration number
    int k = 0;

    //convergence rate
    data_t rate = 1;

if (verbose) std::clog << "Iteration = 0; functional = "<<_func[0]<<"; normalized functional = 1\n";

    data_t gnorm2 = g->norm2();
    if (gnorm2 <= ZERO){
        std::clog << "WARNING: The initial gradient is close to zero. The solver will stop.\n";
        _successful = false;
        return;
    }

    //start the CG loop consisting of 6 steps
    while (k<_niter & rate>_threshold){

        //step 1: p = -g + beta*p
        p->scaleAdd(g, beta, -1);

        //step 2: compute alpha = g^2/p'.L.p
        prob->dres(p);
        alpha = gnorm2 / p->dot(dr);

        //step 3: update model
        m->scaleAdd(p,1,alpha);

        //step 4: update residuals (and gradient)
        beta = gnorm2;
        r->scaleAdd(dr,1,alpha);
        gnorm2 = r->norm2();
        
        //store the residuals norm2
        data_t f = prob->getFunc();
if (verbose) std::clog << "Iteration = "<<k+1<<"; functional = "<<f<<"; ";
if (verbose) std::clog <<"normalized functional = "<<f / _func[0]<<"\n";
        _func.push_back(f);

        //compute the rate of convergence
        if (_func[k+1] <= ZERO){
            std::clog << "WARNING: The objective function is too small. The solver will stop.\n";
            _successful = false;
            k++;
            break;
        }
        if(k>0){
            rate = (_func[k]-_func[k+1])/(std::abs(_func[k])+std::abs(_func[k+1]));
        }

        //step 5: compute beta
        if (gnorm2 <= ZERO){
            std::clog << "WARNING: The gradient is close to zero. The solver will stop.\n";
            _successful = false;
            k++;
            break;
        }
        beta = gnorm2 / beta;

        //step 6: iterate
        k++;
    }

if (verbose) std::clog << "Total number of CG iterations: "<<k<<"\n";
}