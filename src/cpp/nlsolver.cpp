#include <iostream>
#include "nlsolver.hpp"
#include "IO.hpp"
#include "misc.hpp"

#define ZERO 1e-10
#define STOL 1e-3

data_t quadratic(data_t f, data_t df, data_t stp0, data_t f0){
    data_t ans=0;
    if ( f0 <= f + df*stp0) ans=stp0*1.1;
    else ans = -df*stp0*stp0/(2*(f0-f-df*stp0));
    if ((ans/stp0 < STOL) || (std::abs(ans-stp0)/stp0) < STOL) ans = stp0/2.0;
    return ans;
}

data_t cubic(data_t f, data_t df, data_t stp0, data_t f0, data_t stp1, data_t f1){
    data_t x1 = f1-f-df*stp1;
    data_t x2 = f0-f-df*stp0;
    data_t x3 = stp0*stp0*stp1*stp1*(stp1-stp0);
    data_t a = (stp0*stp0*x1-stp1*stp1*x2)/x3;
    data_t b = (-stp0*stp0*stp0*x1+stp1*stp1*stp1*x2)/x3;
    data_t delta = b*b - 3*a*df;
    data_t step, func, ans;
    if (f1<f0){step=stp1; func=f1;}
    else {step=stp0; func=f0;}
    if (std::abs(a)<ZERO || delta<0) ans = quadratic(f,df,step,func);
    else ans = (-b+sqrt(delta))/(3*a);
    if ((ans/step < STOL) || (std::abs(step-ans)/step) < STOL) ans = step/2.0;
    return ans;
}

bool weak_wolfe::lineSearch(optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose){

    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();
    int n = m->getN123();
    data_t * pm = m->getVals();
    data_t * pm0 = m0->getVals();
    data_t * pp = p->getVals();
    data_t alpha=0,df0=0,f=0;

    int trial=0;
    while (true){
        if (trial > max_trial){
            fprintf(stderr,"==============================\nWARNING: Maximum number of trials reached before finding a satisfactory step length. Solver STOPS\n==============================\n");
            return false;
        }
        
        // first trial
        if (trial == 0){
            df0 = _df;
            _f = prob->getFunc();
            _df = g->dot(p);
            f=_f;

            // first iteration
            if (iter == 0) {
                _stp0 = _a0/gnorm;
                if (_a1>0) _stp0 = std::max(_stp0,_a1*m->norm()/gnorm);
                trial++;
                alpha=_stp0;
                continue;
            }

            // subsequent iterations
            else {
                _stp0 = _stp*df0/_df;
                if (_version==1){
                    _stp0 = 2*(_f - _f0)/_df;
                    _stp0 = std::min(1.0,1.00001*_stp0);
                }
                trial++;
                alpha=_stp0;
                continue;
            }
        }
     
        // second trial
        else if (trial==1){
            // update the model
            #pragma omp parallel for
            for (int i=0; i<n; i++) pm[i] = pm0[i] + alpha*pp[i];

            //update residuals
            prob->res();

            //compute objective function
            feval++;
            f = prob->getFunc();
            _f0=f;

            if (f!=f) throw std::runtime_error("==============================\nERROR: the objective function is NaN.\n==============================\n");
            if (verbose){
                fprintf(stderr,"==============================\n");
                fprintf(stderr,"Trial = %d\n",trial);
                fprintf(stderr,"Step length = %.10f\n",alpha);
                fprintf(stderr,"Function value = %.10f\n",f);
                fprintf(stderr,"Last gradient norm = %f\n",gnorm);
            }

            if (f <= _f + _c1*_stp0*_df || f<=ZERO){
                _stp=_stp0;
                _f0 = _f;
                _f = f;
                if (verbose){
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Found a satisfactory step length\n");
                    if (f<=ZERO) fprintf(stderr,"Function value became too small\n");
                    fprintf(stderr,"==============================\n");
                }
                if (f<=ZERO) return false;
                return true;
            }
            else{
                _stp1 = quadratic(_f,_df,_stp0,_f0);
                alpha = _stp1;
                trial++;
                continue;
            }
        }

        // subsequent trials
        else {
            if (trial>2) _f0 = _f1;
            // update the model
            #pragma omp parallel for
            for (int i=0; i<n; i++) pm[i] = pm0[i] + alpha*pp[i];

            //update residuals
            prob->res();

            //compute objective function
            feval++;
            f = prob->getFunc();
            _f1 = f;

            if (f!=f) throw std::runtime_error("==============================\nERROR: the objective function is NaN.\n==============================\n");
            if (verbose){
                fprintf(stderr,"==============================\n");
                fprintf(stderr,"Trial = %d\n",trial);
                fprintf(stderr,"Step length = %.10f\n",alpha);
                fprintf(stderr,"Function value = %.10f\n",f);
                fprintf(stderr,"Last gradient norm = %f\n",gnorm);
            }

            if (f <= _f + _c1*_stp1*_df || f<=ZERO){
                _f0 = _f;
                _f=f;
                _stp=_stp1;
                if (verbose){
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Found a satisfactory step length\n");
                    if (f<=ZERO) fprintf(stderr,"Function value became too small\n");
                    fprintf(stderr,"==============================\n");
                }
                if (f<=ZERO) return false;
                return true;
            }
            else{
                data_t temp = cubic(_f,_df,_stp0,_f0,_stp1,_f1);
                _stp0=_stp1;
                _stp1=temp;
                alpha=_stp1;
                trial++;
                continue;
            }
        }
    }
}

bool regular_wolfe::lineSearch(optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose){

    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();
    int n = m->getN123();
    data_t * pm = m->getVals();
    data_t * pm0 = m0->getVals();
    data_t * pp = p->getVals();
    data_t alpha=0,df0=0,df1=0,f=0;

    int trial=0;
    // first trial
    df0 = _df;
    _f = prob->getFunc();
    _df = g->dot(p);
    f=_f;
    _stp0=0;

    // first iteration
    if (iter == 0) {
        _f0=_f;
        _stp1 = _a0/gnorm;
        if (_a1>0) _stp1 = std::max(_stp1,_a1*m->norm()/gnorm);
        alpha=_stp1;
    }

    // subsequent iterations
    else {
        _stp1 = _stp*df0/_df;
        if (_version==1){
            _stp1 = 2*(_f - _f0)/_df;
            _stp1 = std::min(1.0,1.00001*_stp1);
        }
        alpha=_stp1;
    }

    while (true){
        if (trial >= max_trial){
            fprintf(stderr,"==============================\nWARNING: Maximum number of trials reached before finding a satisfactory step length. Solver STOPS\n==============================\n");
            return false;
        }

        // Evaluate phi(alpha)
        #pragma omp parallel for
        for (int i=0; i<n; i++) pm[i] = pm0[i] + alpha*pp[i];
        prob->res();
        feval++;
        trial++;
        f = prob->getFunc();
        _f1= f;
        if (f!=f) throw std::runtime_error("==============================\nERROR: the objective function is NaN.\n==============================\n");

        if (f<=ZERO){
            _stp=alpha;
            _f0=_f;
            _f = f;
            if (verbose){
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Trial = %d\n",trial);
                    fprintf(stderr,"Step length = %f\n",alpha);
                    fprintf(stderr,"Function value = %f\n",f);
                    fprintf(stderr,"Last gradient norm = %f\n",gnorm);
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Found a satisfactory step length\n");
                    fprintf(stderr,"Function value became too small\n");
                    fprintf(stderr,"==============================\n");
            }
            return false;
        }
        // check if (_stp0, _stp1) contains the desired step length
        if ((_f1>_f+_c1*alpha*_df) || (_f1>=_f0 && trial>1)){
            if (verbose){
                fprintf(stderr,"==============================\n");
                fprintf(stderr,"Trial = %d\n",trial);
                fprintf(stderr,"Step length = %f\n",alpha);
                fprintf(stderr,"Function value = %f\n",f);
                fprintf(stderr,"Last gradient norm = %f\n",gnorm);
            }
            return zoom(false,prob,m0,p,iter,trial,max_trial,feval,geval,gnorm,verbose);
        }

        // Evaluate phi'(stp1)
        prob->grad();
        geval++;
        gnorm=g->norm();
        df1 = g->dot(p);

        if (verbose){
                fprintf(stderr,"==============================\n");
                fprintf(stderr,"Trial = %d\n",trial);
                fprintf(stderr,"Step length = %f\n",alpha);
                fprintf(stderr,"Function value = %f\n",f);
                fprintf(stderr,"Gradient norm = %f\n",gnorm);
        }

        // check the curvature condition
        if (df1 >= _c2*_df){
            _stp=alpha;
            _f0=_f;
            _f = f;
            if (verbose){
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Found a satisfactory step length\n");
                    fprintf(stderr,"==============================\n");
            }
            return true;
        }

        // Check if phi(stp1) is increasing
        if (df1 >= 0){
            return zoom(true,prob,m0,p,iter,trial,max_trial,feval,geval,gnorm,verbose);
        }

        // Choose another step in (stp1,stp_max)
        _stp0 = _stp1;
        _stp1 = quadratic(_f,_df,_stp1,_f1);
        _stp1 = std::min(_stp1, _stp_max);
        alpha = _stp1;
        _f0=_f1;
    }   
}

bool regular_wolfe::zoom(bool reverse,
                        optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int &trial,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose){

    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();
    int n = m->getN123();
    data_t * pm = m->getVals();
    data_t * pm0 = m0->getVals();
    data_t * pp = p->getVals();
    data_t alpha=0,df1=0,f=0;
    data_t alo=_stp0, ahi=_stp1,flo=_f0,fhi=_f1;
    if (reverse) {alo=_stp1, ahi=_stp0,flo=_f1,fhi=_f0;}

    while (true){
        if (trial >= max_trial){
            fprintf(stderr,"==============================\nWARNING: Maximum number of trials reached before finding a satisfactory step length. Solver STOPS\n==============================\n");
            return false;
        }
        // Interpolate quadratic or cubic between stp0 and stp1
        if (alo==0) alpha =  quadratic(_f,_df,ahi,fhi);
        else if (ahi==0) alpha =  quadratic(_f,_df,alo,flo);
        else alpha = cubic(_f,_df,alo,flo,ahi,fhi);

        // Evaluate phi(alpha)
        #pragma omp parallel for
        for (int i=0; i<n; i++) pm[i] = pm0[i] + alpha*pp[i];
        prob->res();
        feval++;
        trial++;
        f = prob->getFunc();
        _f0=_f1; _f1= f;
        if (f!=f) throw std::runtime_error("==============================\nERROR: the objective function is NaN.\n==============================\n");

        if (f<=ZERO){
            _stp=alpha;
            _f0=_f;
            _f = f;
            if (verbose){
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Trial = %d\n",trial);
                    fprintf(stderr,"Step length = %f\n",alpha);
                    fprintf(stderr,"Function value = %f\n",f);
                    fprintf(stderr,"Last gradient norm = %f\n",gnorm);
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Found a satisfactory step length\n");
                    fprintf(stderr,"Function value became too small\n");
                    fprintf(stderr,"==============================\n");
            }
            return false;
        }
        // check if (alo, alpha) contains the desired step length
        if ((_f1>_f+_c1*alpha*_df) || (_f1>=flo)){
            if (verbose){
                fprintf(stderr,"==============================\n");
                fprintf(stderr,"Trial = %d\n",trial);
                fprintf(stderr,"Step length = %f\n",alpha);
                fprintf(stderr,"Function value = %f\n",f);
                fprintf(stderr,"Last gradient norm = %f\n",gnorm);
            }
            ahi = alpha;
            fhi = _f1;
        }
        else{
            // Evaluate phi'(alpha)
            prob->grad();
            geval++;
            gnorm=g->norm();
            df1 = g->dot(p);

            if (verbose){
                fprintf(stderr,"==============================\n");
                fprintf(stderr,"Trial = %d\n",trial);
                fprintf(stderr,"Step length = %f\n",alpha);
                fprintf(stderr,"Function value = %f\n",f);
                fprintf(stderr,"Gradient norm = %f\n",gnorm);
            }

            // check the curvature condition
            if (df1 >= _c2*_df){
                _stp=alpha;
                _f0=_f;
                _f = f;
                if (verbose){
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Found a satisfactory step length\n");
                    fprintf(stderr,"==============================\n");
                }
                return true;
            }

            if (df1*(ahi-alo)>=0) {
                ahi = alo;
                fhi = flo;
            }
            alo = alpha;
            flo = _f1;
        }
    }
}

bool strong_wolfe::lineSearch(optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose){

    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();
    int n = m->getN123();
    data_t * pm = m->getVals();
    data_t * pm0 = m0->getVals();
    data_t * pp = p->getVals();
    data_t alpha=0,df0=0,df1=0,f=0;

    int trial=0;
    // first trial
    df0 = _df;
    _f = prob->getFunc();
    _df = g->dot(p);
    f=_f;
    _stp0=0;

    // first iteration
    if (iter == 0) {
        _f0=_f;
        _stp1 = _a0/gnorm;
        if (_a1>0) _stp1 = std::max(_stp1,_a1*m->norm()/gnorm);
        alpha=_stp1;
    }

    // subsequent iterations
    else {
        _stp1 = _stp*df0/_df;
        if (_version==1){
            _stp1 = 2*(_f - _f0)/_df;
            _stp1 = std::min(1.0,1.00001*_stp1);
        }
        alpha=_stp1;
    }

    while (true){
        if (trial >= max_trial){
            fprintf(stderr,"==============================\nWARNING: Maximum number of trials reached before finding a satisfactory step length. Solver STOPS\n==============================\n");
            return false;
        }

        // Evaluate phi(alpha)
        #pragma omp parallel for
        for (int i=0; i<n; i++) pm[i] = pm0[i] + alpha*pp[i];
        prob->res();
        feval++;
        trial++;
        f = prob->getFunc();
        _f1= f;
        if (f!=f) throw std::runtime_error("==============================\nERROR: the objective function is NaN.\n==============================\n");

        if (f<=ZERO){
            _stp=alpha;
            _f0=_f;
            _f = f;
            if (verbose){
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Trial = %d\n",trial);
                    fprintf(stderr,"Step length = %f\n",alpha);
                    fprintf(stderr,"Function value = %f\n",f);
                    fprintf(stderr,"Last gradient norm = %f\n",gnorm);
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Found a satisfactory step length\n");
                    fprintf(stderr,"Function value became too small\n");
                    fprintf(stderr,"==============================\n");
            }
            return false;
        }
        // check if (_stp0, _stp1) contains the desired step length
        if ((_f1>_f+_c1*alpha*_df) || (_f1>=_f0 && trial>1)){
            if (verbose){
                fprintf(stderr,"==============================\n");
                fprintf(stderr,"Trial = %d\n",trial);
                fprintf(stderr,"Step length = %f\n",alpha);
                fprintf(stderr,"Function value = %f\n",f);
                fprintf(stderr,"Last gradient norm = %f\n",gnorm);
            }
            return zoom(false,prob,m0,p,iter,trial,max_trial,feval,geval,gnorm,verbose);
        }

        // Evaluate phi'(stp1)
        prob->grad();
        geval++;
        gnorm=g->norm();
        df1 = g->dot(p);

        if (verbose){
                fprintf(stderr,"==============================\n");
                fprintf(stderr,"Trial = %d\n",trial);
                fprintf(stderr,"Step length = %f\n",alpha);
                fprintf(stderr,"Function value = %f\n",f);
                fprintf(stderr,"Gradient norm = %f\n",gnorm);
        }

        // check the curvature condition
        if (std::abs(df1)<= -_c2*_df){
            _stp=alpha;
            _f0=_f;
            _f = f;
            if (verbose){
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Found a satisfactory step length\n");
                    fprintf(stderr,"==============================\n");
            }
            return true;
        }

        // Check if phi(stp1) is increasing
        if (df1 >= 0){
            return zoom(true,prob,m0,p,iter,trial,max_trial,feval,geval,gnorm,verbose);
        }

        // Choose another step in (stp1,stp_max)
        _stp0 = _stp1;
        _stp1 = quadratic(_f,_df,_stp1,_f1);
        _stp1 = std::min(_stp1, _stp_max);
        alpha = _stp1;
        _f0=_f1;
    }   
}

bool strong_wolfe::zoom(bool reverse,
                        optimization * prob,
                        std::shared_ptr<vecReg<data_t> >m0,
                        std::shared_ptr<vecReg<data_t> >p,
                        int iter,
                        int &trial,
                        int max_trial,
                        int &feval,
                        int &geval,
                        data_t &gnorm,
                        bool verbose){

    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();
    int n = m->getN123();
    data_t * pm = m->getVals();
    data_t * pm0 = m0->getVals();
    data_t * pp = p->getVals();
    data_t alpha=0,df1=0,f=0;
    data_t alo=_stp0, ahi=_stp1,flo=_f0,fhi=_f1;
    if (reverse) {alo=_stp1, ahi=_stp0,flo=_f1,fhi=_f0;}

    while (true){
        if (trial >= max_trial){
            fprintf(stderr,"==============================\nWARNING: Maximum number of trials reached before finding a satisfactory step length. Solver STOPS\n==============================\n");
            return false;
        }
        // Interpolate quadratic or cubic between stp0 and stp1
        if (alo==0) alpha =  quadratic(_f,_df,ahi,fhi);
        else if (ahi==0) alpha =  quadratic(_f,_df,alo,flo);
        else alpha = cubic(_f,_df,alo,flo,ahi,fhi);

        // Evaluate phi(alpha)
        #pragma omp parallel for
        for (int i=0; i<n; i++) pm[i] = pm0[i] + alpha*pp[i];
        prob->res();
        feval++;
        trial++;
        f = prob->getFunc();
        _f0=_f1; _f1= f;
        if (f!=f) throw std::runtime_error("==============================\nERROR: the objective function is NaN.\n==============================\n");

        if (f<=ZERO){
            _stp=alpha;
            _f0=_f;
            _f = f;
            if (verbose){
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Trial = %d\n",trial);
                    fprintf(stderr,"Step length = %f\n",alpha);
                    fprintf(stderr,"Function value = %f\n",f);
                    fprintf(stderr,"Last gradient norm = %f\n",gnorm);
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Found a satisfactory step length\n");
                    fprintf(stderr,"Function value became too small\n");
                    fprintf(stderr,"==============================\n");
            }
            return false;
        }
        // check if (alo, alpha) contains the desired step length
        if ((_f1>_f+_c1*alpha*_df) || (_f1>=flo)){
            if (verbose){
                fprintf(stderr,"==============================\n");
                fprintf(stderr,"Trial = %d\n",trial);
                fprintf(stderr,"Step length = %f\n",alpha);
                fprintf(stderr,"Function value = %f\n",f);
                fprintf(stderr,"Last gradient norm = %f\n",gnorm);
            }
            ahi = alpha;
            fhi = _f1;
        }
        else{
            // Evaluate phi'(alpha)
            prob->grad();
            geval++;
            gnorm=g->norm();
            df1 = g->dot(p);

            if (verbose){
                fprintf(stderr,"==============================\n");
                fprintf(stderr,"Trial = %d\n",trial);
                fprintf(stderr,"Step length = %f\n",alpha);
                fprintf(stderr,"Function value = %f\n",f);
                fprintf(stderr,"Gradient norm = %f\n",gnorm);
            }

            // check the curvature condition
            if (std::abs(df1)<= -_c2*_df){
                _stp=alpha;
                _f0=_f;
                _f = f;
                if (verbose){
                    fprintf(stderr,"==============================\n");
                    fprintf(stderr,"Found a satisfactory step length\n");
                    fprintf(stderr,"==============================\n");
                }
                return true;
            }

            if (df1*(ahi-alo)>=0) {
                ahi = alo;
                fhi = flo;
            }
            alo = alpha;
            flo = _f1;
        }
    }
}

void nlsolver::testLinear(const bool verbose){

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

    this->run(&prob,verbose);

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
    std::clog << "model error norm squared = " << error << std::endl;
}

void nlsolver::testRosenbrock(const bool verbose){

    std::shared_ptr<vecReg<data_t> > m = std::make_shared<vecReg<data_t> >(hypercube<data_t>(2));
    data_t x0 = -1.2;
    data_t y0 = 1.0;
    m->getVals()[0]=x0;
    m->getVals()[1]=y0;
    data_t a=1, b=100;
    rosenbrock ros(a,b,m);
    this->run(&ros,verbose);
    data_t err = (m->getVals()[0]-a)*(m->getVals()[0]-a) + (m->getVals()[1]-a*a)*(m->getVals()[1]-a*a);
    err = sqrt(err);

    fprintf(stderr,"True model x=%f ; y=%f\n",a,a*a);
    fprintf(stderr,"Initial guess x=%f ; y=%f\n",x0,y0);
    fprintf(stderr,"Estimated model x=%f ; y=%f\n",m->getVals()[0],m->getVals()[1]);
    fprintf(stderr,"Model error norm err=%f\n",err);
    fprintf(stderr,"Initial functional f=%f\n",_func[0]);
    fprintf(stderr,"Final functional f=%f\n",_func[_func.size()-1]);
    fprintf(stderr,"Last gradient norm g=%f\n",ros.getGrad()->norm());
    fprintf(stderr,"Number of iterations niter=%d\n",_func.size()-1);
}

void nlsolver::testParaboloid(const bool verbose){

    std::shared_ptr<vecReg<data_t> > m = std::make_shared<vecReg<data_t> >(hypercube<data_t>(2));
    data_t x0 = -1;
    data_t y0 = 2;
    m->getVals()[0]=x0;
    m->getVals()[1]=y0;
    data_t a=1, b=1;
    paraboloid para(a,b,m);
    this->run(&para,verbose);

    fprintf(stderr,"True model x=%f ; y=%f\n",0,0);
    fprintf(stderr,"Initial guess x=%f ; y=%f\n",x0,y0);
    fprintf(stderr,"Estimated model x=%f ; y=%f\n",m->getVals()[0],m->getVals()[1]);
    fprintf(stderr,"Model error norm err=%f\n",m->norm());
    fprintf(stderr,"Initial functional f=%f\n",_func[0]);
    fprintf(stderr,"Final functional f=%f\n",_func[_func.size()-1]);
    fprintf(stderr,"Last gradient norm g=%f\n",para.getGrad()->norm());
    fprintf(stderr,"Number of iterations niter=%d\n",_func.size()-1);
}

void nlsd::run(optimization * prob, const bool verbose, std::string output, int isave){
    
    reset();

    prob->initGrad();
    prob->initRes();

    //residuals
    prob->res();

    //model
    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    std::shared_ptr<vecReg<data_t> > m0 = m->clone();
    int n = m->getN123();
    data_t * pm = m->getVals();
    data_t * pm0 = m0->getVals();
    
    //gradient and search direction
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();
    std::shared_ptr<vecReg<data_t> > p = g->clone();
    data_t * pp = p->getVals();

    //objective function
    _feval++;
    data_t f = prob->getFunc();
    _func.push_back(f);

    if ((_func[0] <= ZERO) || (_func[0]!= _func[0])) {
        fprintf(stderr,"==============================\nERROR: the initial objective function is negative or NaN.\n==============================\n");
        return;
    }

    //iteration number
    int k = 0;

    //convergence rate
    data_t rate = 1;
    data_t gnorm = 999;

fprintf(stderr,"#########################################################################\n");
fprintf(stderr,"Iteration = 0; functional = %f; normalized functional = 1\n",f);
fprintf(stderr,"#########################################################################\n");

    //start the SD loop
    bool success = true;
    while (k<_niter && rate>_threshold && gnorm>ZERO && f>ZERO && success){

        // compute the gradient
        if (_lsearch->getFlagG() || k==0){
            prob->grad();
            _geval++;
            gnorm=g->norm();
            if (output!="none" && k==0 && isave!=0){
                sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
                sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
            }
        }

        // update search direction (negative of the gradient)
        p->scaleAdd(g,0,-1);

        // copy the updated model to the backup vector
        if (k>0) memcpy(pm0,pm,n*sizeof(data_t));

        // perform a line search
        success = _lsearch->lineSearch(prob,m0,p,k,_max_trial,_feval,_geval,gnorm,verbose);
        if (!success) break;
        
        f = prob->getFunc();
        _func.push_back(f);

fprintf(stderr,"#########################################################################\n");
fprintf(stderr,"Iteration = %d; functional = %f; normalized functional = %f\n",k+1,f,f / _func[0]);
fprintf(stderr,"#########################################################################\n");

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
        }

        //iterate
        k++;

        if (output!="none" && k % isave==0){
            sepWrite(m,output+"model_iter_"+std::to_string(k)+".H");
            sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
            sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
        }
    }
    if (output!="none" && k % isave!=0){
        sepWrite(m,output+"model_iter_"+std::to_string(k)+".H");
        sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
        sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
    }

    fprintf(stderr,"\n==============================\n");
    fprintf(stderr,"Total number of NLSD iterations: %d\n",k);
    fprintf(stderr,"Total number of function evaluations: %d\n",_feval);
    fprintf(stderr,"Total number of gradient evaluations: %d\n",_geval);
    fprintf(stderr,"==============================\n\n");
}

void nlcg::run(optimization * prob, const bool verbose, std::string output, int isave){
    
    reset();

    prob->initGrad();
    prob->initRes();

    //residuals
    prob->res();

    //model
    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    std::shared_ptr<vecReg<data_t> > m0 = m->clone();
    int n = m->getN123();
    data_t * pm = m->getVals();
    data_t * pm0 = m0->getVals();
    
    //gradient and search direction
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();
    std::shared_ptr<vecReg<data_t> > p = g->clone();
    std::shared_ptr<vecReg<data_t> > g0;
    data_t * pg0;
    if (_method>0) {g0 = g->clone(); pg0 = g0->getVals();}
    data_t * pg = g->getVals();
    data_t * pp = p->getVals();

    //objective function
    _feval++;
    data_t f = prob->getFunc();
    _func.push_back(f);

    if ((_func[0] <= ZERO) || (_func[0] != _func[0])) {
        fprintf(stderr,"==============================\nERROR: the initial objective function is negative or NaN.\n==============================\n");
        return;
    }

    //direction weight
    data_t beta=0;

    //iteration number
    int k = 0;

    //convergence rate
    data_t rate = 1;
    data_t gnorm = 999;
    data_t gnorm0 = -1;

fprintf(stderr,"#########################################################################\n");
fprintf(stderr,"Iteration = 0; functional = %f; normalized functional = 1\n",f);
fprintf(stderr,"#########################################################################\n");

    //start the CG loop
    bool success=true;
    while (k<_niter && rate>_threshold && gnorm>ZERO && f>ZERO && success){

        // compute the gradient
        if (_lsearch->getFlagG() || k==0){
            prob->grad();
            _geval++;
            gnorm=g->norm();
            if (output!="none" && k==0 && isave!=0){
                sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
                sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
            }
        }

        // compute beta
        if (k>0){
            if (_method == 0) beta = gnorm*gnorm/(gnorm0*gnorm0);
            else if (_method == 1) beta = (gnorm*gnorm - g->dot(g0))/(gnorm0*gnorm0);
            else if (_method == 2) beta = std::max((data_t)0.0, (gnorm*gnorm - g->dot(g0))/(gnorm0*gnorm0));
        }

        // update search direction
        p->scaleAdd(g,beta,-1);

        // check that p is a descent direction
        if (p->dot(g) >= 0) {
            fprintf(stderr,"==============================\nWARNING: Search direction is not a descent direction. Solver STOPS\n==============================\n");
            break;            
        }

        // copy the updated model to the backup vector
        if (k>0) memcpy(pm0,pm,n*sizeof(data_t));

        // store the previous gradient
        if (_method>0) memcpy(pg0,pg,n*sizeof(data_t));
        gnorm0=gnorm;

        // perform a line search
        success=_lsearch->lineSearch(prob,m0,p,k,_max_trial,_feval,_geval,gnorm,verbose);
        if (!success) break;

        f = prob->getFunc();
        _func.push_back(f);

fprintf(stderr,"#########################################################################\n");
fprintf(stderr,"Iteration = %d; functional = %f; normalized functional = %f\n",k+1,f,f / _func[0]);
fprintf(stderr,"#########################################################################\n");

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
        }

        //iterate
        k++;

        if (output!="none" && k % isave==0){
            sepWrite(m,output+"model_iter_"+std::to_string(k)+".H");
            sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
            sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
        }
    }
    if (output!="none" && k % isave!=0){
        sepWrite(m,output+"model_iter_"+std::to_string(k)+".H");
        sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
        sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
    }

    fprintf(stderr,"\n==============================\n");
    fprintf(stderr,"Total number of NLCG iterations: %d\n",k);
    fprintf(stderr,"Total number of function evaluations: %d\n",_feval);
    fprintf(stderr,"Total number of gradient evaluations: %d\n",_geval);
    fprintf(stderr,"==============================\n\n");
}

void bfgs::run(optimization * prob, const bool verbose, std::string output, int isave){
    
    reset();

    prob->initGrad();
    prob->initRes();

    //model
    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    std::shared_ptr<vecReg<data_t> > m0 = m->clone();
    int n = m->getN123();
    data_t * pm = m->getVals();
    data_t * pm0 = m0->getVals();
    
    //gradient and search direction
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();
    std::shared_ptr<vecReg<data_t> > p = g->clone();
    std::shared_ptr<vecReg<data_t> > g0 = g->clone();
    data_t * pg = g->getVals();
    data_t * pg0 = g0->getVals();
    data_t * pp = p->getVals();

    //residuals
    prob->res();

    //objective function
    _feval++;
    data_t f = prob->getFunc();
    _func.push_back(f);

    if ((_func[0] <= ZERO) || (_func[0] != _func[0])) {
        fprintf(stderr,"==============================\nERROR: the initial objective function is negative or NaN.\n==============================\n");
        return;
    }

    //initialize the inverse hessian H
    if (_H == nullptr) _H = new matrix(n,1);
    else {
        successCheck(sqrt(_H->getMat()->getN123())==n,__FILE__,__LINE__,"The provided inverse Hessian has an incorrect size\n");
    }

    //model and gradient difference vectors for inverse Hessian update
    std::shared_ptr<vecReg<data_t> > s = m->clone();
    std::shared_ptr<vecReg<data_t> > y = g->clone();
    data_t * ps = s->getVals();
    data_t * py = y->getVals();

    //iteration number
    int k = 0;

    //convergence rate
    data_t rate = 1;
    data_t gnorm = 999;

fprintf(stderr,"#########################################################################\n");
fprintf(stderr,"Iteration = 0; functional = %f; normalized functional = 1\n",f);
fprintf(stderr,"#########################################################################\n");

    //start the BFGS loop
    bool success = true;
    while (k<_niter && rate>_threshold && gnorm>ZERO && f>ZERO && success){

        // compute the gradient
        if (_lsearch->getFlagG() || k==0){
            prob->grad();
            _geval++;
            gnorm=g->norm();
            if (output!="none" && k==0 && isave!=0){
                sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
                sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
            }
        }

        //update the inverse Hessian
        if (k>0) {
            #pragma omp parallel for
            for (int i=0; i<n; i++){
                ps[i] = pm[i] - pm0[i];
                py[i] = pg[i] - pg0[i];
            }
            updateH(_H,s,y);
        }

        // update search direction
        _H->apply_forward(false,pg,pp);
        p->scale(-1);

        // check that p is a descent direction
        if (p->dot(g) >= 0) {
            fprintf(stderr,"==============================\nWARNING: Search direction is not a descent direction. Solver STOPS\n==============================\n");
            break;            
        }

        // copy the updated model to the backup vector
        if (k>0) memcpy(pm0,pm,n*sizeof(data_t));

        // store the previous gradient
        memcpy(pg0,pg,n*sizeof(data_t));

        // perform a line search
        success = _lsearch->lineSearch(prob,m0,p,k,_max_trial,_feval,_geval,gnorm,verbose);
        if (!success) break;
        
        f = prob->getFunc();
        _func.push_back(f);

fprintf(stderr,"#########################################################################\n");
fprintf(stderr,"Iteration = %d; functional = %f; normalized functional = %f\n",k+1,f,f / _func[0]);
fprintf(stderr,"#########################################################################\n");

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
        }

        //iterate
        k++;

        if (output!="none" && k % isave==0){
            sepWrite(m,output+"model_iter_"+std::to_string(k)+".H");
            sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
            sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
        }
    }
    if (output!="none" && k % isave!=0){
        sepWrite(m,output+"model_iter_"+std::to_string(k)+".H");
        sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
        sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
    }

    fprintf(stderr,"\n==============================\n");
    fprintf(stderr,"Total number of BFGS iterations: %d\n",k);
    fprintf(stderr,"Total number of function evaluations: %d\n",_feval);
    fprintf(stderr,"Total number of gradient evaluations: %d\n",_geval);
    fprintf(stderr,"==============================\n\n");
}

bool bfgs::testSPD(matrix * M){

    std::shared_ptr<vecReg<data_t> > x = std::make_shared<vecReg<data_t> > (*M->getDomain());
    std::shared_ptr<vecReg<data_t> > y = std::make_shared<vecReg<data_t> > (*M->getDomain());
    x->zero();
    y->zero();
    data_t * pmat = M->getMat()->getVals();

    bool ans = true;
    int n = x->getN123();

    int i=0; int j=1;
    while (ans==true && i<n){
        j=i+1;
        while(ans==true && j<n){
            if (std::abs(pmat[i*n+j] - pmat[j*n+i])>1e-7) ans=false;
            j++;
        }
        i++;
    }

    i=0;
    while (ans == true && i<100){
        x->random();
        M->forward(false,x,y);
        if (x->dot(y)<=0) ans=false;
        i++;
    }

    if (ans==false){
        fprintf(stderr,"\n\n\n WARNING: the matrix H is not SPD!!! \n\n\n");
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                std::clog << pmat[i*n+j] << "\t";
            }
            std::clog << "\n";
        }   
    }

    return ans;
}

void bfgs::updateH(matrix * H ,std::shared_ptr<vecReg<data_t> > s, std::shared_ptr<vecReg<data_t> > y){
    
    int n = s->getN123();

    std::shared_ptr<vecReg<data_t> > x = s->clone();
    data_t * ps = s->getVals();
    data_t * py = y->getVals();
    data_t * px = x->getVals();
    data_t * pmat = H->getMat()->getVals();

    H->apply_forward(false,py,px);
    data_t rho = 1.0/s->dot(y);
    successCheck(rho > 0,__FILE__,__LINE__,"The value of rho=1/(s'y) must be > 0 when updating the Hessian matrix\n");
    data_t a = y->dot(x);
    data_t ui;

    #pragma omp parallel for private(ui)
    for (int i=0; i<n; i++){
        ui = rho*((1+rho*a)*ps[i] - px[i]);
        for (int j=i; j<n; j++){
            pmat[i*n+j] = pmat[i*n+j] + ui*ps[j] - rho*ps[i]*px[j];
            pmat[j*n+i] = pmat[i*n+j];
        }
    }
    //testSPD(H);
}

void lbfgs::run(optimization * prob, const bool verbose, std::string output, int isave){
    
    reset();

    prob->initGrad();
    prob->initRes();

    //model
    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    std::shared_ptr<vecReg<data_t> > m0 = m->clone();
    int n = m->getN123();
    data_t * pm = m->getVals();
    data_t * pm0 = m0->getVals();
    
    //gradient and search direction
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();
    std::shared_ptr<vecReg<data_t> > p = g->clone();
    std::shared_ptr<vecReg<data_t> > g0 = g->clone();
    data_t * pg = g->getVals();
    data_t * pg0 = g0->getVals();
    data_t * pp = p->getVals();

    //residuals
    prob->res();

    //objective function
    _feval++;
    data_t f = prob->getFunc();
    _func.push_back(f);

    if ((_func[0] <= ZERO) || (_func[0] != _func[0])) {
        fprintf(stderr,"==============================\nERROR: the initial objective function is negative or NaN.\n==============================\n");
        return;
    }

    //initialize the s,y pairs and the inverse hessian H
    if (_H0 == nullptr) {_H0 = std::make_shared<vecReg<data_t> >(*m->getHyper()); _H0->set(1);}
    else successCheck(_H0->getN123()==n,__FILE__,__LINE__,"The provided inverse Hessian has an incorrect size\n");
    std::vector<std::shared_ptr<vecReg<data_t> > > s(_m);
    std::vector<std::shared_ptr<vecReg<data_t> > > y(_m);
    std::vector<data_t> rho(_m);
    for (int i=0; i<_m; i++){
        s[i] = std::make_shared<vecReg<data_t> >(*m->getHyper());
        y[i] = std::make_shared<vecReg<data_t> >(*g->getHyper());
    }

    //iteration number
    int k = 0;

    //convergence rate
    data_t rate = 1;
    data_t gnorm = 999;

fprintf(stderr,"#########################################################################\n");
fprintf(stderr,"Iteration = 0; functional = %f; normalized functional = 1\n",f);
fprintf(stderr,"#########################################################################\n");

    //start the l-BFGS loop
    bool success = true;
    while (k<_niter && rate>_threshold && gnorm>ZERO && f>ZERO && success){

        // compute the gradient
        if (_lsearch->getFlagG() || k==0){
            prob->grad();
            _geval++;
            gnorm=g->norm();
            if (output!="none" && k==0 && isave!=0){
                sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
                sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
            }
        }

        //update the inverse Hessian s,y pairs
        updateSY(s,y,rho,pm,pm0,pg,pg0,k);

        // update search direction
        computeHg(p,g,s,y,rho,k);
        p->scale(-1);

        // check that p is a descent direction
        if (p->dot(g) >= 0) {
            fprintf(stderr,"==============================\nWARNING: Search direction is not a descent direction. Solver STOPS\n==============================\n");
            break;            
        }

        // copy the updated model to the backup vector
        if (k>0) memcpy(pm0,pm,n*sizeof(data_t));

        // store the previous gradient
        memcpy(pg0,pg,n*sizeof(data_t));

        // perform a line search
        success = _lsearch->lineSearch(prob,m0,p,k,_max_trial,_feval,_geval,gnorm,verbose);
        if (!success) break;

        f = prob->getFunc();
        _func.push_back(f);

fprintf(stderr,"#########################################################################\n");
fprintf(stderr,"Iteration = %d; functional = %f; normalized functional = %f\n",k+1,f,f / _func[0]);
fprintf(stderr,"#########################################################################\n");

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
        }

        //iterate
        k++;

        if (output!="none" && k % isave==0){
            sepWrite(m,output+"model_iter_"+std::to_string(k)+".H");
            sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
            sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
        }
    }
    if (output!="none" && k % isave!=0){
        sepWrite(m,output+"model_iter_"+std::to_string(k)+".H");
        sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
        sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
    }

    fprintf(stderr,"\n==============================\n");
    fprintf(stderr,"Total number of l-BFGS iterations: %d\n",k);
    fprintf(stderr,"Total number of function evaluations: %d\n",_feval);
    fprintf(stderr,"Total number of gradient evaluations: %d\n",_geval);
    fprintf(stderr,"==============================\n\n");
}

void lbfgs::updateSY(std::vector<std::shared_ptr<vecReg<data_t> > > &s,
                    std::vector<std::shared_ptr<vecReg<data_t> > > &y,
                    std::vector<data_t> &rho,
                    data_t * pm, data_t * pm0,
                    data_t * pg, data_t * pg0,
                    int k){

    if (k > 0 && k<=_m){
        data_t * ps = s[k-1]->getVals();
        data_t * py = y[k-1]->getVals();
        int n = s[0]->getN123();
        #pragma omp parallel for
        for (int i=0; i<n; i++){
            ps[i] = pm[i] - pm0[i];
            py[i] = pg[i] - pg0[i];
        }
        rho[k-1] = s[k-1]->dot(y[k-1]);
    }
    else if (k>_m){
        std::shared_ptr<vecReg<data_t> > temp_s = s[0];
        std::shared_ptr<vecReg<data_t> > temp_y = y[0];
        for (int i=0; i<_m-1; i++){
            s[i] = s[i+1];
            y[i] = y[i+1];
            rho[i] = rho[i+1];
        }
        s[_m-1] = temp_s;
        y[_m-1] = temp_y;
        data_t * ps = s[_m-1]->getVals();
        data_t * py = y[_m-1]->getVals();
        int n = s[0]->getN123();
        #pragma omp parallel for
        for (int i=0; i<n; i++){
            ps[i] = pm[i] - pm0[i];
            py[i] = pg[i] - pg0[i];
        }
        rho[_m-1] = s[_m-1]->dot(y[_m-1]);
    }
}

void lbfgs::computeHg(std::shared_ptr<vecReg<data_t> > &p, std::shared_ptr<vecReg<data_t> > g,
                                std::vector<std::shared_ptr<vecReg<data_t> > > s, std::vector<std::shared_ptr<vecReg<data_t> > > y, std::vector<data_t> rho, int k){

    std::shared_ptr<vecReg<data_t> > q = g->clone();
    int m = std::min(k,_m);
    data_t b, gamma=1.0;
    std::vector<data_t> a(m,0);
    data_t * pH0 = _H0->getVals();
    data_t * pp = p->getVals();
    data_t * pq = q->getVals();
    int n = p->getN123();

    if (m>0){
        for (int i=m-1; i>=0; i--){
            a[i] = s[i]->dot(q)/rho[i];
            q->scaleAdd(y[i],1,-a[i]);
        }
        gamma = rho[m-1]/y[m-1]->norm2();
    }
    #pragma omp parallel for
    for (int i=0; i<n; i++) pp[i] = gamma*pH0[i]*pq[i];
    if (m>0){
        for (int i=0; i<m; i++){
            b = (y[i]->dot(p))/rho[i];
            p->scaleAdd(s[i],1,a[i]-b);
        }
    }
}

void newton::run(optimization * prob, const bool verbose, std::string output, int isave){
    
    reset();

    prob->initGrad();
    prob->initRes();

    //model
    std::shared_ptr<vecReg<data_t> > m = prob->getMod();
    std::shared_ptr<vecReg<data_t> > m0 = m->clone();
    int n = m->getN123();
    data_t * pm = m->getVals();
    data_t * pm0 = m0->getVals();
    
    //gradient and search direction
    std::shared_ptr<vecReg<data_t> > g = prob->getGrad();
    std::shared_ptr<vecReg<data_t> > p = g->clone();
    std::shared_ptr<vecReg<data_t> > g0 = g->clone();
    data_t * pg = g->getVals();
    data_t * pg0 = g0->getVals();
    data_t * pp = p->getVals();

    //residuals
    prob->res();

    //objective function
    _feval++;
    data_t f = prob->getFunc();
    _func.push_back(f);

    if ((_func[0] <= ZERO) || (_func[0] != _func[0])) {
        fprintf(stderr,"==============================\nERROR: the initial objective function is negative or NaN.\n==============================\n");
        return;
    }

    //initialize the hessian matrix B
    std::shared_ptr<vecReg<data_t> > B = std::make_shared<vecReg<data_t> >(hypercube<data_t>(n*n));

    //iteration number
    int k = 0;

    //convergence rate
    data_t rate = 1;
    data_t gnorm = 999;

fprintf(stderr,"#########################################################################\n");
fprintf(stderr,"Iteration = 0; functional = %f; normalized functional = 1\n",f);
fprintf(stderr,"#########################################################################\n");

    //start the Newton loop
    bool success = true;
    while (k<_niter && rate>_threshold && gnorm>ZERO && f>ZERO && success){

        // compute the gradient
        if (_lsearch->getFlagG() || k==0){
            prob->grad();
            _geval++;
            gnorm=g->norm();
            if (output!="none" && k==0 && isave!=0){
                sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
                sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
            }
        }

        //compute the Hessian
        prob->hessian(B);

        // update search direction s.t. Bp = g
        solve_Axb<data_t>(B->getVals(), pp, pg, n);
        p->scale(-1);

        // check that p is a descent direction
        if (p->dot(g) >= 0) {
            fprintf(stderr,"==============================\nWARNING: Search direction is not a descent direction. Solver STOPS\n==============================\n");
            break;            
        }

        // copy the updated model to the backup vector
        if (k>0) memcpy(pm0,pm,n*sizeof(data_t));

        // store the previous gradient
        memcpy(pg0,pg,n*sizeof(data_t));

        // perform a line search
        success = _lsearch->lineSearch(prob,m0,p,k,_max_trial,_feval,_geval,gnorm,verbose);
        if (!success) break;
        
        f = prob->getFunc();
        _func.push_back(f);

fprintf(stderr,"#########################################################################\n");
fprintf(stderr,"Iteration = %d; functional = %f; normalized functional = %f\n",k+1,f,f / _func[0]);
fprintf(stderr,"#########################################################################\n");

        //compute the rate of convergence
        if(k>0){
            rate = (_func[k]-_func[k+1])/_func[k];
        }

        //iterate
        k++;

        if (output!="none" && k % isave==0){
            sepWrite(m,output+"model_iter_"+std::to_string(k)+".H");
            sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
            sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
        }
    }
    if (output!="none" && k % isave!=0){
        sepWrite(m,output+"model_iter_"+std::to_string(k)+".H");
        sepWrite(g,output+"gradient_iter_"+std::to_string(k)+".H");
        sepWrite(prob->getRes(),output+"residual_iter_"+std::to_string(k)+".H");
    }

    fprintf(stderr,"\n==============================\n");
    fprintf(stderr,"Total number of Newton iterations: %d\n",k);
    fprintf(stderr,"Total number of function evaluations: %d\n",_feval);
    fprintf(stderr,"Total number of gradient evaluations: %d\n",_geval);
    fprintf(stderr,"==============================\n\n");
}


#undef ZERO
#undef STOL