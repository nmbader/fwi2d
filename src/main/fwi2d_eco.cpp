#ifdef DOUBLE_PRECISION
    typedef double data_t;
#else
    typedef float data_t;
#endif

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

// Executable to run 2D FWI

int main(int argc, char **argv){
    
    initpar(argc,argv);

// Read parameters for wave propagation and inversion
    param par;
    readParameters(argc, argv, par);

// Read inputs/outputs files
    std::string source_file="none", model_file="none", data_file="none", output_file="none", ioutput_file="none", obj_func_file="none", inverse_diagonal_hessian_file="none";
    readParam<std::string>(argc, argv, "source", source_file);
    readParam<std::string>(argc, argv, "model", model_file);
    readParam<std::string>(argc, argv, "data", data_file);
    readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "ioutput", ioutput_file);
    readParam<std::string>(argc, argv, "obj_func", obj_func_file);

    successCheck(source_file!="none",__FILE__,__LINE__,"Source wavelet is not provided\n");
    successCheck(model_file!="none",__FILE__,__LINE__,"Earth model is not provided\n");
    successCheck(data_file!="none",__FILE__,__LINE__,"Data to be inverted is not provided\n");

    std::shared_ptr<vec> src = sepRead<data_t>(source_file);
    std::shared_ptr<vec> data = sepRead<data_t>(data_file);
    std::shared_ptr<vec> model = sepRead<data_t>(model_file);

    std::shared_ptr<vec> gmask = nullptr;
    std::shared_ptr<vec> w = nullptr;
    std::shared_ptr<vec> invDiagH = nullptr;
    std::shared_ptr<vec> prior = nullptr;
    if (par.mask_file!="none") gmask = sepRead<data_t>(par.mask_file);
    if (par.weights_file!="none") w = sepRead<data_t>(par.weights_file);
    if (par.inverse_diagonal_hessian_file!="none") invDiagH = sepRead<data_t>(par.inverse_diagonal_hessian_file);
    if (par.prior_file!="none") prior = sepRead<data_t>(par.prior_file);

// Analyze the inputs and parameters and modify if necessary
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0);
    analyzeBsplines(*model->getHyper(),par);
    analyzeNLInversion(par);

// Build model precon if B-splines are activated
// ----------------------------------------------------------------------------------------//
    std::shared_ptr<vecReg<data_t> > bsmodel = model;
    std::shared_ptr<vecReg<data_t> > bsmask = gmask;
    std::shared_ptr<vecReg<data_t> > bsinvDiagH = invDiagH;
    std::shared_ptr<vecReg<data_t> > bsprior = prior;
    chainLOper * BD;

if (par.bsplines)
{
    std::vector<data_t> kx;
    std::vector<data_t> kz;
    setKnot(kx,par.bs_controlx,par.bs_mx);
    setKnot(kz,par.bs_controlz,par.bs_mz);

    std::vector<axis<data_t> > axes = model->getHyper()->getAxes();
    axes[0].n=par.bs_nz; axes[1].n=par.bs_nx;
    bsmodel = std::make_shared<vecReg<data_t> >(vecReg<data_t>(hypercube<data_t>(axes)));
    fillin(bsmodel,model,par.bs_controlx,par.bs_controlz);
    if (par.mask_file != "none"){
        bsmask = std::make_shared<vecReg<data_t> >(vecReg<data_t>(hypercube<data_t>(axes)));
        fillin(bsmask,gmask,par.bs_controlx,par.bs_controlz);
    }
    if (par.inverse_diagonal_hessian_file != "none"){
        bsinvDiagH = std::make_shared<vecReg<data_t> >(vecReg<data_t>(hypercube<data_t>(axes)));
        fillin(bsinvDiagH,invDiagH,par.bs_controlx,par.bs_controlz);
    }
    if (par.prior_file != "none"){
        bsprior = std::make_shared<vecReg<data_t> >(vecReg<data_t>(hypercube<data_t>(axes)));
        fillin(bsprior,prior,par.bs_controlx,par.bs_controlz);
    }

    duplicate D(*bsmodel->getHyper(),par.bs_mx,par.bs_mz);
    bsplines3 B(*D.getRange(),*model->getHyper(),kx,kz);
    BD = new chainLOper(&B,&D);
}    
// ----------------------------------------------------------------------------------------//

// Build model precon if soft clipping is activated
// ----------------------------------------------------------------------------------------//
    emodelSoftClip * S;
    if (par.soft_clip) {
        S = new emodelSoftClip(*model->getHyper(), par.vpmin, par.vpmax, par.vsmin, par.vsmax, par.rhomin, par.rhomax, 1/sqrt(2.00001), 9, 9);
        if (par.verbose>0) fprintf(stderr,"Soft clipping is added to the inversion. It overrides the hard clipping\n");
    }
// ----------------------------------------------------------------------------------------//

    nloper * op = nullptr;
    nl_we_op * L;
    if (par.nmodels==2) L=new nl_we_op_a(*model->getHyper(),allsrc,par);
    else if (par.nmodels==3) L=new nl_we_op_e(*model->getHyper(),allsrc,par);
    else if (par.nmodels==5) L=new nl_we_op_vti(*model->getHyper(),allsrc,par);

    if (par.bsplines)
    {
        if (par.soft_clip)
        {
            op  = new chainNLOper(S,BD);
        }
        else{
            op  = BD;
        }
    }
    else
    {
        if (par.soft_clip)
        {
            op  = S;
        }
    }

    nloper * D = nullptr;
    loper * R;
    if (par.regularization==0) R = new identity (*model->getHyper());
    else if (par.regularization==1) R = new gradient2d(*model->getHyper());
    else if (par.regularization==2) R = new laplacian2d(*model->getHyper());
    else R = nullptr;
    if (R!=nullptr){
        if (op==nullptr) D = R->clone();
        else D = new chainNLOper(R, op);
        delete R;
    }

    nlls_fwi_eco * prob;
    if (D != nullptr) {
        prob = new nlls_fwi_reg(L, D, bsmodel, data, par.lambda, bsprior, op, bsmask, w, par.normalize, par.integrate, par.envelop);
    }
    else prob = new nlls_fwi_eco(L, bsmodel, data, op, bsmask, w, par.normalize, par.integrate, par.envelop);

    lsearch * ls;
    if (par.lsearch=="weak_wolfe") ls = new weak_wolfe(par.ls_c1, par.ls_a0, par.ls_a1, par.ls_version);
    else if(par.lsearch=="strong_wolfe") ls = new strong_wolfe(par.ls_c1, par.ls_c2, par.ls_a0, par.ls_a1, par.ls_max_step, par.ls_version);
    else ls = new regular_wolfe(par.ls_c1, par.ls_c2, par.ls_a0, par.ls_a1, par.ls_max_step, par.ls_version);

    nlsolver * solver;
    if (par.nlsolver=="nlsd") solver = new nlsd(par.niter, par.max_trial, par.threshold, ls); 
    else if (par.nlsolver=="nlcg") solver = new nlcg(par.niter, par.max_trial, par.threshold, ls); 
    else if (par.nlsolver=="bfgs") solver = new bfgs(par.niter, par.max_trial, par.threshold, ls); 
    else solver = new lbfgs(par.niter, par.max_trial, par.threshold, ls, bsinvDiagH,5); 
    
    solver->run(prob, par.verbose>0, ioutput_file, par.isave);
    
    if (D != nullptr) delete D;
    if (op != nullptr)
    {
        op->forward(false, bsmodel, model);
        delete op;
    }

    if (output_file!="none") sepWrite<data_t>(model, output_file);
    if (obj_func_file!="none") {
        std::shared_ptr<vecReg<data_t> > func = std::make_shared<vecReg<data_t> > (hypercube<data_t>(solver->_func.size()));
        memcpy(func->getVals(), solver->_func.data(), solver->_func.size()*sizeof(data_t));
        sepWrite(func,obj_func_file);
        if (D != nullptr){
            std::shared_ptr<vecReg<data_t> > dfunc = std::make_shared<vecReg<data_t> > (hypercube<data_t>(prob->_dfunc.size()));
            std::shared_ptr<vecReg<data_t> > mfunc = std::make_shared<vecReg<data_t> > (hypercube<data_t>(prob->_mfunc.size()));
            memcpy(dfunc->getVals(), prob->_dfunc.data(), prob->_dfunc.size()*sizeof(data_t));
            memcpy(mfunc->getVals(), prob->_mfunc.data(), prob->_mfunc.size()*sizeof(data_t));
            sepWrite(dfunc,obj_func_file+".d");
            sepWrite(mfunc,obj_func_file+".m");
        }
    }

    delete L;
    delete prob;

return 0;
}