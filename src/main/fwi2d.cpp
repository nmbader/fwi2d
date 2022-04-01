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

    successCheck(source_file!="none",__FILE__,__LINE__,"Source wavelet is not provided\n");
    successCheck(model_file!="none",__FILE__,__LINE__,"Earth model is not provided\n");
    successCheck(data_file!="none",__FILE__,__LINE__,"Data to be inverted is not provided\n");

    std::shared_ptr<vec> src = read<data_t>(source_file, par.format);
    std::shared_ptr<vec> data = read<data_t>(data_file, par.format);
    std::shared_ptr<vec> model = read<data_t>(model_file, par.format);

    std::shared_ptr<vec> hrz = nullptr;
    std::shared_ptr<vec> gmask = nullptr;
    std::shared_ptr<vec> w = nullptr;
    std::shared_ptr<vec> invDiagH = nullptr;
    std::shared_ptr<vec> prior = nullptr;
    if (par.mask_file!="none") {gmask = read<data_t>(par.mask_file, par.format); successCheck(gmask->getN123()==model->getN123(),__FILE__,__LINE__,"Gradient mask must have the same number of samples as the model\n");}
    if (par.weights_file!="none") {w = read<data_t>(par.weights_file, par.format); successCheck(gmask->getN123()==model->getN123(),__FILE__,__LINE__,"Data weights must have the same number of samples as the data\n");}
    if (par.inverse_diagonal_hessian_file!="none") {invDiagH = read<data_t>(par.inverse_diagonal_hessian_file, par.format); successCheck(gmask->getN123()==model->getN123(),__FILE__,__LINE__,"Inverse diagonal Hessian must have the same number of samples as the model\n");}
    if (par.prior_file!="none") {prior = read<data_t>(par.prior_file, par.format); successCheck(gmask->getN123()==model->getN123(),__FILE__,__LINE__,"Prior model must have the same number of samples as the model\n");}

// Analyze the inputs and parameters and modify if necessary
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0);
    analyzeBsplines(*model->getHyper(),par);
    analyzeNLInversion(par);

// Build model precon if 1D inversion is activated
// ----------------------------------------------------------------------------------------//
    std::shared_ptr<vec> model1d = model;
    std::shared_ptr<vec> gmask1d = gmask;
    std::shared_ptr<vec> invDiagH1d = invDiagH;
    std::shared_ptr<vec> prior1d = prior;
    loper * ext1d;

if (par.inversion1d)
{
    if (par.horizon_file!="none") {
        hrz = read<data_t>(par.horizon_file, par.format);
        successCheck(hrz->getN123()==model->getHyper()->getAxis(2).n,__FILE__,__LINE__,"The provided horizon must have the same size as the x-dimension of the model\n");
    }
    else{
        hrz = std::make_shared<vec>(hyper(model->getHyper()->getAxis(2)));
        hrz->zero();
    }
    std::vector<ax> axes = model->getHyper()->getAxes();
    axes.erase(axes.begin()+1);
    model1d = std::make_shared<vec> (hyper(axes));
    model1d->zero();
    ext1d = new extrapolator1d2d(*model1d->getHyper(),hrz);
    ext1d->inverse(false,model1d,model);

    if (par.mask_file != "none") {
        gmask1d = std::make_shared<vec> (hyper(axes));
        gmask1d->zero();
        ext1d->inverse(false,gmask1d,gmask);
    }
    if (par.inverse_diagonal_hessian_file != "none"){
        invDiagH1d = std::make_shared<vec> (hyper(axes));
        invDiagH1d->zero();
        ext1d->inverse(false,invDiagH1d,invDiagH);
    }
    if (par.prior_file != "none"){
        prior1d = std::make_shared<vec> (hyper(axes));
        prior1d->zero();
        ext1d->inverse(false,prior1d,prior);
    }
}
// ----------------------------------------------------------------------------------------//

// Build model precon if B-splines are activated
// ----------------------------------------------------------------------------------------//
    std::shared_ptr<vec> bsmodel = model1d;
    std::shared_ptr<vec> bsmask = gmask1d;
    std::shared_ptr<vec> bsinvDiagH = invDiagH1d;
    std::shared_ptr<vec> bsprior = prior1d;
    loper * BD;

if (par.bsplines)
{
    std::vector<data_t> kx;
    std::vector<data_t> kz;
    setKnot(kx,par.bs_controlx,par.bs_mx);
    setKnot(kz,par.bs_controlz,par.bs_mz);

    std::vector<ax > axes = model->getHyper()->getAxes();
    axes[0].n=par.bs_nz; axes[1].n=par.bs_nx;
    if (par.inversion1d) axes.erase(axes.begin()+1);
    bsmodel = std::make_shared<vec>(vec(hyper(axes)));
    if (par.inversion1d) fillin1d(bsmodel,model1d,par.bs_controlz);
    else fillin(bsmodel,model,par.bs_controlx,par.bs_controlz);
    if (par.mask_file != "none"){
        bsmask = std::make_shared<vec>(vec(hyper(axes)));
        if (par.inversion1d) fillin1d(bsmask,gmask1d,par.bs_controlz);
        else fillin(bsmask,gmask,par.bs_controlx,par.bs_controlz);
    }
    if (par.inverse_diagonal_hessian_file != "none"){
        bsinvDiagH = std::make_shared<vec>(vec(hyper(axes)));
        if (par.inversion1d) fillin1d(bsinvDiagH,invDiagH1d,par.bs_controlz);
        else fillin(bsinvDiagH,invDiagH,par.bs_controlx,par.bs_controlz);
    }
    if (par.prior_file != "none"){
        bsprior = std::make_shared<vec>(vec(hyper(axes)));
        if (par.inversion1d) fillin1d(bsprior,prior1d,par.bs_controlz);
        else fillin(bsprior,prior,par.bs_controlx,par.bs_controlz);
    }

    if (par.inversion1d){
        duplicate1d D(*bsmodel->getHyper(),par.bs_mz);
        bsplines31d B(*D.getRange(),*model1d->getHyper(),kz);
        BD = new chainLOper(&B,&D);
    }
    else{
        duplicate D(*bsmodel->getHyper(),par.bs_mx,par.bs_mz);
        bsplines3 B(*D.getRange(),*model->getHyper(),kx,kz);
        BD = new chainLOper(&B,&D);
    }
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
        if (par.soft_clip)
        {
            if (par.inversion1d) {
                chainLOper EBD(ext1d,BD);
                op = new chainNLOper(S,&EBD);
            }
            else op  = new chainNLOper(S,BD);
        }
        else{
            if (par.inversion1d) op = new chainLOper(ext1d,BD);
            else op  = BD->clone();
        }
    }
    else
    {
        if (par.soft_clip)
        {
            if (par.inversion1d) op = new chainNLOper(S,ext1d);
            else op  = S->clone();
        }
        else if (par.inversion1d) op = ext1d->clone();
    }
    if (par.bsplines) delete BD;
    if (par.inversion1d) delete ext1d;
    if (par.soft_clip) delete S;

    nloper * D = nullptr;
    loper * R;
    if (par.regularization==0) R = new identity (*model->getHyper());
    else if (par.regularization==1) R = new gradient2d(*model->getHyper(),par.reg_xweight,par.reg_zweight);
    else if (par.regularization==2) R = new laplacian2d(*model->getHyper(),par.reg_xweight,par.reg_zweight);
    else R = nullptr;
    if (R!=nullptr){
        if (op==nullptr) D = R->clone();
        else D = new chainNLOper(R, op);
        delete R;
    }

    nlls_fwi_eco * prob;
    if (D != nullptr) {
        prob = new nlls_fwi_reg(L, D, bsmodel, data, par.lambda, bsprior, op, bsmask, w);
    }
    else prob = new nlls_fwi_eco(L, bsmodel, data, op, bsmask, w);

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
    
    if (D != nullptr) delete D;
    if (op != nullptr)
    {
        op->forward(false, bsmodel, model);
        delete op;
    }

    if (rank==0 && output_file!="none") write<data_t>(model, output_file, par.format, par.datapath);
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

    delete L;
    delete prob;

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif

return 0;
}