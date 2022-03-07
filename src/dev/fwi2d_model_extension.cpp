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

// structure to contain sources location, control points and knot vectors for B-splines model extension
struct bs_sources{
    std::vector<std::vector<data_t> > controlx; // control points
    std::vector<std::vector<data_t> > controlz;
    std::vector<std::vector<int> > mx; // multiplicity of control points
    std::vector<std::vector<int> > mz;
    std::vector<std::vector<data_t> > kx; // knots
    std::vector<std::vector<data_t> > kz;
};

void set_bs_sources(bs_sources &bss, const hyper &domain, std::vector<std::vector<data_t> > &sxz, int gx, int gz){
    
    int ns=sxz.size();
    ax Z=domain.getAxis(1);
    ax X=domain.getAxis(2);

    bss.controlx.resize(ns);
    bss.controlz.resize(ns);
    bss.mx.resize(ns);
    bss.mz.resize(ns);
    bss.kx.resize(ns);
    bss.kz.resize(ns);

    // make the gap always odd
    int gapx = floor(gx/2.0);
    int gapz = floor(gz/2.0);
    gapx = 2*gapx + 1;
    gapz = 2*gapz + 1;
    if (gx==0) gapx=0;
    if (gz==0) gapz=0;

    for (int s=0; s<ns; s++){

        bss.controlx[s].resize(X.n-gapx);
        bss.controlz[s].resize(Z.n-gapz);
        bss.mx[s].resize(X.n-gapx);
        bss.mz[s].resize(Z.n-gapz);

        int ixs = round((sxz[s][0] - X.o)/X.d);
        int izs = round((sxz[s][1] - Z.o)/Z.d);

        for (int i=0; i<ixs-floor(gapx/2.0); i++){
            bss.controlx[s][i] = X.o + i*X.d;
            bss.mx[s][i]=1;
        }
        for (int i=ixs+ceil(gapx/2.0); i<X.n; i++){
            bss.controlx[s][i-gapx] = X.o + i*X.d;
            bss.mx[s][i-gapx]=1;
        }
        bss.mx[s][0]=2;
        bss.mx[s][X.n-1-gapx]=2;
        setKnot(bss.kx[s],bss.controlx[s],bss.mx[s]);

        for (int i=0; i<izs-floor(gapz/2.0); i++){
            bss.controlz[s][i] = Z.o + i*Z.d;
            bss.mz[s][i]=1;
        }
        for (int i=izs+ceil(gapz/2.0); i<Z.n; i++){
            bss.controlz[s][i-gapz] = Z.o + i*Z.d;
            bss.mz[s][i-gapz]=1;
        }
        bss.mz[s][0]=2;
        bss.mz[s][Z.n-1-gapz]=2;
        setKnot(bss.kz[s],bss.controlz[s],bss.mz[s]);
    }
}

// operator that maps from several B-spline models (one for each source) to the physical space, with optional damping (cosine squared) around the sources location
// Op = sum(Ws.Ss.ps) where summation is over sources, 'ps' is the B-spline model for source 's', 'Ss' is the B-spline back-projection operator, 'Ws' is the normalized cosine squared weighting operator
class model_extension : public loper {
protected:
    data_t _dpower; // damping strength for cosine^2
    data_t _dwidthx; // damping width
    data_t _dwidthz;
    //std::vector<loper *> _BDs; // all B-splines operators (one per source)
    std::shared_ptr<vec> _nW; // normalization for the weighting (sum of all weighting operators)
    std::vector<std::vector<data_t> > _sxz; // sources location

public:
    std::vector<loper *> _BDs; // all B-splines operators (one per source)
    model_extension(){}
    ~model_extension(){
        int ns = _sxz.size();
        for (int s=0; s<ns; s++) delete _BDs[s];
    }
    model_extension(const hypercube<data_t> &range, bs_sources &bss, std::vector<std::vector<data_t> > &sxz, data_t dwidthx, data_t dwidthz, data_t dpower){
        successCheck(range.getNdim() == 3,__FILE__,__LINE__,"Range must contain 3 dimensions\n");
        successCheck(sxz.size() == bss.mx.size(),__FILE__,__LINE__,"Number of sources must be the same in 'sxz' and 'bss'\n");

        std::vector<ax > axes = range.getAxes();
        axes[0].n=bss.mz[0].size(); axes[1].n=bss.mx[0].size();
        ax S(sxz.size(),0,1);
        axes.push_back(S);
        
        _domain = hyper(axes);
        _range = range;
        _dwidthx = dwidthx;
        _dwidthz = dwidthz;
        _dpower = dpower;
        _sxz=sxz;

        // build the B-spline operators and the normalization weights
        ax Z = range.getAxis(1);
        ax X = range.getAxis(2);
        ax C = range.getAxis(3);
        int ns = sxz.size();
        _nW = std::make_shared<vec> (vec(hyper(Z,X)));
        _nW->zero();
        data_t (*pw) [Z.n] = (data_t (*)[Z.n]) _nW->getVals();

        for (int s=0; s<ns; s++){

            // B-spline operator
            hyper hyp(ax(bss.mz[s].size(),0,1),ax(bss.mx[s].size(),0,1),C);
            duplicate D(hyp,bss.mx[s],bss.mz[s]);
            bsplines3 B(*D.getRange(),range,bss.kx[s],bss.kz[s]);
            loper * BD = new chainLOper(&B,&D);
            _BDs.push_back(BD);

            // normalization weight
            int ixmin = ceil((sxz[s][0]-dwidthx-X.o)/X.d);
            int ixmax = ceil((sxz[s][0]+dwidthx-X.o)/X.d);
            int izmin = ceil((sxz[s][1]-dwidthz-Z.o)/Z.d);
            int izmax = ceil((sxz[s][1]+dwidthz-Z.o)/Z.d);

            for (int ix=ixmin; ix<ixmax; ix++){
                data_t x = X.o + ix*X.d;
                data_t rx = (x-_sxz[s][0])/_dwidthx;
                rx*=rx;
                for (int iz=izmin; iz<izmax; iz++){
                    data_t z = Z.o + iz*Z.d;
                    data_t rz = (z-_sxz[s][1])/_dwidthz;
                    rz*=rz;
                    data_t d = sqrt(rx+rz);
                    data_t val = 1;
                    if (d<1) val = cos(_dpower*0.5*M_PI*(1-d));
                    pw[ix][iz] += val*val;
                }
            }

            for (int ix=0; ix<ixmin; ix++){
                for (int iz=0; iz<Z.n; iz++) pw[ix][iz] += 1;
            }
            for (int ix=ixmax; ix<X.n; ix++){
                for (int iz=0; iz<Z.n; iz++) pw[ix][iz] += 1;
            }
            for (int ix=ixmin; ix<ixmax; ix++){
                for (int iz=0; iz<izmin; iz++) pw[ix][iz] += 1;
                for (int iz=izmax; iz<Z.n; iz++) pw[ix][iz] += 1;
            }
        }
    }

    model_extension * clone() const {
        model_extension * op = new model_extension();
        op->_dpower = _dpower;
        op->_dwidthx = _dwidthx;
        op->_dwidthz = _dwidthz;
        op->_domain = _domain;
        op->_range = _range;
        op->_sxz = _sxz;
        op->_nW = _nW->clone();
        int ns = _sxz.size();
        for (int s=0; s<ns; s++) {
            loper * BD = _BDs[s]->clone();
            op->_BDs.push_back(BD);
        }
        return op;
    }    
    void apply_forward(bool add, const data_t * pmod, data_t * pdat){

        if (add == false) memset(pdat, 0, _range.getN123()*sizeof(data_t));

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        int ns = _sxz.size();
        int nm = _domain.getN123()/ns;
        data_t (*pd) [X.n][Z.n] = (data_t (*)[X.n][Z.n]) pdat;
        const data_t (*pw) [Z.n] = (const data_t (*)[Z.n]) _nW->getVals();

        data_t * vec = new data_t[_range.getN123()];
        memset(vec, 0, _range.getN123()*sizeof(data_t));
        data_t (*pvec) [X.n][Z.n] = (data_t (*)[X.n][Z.n]) vec;

        for (int s=0; s<ns; s++){

            // Ss.ps
            _BDs[s]->apply_forward(false, pmod+s*nm, vec);

            // Ws.Ss.ps
            int ixmin = ceil((_sxz[s][0]-_dwidthx-X.o)/X.d);
            int ixmax = ceil((_sxz[s][0]+_dwidthx-X.o)/X.d);
            int izmin = ceil((_sxz[s][1]-_dwidthz-Z.o)/Z.d);
            int izmax = ceil((_sxz[s][1]+_dwidthz-Z.o)/Z.d);
            for (int c=0; c<C.n; c++){

                #pragma omp parallel for
                for (int ix=ixmin; ix<ixmax; ix++){
                    data_t x = X.o + ix*X.d;
                    data_t rx = (x-_sxz[s][0])/_dwidthx;
                    rx*=rx;
                    for (int iz=izmin; iz<izmax; iz++){
                        data_t z = Z.o + iz*Z.d;
                        data_t rz = (z-_sxz[s][1])/_dwidthz;
                        rz*=rz;
                        data_t d = sqrt(rx+rz);
                        data_t val = 1;
                        if (d<1) val = cos(_dpower*0.5*M_PI*(1-d));
                        pvec[c][ix][iz] *= val*val;
                    }
                }

                #pragma omp parallel for
                for (int ix=0; ix<X.n; ix++){
                    for (int iz=0; iz<Z.n; iz++) pd[c][ix][iz] += pvec[c][ix][iz]/pw[ix][iz];
                }
            }            
        }
        delete [] vec;
    }

    void apply_adjoint(bool add, data_t * pmod, const data_t * pdat){

        if (add == false) memset(pmod, 0, _domain.getN123()*sizeof(data_t));

        ax Z = _range.getAxis(1);
        ax X = _range.getAxis(2);
        ax C = _range.getAxis(3);
        int ns = _sxz.size();
        int nm = _domain.getN123()/ns;
        const data_t (*pw) [Z.n] = (const data_t (*)[Z.n]) _nW->getVals();

        data_t * vec = new data_t[_range.getN123()];
        memset(vec, 0, _range.getN123()*sizeof(data_t));
        data_t (*pvec) [X.n][Z.n] = (data_t (*)[X.n][Z.n]) vec;

        for (int s=0; s<ns; s++){

            memcpy(vec, pdat, _range.getN123()*sizeof(data_t));

            // Ws'.m
            int ixmin = ceil((_sxz[s][0]-_dwidthx-X.o)/X.d);
            int ixmax = ceil((_sxz[s][0]+_dwidthx-X.o)/X.d);
            int izmin = ceil((_sxz[s][1]-_dwidthz-Z.o)/Z.d);
            int izmax = ceil((_sxz[s][1]+_dwidthz-Z.o)/Z.d);
            for (int c=0; c<C.n; c++){

                #pragma omp parallel for
                for (int ix=ixmin; ix<ixmax; ix++){
                    data_t x = X.o + ix*X.d;
                    data_t rx = (x-_sxz[s][0])/_dwidthx;
                    rx*=rx;
                    for (int iz=izmin; iz<izmax; iz++){
                        data_t z = Z.o + iz*Z.d;
                        data_t rz = (z-_sxz[s][1])/_dwidthz;
                        rz*=rz;
                        data_t d = sqrt(rx+rz);
                        data_t val = 1;
                        if (d<1) val = cos(_dpower*0.5*M_PI*(1-d));
                        pvec[c][ix][iz] *= val*val;
                    }
                }

                #pragma omp parallel for
                for (int ix=0; ix<X.n; ix++){
                    for (int iz=0; iz<Z.n; iz++) pvec[c][ix][iz] /= pw[ix][iz];
                }
            }

            // Ss'.Ws'.m
            _BDs[s]->apply_adjoint(add, pmod+s*nm, vec);
        }
        delete [] vec;
    }
};

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

// Read some parameters
    int gapx=0, gapz=0;
    data_t dwidthx=0, dwidthz=0, dpower=0;
    readParam<int>(argc, argv, "damping_gapx", gapx);
    readParam<int>(argc, argv, "damping_gapz", gapz);
    readParam<data_t>(argc, argv, "damping_widthx", dwidthx);
    readParam<data_t>(argc, argv, "damping_widthz", dwidthz);
    readParam<data_t>(argc, argv, "damping_power", dpower);

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
    if (par.mask_file!="none") gmask = read<data_t>(par.mask_file, par.format);
    if (par.weights_file!="none") w = read<data_t>(par.weights_file, par.format);
    if (par.inverse_diagonal_hessian_file!="none") invDiagH = read<data_t>(par.inverse_diagonal_hessian_file, par.format);
    if (par.prior_file!="none") prior = read<data_t>(par.prior_file, par.format);

// Analyze the inputs and parameters and modify if necessary
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0);
    analyzeBsplines(*model->getHyper(),par);
    analyzeNLInversion(par);

// Build model precon for model extension along sources
// ----------------------------------------------------------------------------------------//
    std::shared_ptr<vec> bsmodel = model;
    std::shared_ptr<vec> bsmask = gmask;
    std::shared_ptr<vec> bsinvDiagH = invDiagH;
    std::shared_ptr<vec> bsprior = prior;
    bs_sources bss;
    set_bs_sources(bss, *model->getHyper(), par.sxz, gapx, gapz);
    loper * E = new model_extension(*model->getHyper(), bss, par.sxz, dwidthx, dwidthz, dpower);

    bsmodel = std::make_shared<vec>(vec(*E->getDomain()));
    int n123 = bsmodel->getN123()/par.sxz.size();  
    if (par.mask_file != "none") bsmask = std::make_shared<vec>(vec(*E->getDomain()));
    if (par.inverse_diagonal_hessian_file != "none") bsinvDiagH = std::make_shared<vec>(vec(*E->getDomain()));
    if (par.prior_file != "none") bsprior = std::make_shared<vec>(vec(*E->getDomain()));

    for (int s=0; s<par.sxz.size(); s++){
        bsfillin F(*model->getHyper(),bss.controlx[s],bss.controlz[s]);
        F.apply_forward(false,model->getCVals(),bsmodel->getVals()+s*n123);
        if (par.mask_file != "none") F.apply_forward(false,gmask->getCVals(),bsmask->getVals()+s*n123);
        if (par.inverse_diagonal_hessian_file != "none") F.apply_forward(false,invDiagH->getCVals(),bsinvDiagH->getVals()+s*n123);
        if (par.prior_file != "none") F.apply_forward(false,prior->getCVals(),bsprior->getVals()+s*n123);
    }

// ----------------------------------------------------------------------------------------//

// Build model precon if B-splines are activated
// ----------------------------------------------------------------------------------------//
    loper * BDF;

if (par.bsplines)
{
    std::vector<data_t> kx;
    std::vector<data_t> kz;
    setKnot(kx,par.bs_controlx,par.bs_mx);
    setKnot(kz,par.bs_controlz,par.bs_mz);

    bsfillin F(*E->getRange(),par.bs_controlx,par.bs_controlz);
    duplicate D(*F.getRange(),par.bs_mx,par.bs_mz);
    bsplines3 B(*D.getRange(),*model->getHyper(),kx,kz);
    chainLOper DF(&D,&F);
    BDF = new chainLOper(&B,&DF);
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
            chainLOper BDFE(BDF,E);
            op = new chainNLOper(S,&BDFE);
        }
        else{
            op = new chainLOper(BDF,E);
        }
    }
    else
    {
        if (par.soft_clip)
        {
            op = new chainNLOper(S,E);
        }
        else op = E->clone();
    }
    if (par.bsplines) delete BDF;
    delete E;
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