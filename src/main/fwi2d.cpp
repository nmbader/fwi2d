#ifdef DOUBLE_PRECISION
    typedef double data_t;
#else
    typedef float data_t;
#endif

#include <string.h>
#include "we_op.hpp"
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
    std::string source_file="none", model_file="none", data_file="none", output_file="none", ioutput_file="none", obj_func_file="none";
    std::string mask_file="none", weights_file="none";
    readParam<std::string>(argc, argv, "source", source_file);
    readParam<std::string>(argc, argv, "model", model_file);
    readParam<std::string>(argc, argv, "data", data_file);
    readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "ioutput", ioutput_file);
    readParam<std::string>(argc, argv, "obj_func", obj_func_file);
    readParam<std::string>(argc, argv, "mask", mask_file);
    readParam<std::string>(argc, argv, "weights", weights_file);

    successCheck(source_file!="none",__FILE__,__LINE__,"Source wavelet is not provided\n");
    successCheck(model_file!="none",__FILE__,__LINE__,"Earth model is not provided\n");
    successCheck(data_file!="none",__FILE__,__LINE__,"Data to be inverted is not provided\n");

    std::shared_ptr<vec> src = sepRead<data_t>(source_file);
    std::shared_ptr<vec> data = sepRead<data_t>(data_file);
    std::shared_ptr<vec> model = sepRead<data_t>(model_file);

    std::shared_ptr<vec> gmask = nullptr;
    std::shared_ptr<vec> w = nullptr;
    if (mask_file!="none") gmask = sepRead<data_t>(mask_file);
    if (weights_file!="none") w = sepRead<data_t>(weights_file);
    
    // Analyze the input source time function and duplicate if necessary
    nl_we_op_e op(*model->getHyper(),src,par);
    nlls_fwi prob(&op, model, data, gmask, w);

    analyzeNLInversion(par);

    lsearch * ls;
    if (par.lsearch=="weak_wolfe") ls = new weak_wolfe();
    else if(par.lsearch=="strong_wolfe") ls = new strong_wolfe();
    else ls = new regular_wolfe();

    nlsolver * solver;
    if (par.nlsolver=="nlsd") solver = new nlsd(par.niter, par.max_trial, par.threshold, ls); 
    else if (par.nlsolver=="nlcg") solver = new nlcg(par.niter, par.max_trial, par.threshold, ls); 
    else if (par.nlsolver=="bfgs") solver = new bfgs(par.niter, par.max_trial, par.threshold, ls); 
    else solver = new lbfgs(par.niter, par.max_trial, par.threshold, ls); 
    
    solver->run(&prob, par.solver_verbose, ioutput_file, par.isave);

    if (output_file!="none") sepWrite<data_t>(model, output_file);
    if (obj_func_file!="none") {
        std::shared_ptr<vecReg<data_t> > func = std::make_shared<vecReg<data_t> > (hypercube<data_t>(solver->_func.size()));
        memcpy(func->getVals(), solver->_func.data(), solver->_func.size()*sizeof(data_t));
        sepWrite(func,obj_func_file);
    }

return 0;
}