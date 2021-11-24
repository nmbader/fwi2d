#ifdef DOUBLE_PRECISION
    typedef double data_t;
#else
    typedef float data_t;
#endif

#include <string.h>
#include "we_op.hpp"
#include "IO.hpp"
#include "seplib.h"

typedef vecReg<data_t> vec;
typedef cvecReg<data_t> cvec;
typedef axis<data_t> ax;
typedef hypercube<data_t> hyper;

// Executable to model 2D seismic data and optionally save the full wavefield

int main(int argc, char **argv){
    
    initpar(argc,argv);

    // Read parameters for wave propagation
    param par;
    readParameters(argc, argv, par);

    // Read inputs/outputs files
    std::string source_file="none", model_file="none", wavefield_file="none", output_file="none";
    readParam<std::string>(argc, argv, "source", source_file);
    readParam<std::string>(argc, argv, "model", model_file);
    readParam<std::string>(argc, argv, "wavefield", wavefield_file);
    readParam<std::string>(argc, argv, "output", output_file);

    successCheck(source_file!="none",__FILE__,__LINE__,"Source wavelet is not provided\n");
    successCheck(model_file!="none",__FILE__,__LINE__,"Earth model is not provided\n");

    std::shared_ptr<vec> src = sepRead<data_t>(source_file);
    std::shared_ptr<vec> model = sepRead<data_t>(model_file);
    
    // Analyze the input source time function and duplicate if necessary, analyze geometry
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0);
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);

    // If more than one shot is modeled, don't save the wavefield
    // if (par.ns>1) par.sub=0;

    // Build the appropriate wave equation operator
    nl_we_op_e * op;
    if (par.nmodels==3) op=new nl_we_op_e(*model->getHyper(),allsrc,par);
    else if (par.nmodels==5) op=new nl_we_op_vti(*model->getHyper(),allsrc,par);

    // Run the forward modeling
    std::shared_ptr<vec> allrcv = std::make_shared<vec> (*op->getRange());
    op->forward(false,model,allrcv);

    if ((wavefield_file!="none") && (op->_par.sub)>0) sepWrite<data_t>(op->_full_wfld, wavefield_file);
    if (output_file!="none") sepWrite<data_t>(allrcv, output_file);

    delete op;

return 0;
}