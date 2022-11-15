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

    int rank=0, size=0;
#ifdef ENABLE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    fprintf (stderr,"\n====================\nSize of MPI communicator = %d ; current rank = %d\n====================\n",size,rank);
#endif

    initpar(argc,argv);

// Read parameters for wave propagation
    param par;
    readParameters(argc, argv, par);
    int verbose=par.verbose;
    if (rank>0) par.verbose=0;
    par.device+=rank;

// Set the maximum number of threads
    if (par.nthreads>0) omp_set_num_threads(par.nthreads);

// Read inputs/outputs files
    std::string source_file="none", model_file="none", data_file="none", wavefield_file="none", output_file="none";
    readParam<std::string>(argc, argv, "source", source_file);
    readParam<std::string>(argc, argv, "model", model_file);
    readParam<std::string>(argc, argv, "data", data_file);
    readParam<std::string>(argc, argv, "wavefield", wavefield_file);
    readParam<std::string>(argc, argv, "output", output_file);

    successCheck(source_file!="none",__FILE__,__LINE__,"Source wavelet is not provided\n");
    successCheck(model_file!="none",__FILE__,__LINE__,"Earth model is not provided\n");

    std::shared_ptr<vec> src = read<data_t>(source_file, par.format);
    std::shared_ptr<vec> model = read<data_t>(model_file, par.format);
    std::shared_ptr<vec> data = nullptr;
    if (data_file != "none") data = read<data_t>(data_file, par.format);

// Extract the reflectivity
    int nx = model->getHyper()->getAxis(2).n;
    int nz = model->getHyper()->getAxis(1).n;
    std::shared_ptr<vec> refl = std::make_shared<vec>(hyper(model->getHyper()->getAxis(1),model->getHyper()->getAxis(2)));
    memcpy(refl->getVals(),model->getCVals()+nx*nz, nx*nz*sizeof(data_t));
    
// Analyze the input source time function and duplicate if necessary, analyze geometry, analyze model
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0);
    analyzeModel(*allsrc->getHyper(),model,par);

// Build the Born operator and the non-linear forward operator
    born_op_a op(model,allsrc,par);
    nl_we_op_a nlop(*model->getHyper(),allsrc,par);
    if (rank>0) par.verbose=verbose;
    std::shared_ptr<vec> allrcv = std::make_shared<vec> (*op.getRange());
    allrcv->zero();
    std::shared_ptr<vec> output = nullptr;

// Run the BORN forward modeling if applicable
    if (data == nullptr)
    {
        op.forward(false,refl,allrcv);
        output = allrcv;
    }

// Run BORN adjoint modeling (migration) if applicable
    else
    {
        successCheck(data->getN123()==allrcv->getN123(),__FILE__,__LINE__,"The provided data is incompatible\n");

        model->set(1,nx*nz);
        refl->zero();
        nlop.forward(false,model,allrcv);
        op._background_wfld = nlop._full_wfld;

        data->scale(-1); // data needs to be scaled by -1
        op.adjoint(false,refl,data);
        output = refl;
    }

    if ((rank==0) && (wavefield_file!="none") && (op._par.sub)>0) write<data_t>(op._background_wfld, wavefield_file, par.format, par.datapath);
    if ((rank==0) && (output_file!="none")) write<data_t>(output, output_file, par.format, par.datapath);

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif

return 0;
}