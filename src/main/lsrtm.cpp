#include <string.h>
#include "we_op.hpp"
#include "IO.hpp"
#include "optimization.hpp"
#include "lsolver.hpp"
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
    std::string source_file="none", model_file="none", data_file="none", output_file="none", solver="cgls", obj_func_file="none";
    int niter=10;
    data_t threshold=0;
    readParam<std::string>(argc, argv, "source", source_file);
    readParam<std::string>(argc, argv, "model", model_file);
    readParam<std::string>(argc, argv, "data", data_file);
    readParam<std::string>(argc, argv, "output", output_file);
    readParam<std::string>(argc, argv, "solver", solver);
    readParam<std::string>(argc, argv, "obj_func", obj_func_file);
    readParam<int>(argc, argv, "niter", niter);
    readParam<data_t>(argc, argv, "threshold", threshold);

    successCheck(source_file!="none",__FILE__,__LINE__,"Source wavelet is not provided\n");
    successCheck(model_file!="none",__FILE__,__LINE__,"Earth model is not provided\n");
    successCheck(data_file!="none",__FILE__,__LINE__,"Data is not provided\n");

    std::shared_ptr<vec> src = read<data_t>(source_file, par.format);
    std::shared_ptr<vec> model = read<data_t>(model_file, par.format);
    std::shared_ptr<vec> data = read<data_t>(data_file, par.format);

// Initialize the reflectivity
    int nx = model->getHyper()->getAxis(2).n;
    int nz = model->getHyper()->getAxis(1).n;
    std::shared_ptr<vec> refl = std::make_shared<vec>(hyper(model->getHyper()->getAxis(1),model->getHyper()->getAxis(2)));
    refl->zero();

// Analyze the input source time function and duplicate if necessary, analyze geometry, analyze model
    analyzeGeometry(*model->getHyper(),par, par.verbose>0);
    std::shared_ptr<vec> allsrc = analyzeWavelet(src, par, par.verbose>0);
    analyzeModel(*allsrc->getHyper(),model,par);

// Build the Born operator
    born_op_a op(model,allsrc,par);
    successCheck(data->getN123()==op.getRange()->getN123(),__FILE__,__LINE__,"The provided data is incompatible\n");

// set the least-squares problem
    llsq prob(&op, refl, data);

// Set the linear solver
    lsolver * sol;
    if (solver == "sdls") sol = new sdls(niter,threshold);
    else sol = new cgls(niter,threshold);

// Invert the data
    sol->run(&prob, par.verbose>0);

    if ((rank==0) && (output_file!="none")) write<data_t>(refl, output_file, par.format, par.datapath);
    if ((rank==0) && (obj_func_file!="none")) {
        std::shared_ptr<vec > func = std::make_shared<vec > (hyper(sol->_func.size()));
        memcpy(func->getVals(), sol->_func.data(), sol->_func.size()*sizeof(data_t));
        write(func,obj_func_file, par.format, par.datapath);
    }

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif

return 0;
}