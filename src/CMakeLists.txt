#Add directories containing the header files for the code
include_directories(include cuda)
set(INC_DIR include)
set(CUDA_INC_DIR cuda)

#set the directory to source code
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp)
set(CUDA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cuda)
set(PYTHON_DIR ${CMAKE_CURRENT_SOURCE_DIR}/python)
set(MAIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/main)
set(DEV_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dev)

set(SRC_FWI2D ${SRC_DIR}/injector.cpp
                ${SRC_DIR}/lsolver.cpp
                ${SRC_DIR}/misc.cpp
                ${SRC_DIR}/nlsolver.cpp
                ${SRC_DIR}/operator.cpp
                ${SRC_DIR}/param.cpp
                ${SRC_DIR}/spatial_operators.cpp
                ${SRC_DIR}/we_op.cpp
	)

if(ENABLE_CUDA)
        set(SRC_FWI2D ${SRC_FWI2D} ${CUDA_DIR}/cudaKernels.cu
        ${CUDA_DIR}/cudaWeOp.cpp
        )
else()
        set(SRC_FWI2D ${SRC_FWI2D} ${SRC_DIR}/dummy.cpp)
endif()

set(HEADERS ${INC_DIR}/injector.hpp
                ${INC_DIR}/IO.hpp
                ${INC_DIR}/lsolver.hpp
                ${INC_DIR}/misc.hpp
                ${INC_DIR}/nlsolver.hpp
                ${INC_DIR}/operator.hpp
                ${INC_DIR}/optimization.hpp
                ${INC_DIR}/param.hpp
                ${INC_DIR}/spatial_operators.hpp
                ${INC_DIR}/vecReg.hpp
                ${INC_DIR}/we_op.hpp
        )

set(PYTHON ${PYTHON_DIR}/seppy.py
                ${PYTHON_DIR}/sep2numpy.py)

if (ENABLE_CUDA)
        set(HEADERS ${HEADERS} ${CUDA_INC_DIR}/cudaCppWrappers.h
                ${CUDA_INC_DIR}/cudaMisc.h
        )
endif()

if(ENABLE_DOUBLE_PRECISION)
        add_custom_command(OUTPUT spatial_operators_ISPC.o
                        COMMAND ${ISPC_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cpp/spatial_operators.ispc -DDOUBLE_PRECISION -I ${INC_DIR} -o spatial_operators_ISPC.o --target=avx2-i32x8 --wno-perf --pic
                        DEPENDS ${SRC_DIR}/spatial_operators.ispc
                        )
else()
        add_custom_command(OUTPUT spatial_operators_ISPC.o
                        COMMAND ${ISPC_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cpp/spatial_operators.ispc -I ${INC_DIR} -o spatial_operators_ISPC.o --target=avx2-i32x8 --wno-perf --pic
                        DEPENDS ${SRC_DIR}/spatial_operators.ispc
                        )
endif()

#building the project library
add_library(fwi2d spatial_operators_ISPC.o ${SRC_FWI2D})
target_link_libraries(fwi2d fftw3f fftw3)
if (MPI_FOUND)
        target_link_libraries(fwi2d ${MPI_CXX_LIBRARIES} )
endif()

#building the main executables
add_executable(WE_MODELING.x ${MAIN_DIR}/we_modeling.cpp)
target_link_libraries(WE_MODELING.x fwi2d sep)

add_executable(BORN_MODELING.x ${MAIN_DIR}/born_modeling.cpp)
target_link_libraries(BORN_MODELING.x fwi2d sep)

add_executable(FWI2D.x ${MAIN_DIR}/fwi2d.cpp)
target_link_libraries(FWI2D.x fwi2d sep)


add_executable(BSPLINES.x ${MAIN_DIR}/bsplines.cpp)
target_link_libraries(BSPLINES.x fwi2d sep)

add_executable(FX.x ${MAIN_DIR}/fx.cpp)
target_link_libraries(FX.x fwi2d sep)

add_executable(FK.x ${MAIN_DIR}/fk.cpp)
target_link_libraries(FK.x fwi2d sep)

add_executable(GENERATE_WAVELET.x ${MAIN_DIR}/generate_wavelet.cpp)
target_link_libraries(GENERATE_WAVELET.x fwi2d sep)

add_executable(DERIVATIVE.x ${MAIN_DIR}/derivative.cpp)
target_link_libraries(DERIVATIVE.x fwi2d sep)

add_executable(INTEGRAL.x ${MAIN_DIR}/integral.cpp)
target_link_libraries(INTEGRAL.x fwi2d sep)

add_executable(HILBERT.x ${MAIN_DIR}/hilbert.cpp)
target_link_libraries(HILBERT.x fwi2d sep)

add_executable(NOISE.x ${MAIN_DIR}/noise.cpp)
target_link_libraries(NOISE.x fwi2d sep)

add_executable(RESAMPLE.x ${MAIN_DIR}/resample.cpp)
target_link_libraries(RESAMPLE.x fwi2d sep)

add_executable(FK_FILTER.x ${MAIN_DIR}/fk_filter.cpp)
target_link_libraries(FK_FILTER.x fwi2d sep)

add_executable(FX_FILTER.x ${MAIN_DIR}/fx_filter.cpp)
target_link_libraries(FX_FILTER.x fwi2d sep)

add_executable(MATCHING.x ${MAIN_DIR}/matching.cpp)
target_link_libraries(MATCHING.x fwi2d sep)

add_executable(SEP2BIN.x ${MAIN_DIR}/sep2bin.cpp)
target_link_libraries(SEP2BIN.x fwi2d sep)

add_executable(TR2TRDECON.x ${MAIN_DIR}/tr2trDecon.cpp)
target_link_libraries(TR2TRDECON.x fwi2d sep)

#installing libraries and executables
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS fwi2d DESTINATION lib)
install(FILES ${PYTHON} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION bin)
install(TARGETS WE_MODELING.x
                BORN_MODELING.x
                FWI2D.x
                BSPLINES.x
                FX.x
                FK.x
                GENERATE_WAVELET.x
                DERIVATIVE.x
                INTEGRAL.x
                HILBERT.x
                NOISE.x
                RESAMPLE.x
                FK_FILTER.x
                FX_FILTER.x
                MATCHING.x
                SEP2BIN.x
                TR2TRDECON.x
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION bin)

# For development purposes
if (ENABLE_DEV)
        #add_executable(TEST.x ${MAIN_DIR}/test.cpp)
        #target_link_libraries(TEST.x fwi2d sep)

        #add_executable(TEST_FD.x ${MAIN_DIR}/test_fd.cpp)
        #target_link_libraries(TEST_FD.x fwi2d sep)

        add_executable(FWI2D_MODEL_EXTENSION.x ${DEV_DIR}/fwi2d_model_extension.cpp)
        target_link_libraries(FWI2D_MODEL_EXTENSION.x fwi2d sep)

        add_executable(FWI2D_MODEL_EXTENSION2.x ${DEV_DIR}/fwi2d_model_extension2.cpp)
        target_link_libraries(FWI2D_MODEL_EXTENSION2.x fwi2d sep)

        add_executable(LEAKY_TUBE_WAVE.x ${DEV_DIR}/leaky_tw.cpp)
        target_link_libraries(LEAKY_TUBE_WAVE.x fwi2d sep)

        install(TARGETS FWI2D_MODEL_EXTENSION.x
                FWI2D_MODEL_EXTENSION2.x
                LEAKY_TUBE_WAVE.x
                #TEST_FD.x
                #TEST.x
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION bin)
endif()