#Add directories containing the header files for the code
include_directories(include ../external/nlsolvers)
set(INC_DIR include)

#set the directory to source code
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp)
set(MAIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/main)

set(SRC_FWI2D ${SRC_DIR}/injector.cpp
                ${SRC_DIR}/lsolver.cpp
                ${SRC_DIR}/misc.cpp
                ${SRC_DIR}/nlsolver.cpp
                ${SRC_DIR}/operator.cpp
                ${SRC_DIR}/param.cpp
                ${SRC_DIR}/spatial_operators.cpp
                ${SRC_DIR}/we_op.cpp
	)

set(SRC_NLSOLVER ${CMAKE_CURRENT_SOURCE_DIR}/../external/nlsolvers/lbfgs_c.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../external/nlsolvers/nlcg.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../external/nlsolvers/steepest.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../external/nlsolvers/lsearch.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/../external/nlsolvers/la.cpp
)

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

if(DOUBLE_PRECISION)
        add_custom_command(OUTPUT spatial_operators_ISPC.o
                        COMMAND ispc ${CMAKE_CURRENT_SOURCE_DIR}/cpp/spatial_operators.ispc -DDOUBLE_PRECISION -I ${INC_DIR} -o spatial_operators_ISPC.o --target=avx2-i32x8 --wno-perf
                        DEPENDS ${SRC_DIR}/spatial_operators.ispc
                        )
else()
        add_custom_command(OUTPUT spatial_operators_ISPC.o
                        COMMAND ispc ${CMAKE_CURRENT_SOURCE_DIR}/cpp/spatial_operators.ispc -I ${INC_DIR} -o spatial_operators_ISPC.o --target=avx2-i32x8 --wno-perf
                        DEPENDS ${SRC_DIR}/spatial_operators.ispc
                        )
endif()

#building the project library
add_library(fwi2d spatial_operators_ISPC.o ${SRC_FWI2D})
target_link_libraries(fwi2d fftw3f fftw3)

#add_library(nlsolvers ${SRC_NLSOLVER})

#building the executables
#add_executable(TEST.x ${MAIN_DIR}/test.cpp)
#target_link_libraries(TEST.x fwi2d sep)

#add_executable(TEST_FD.x ${MAIN_DIR}/test_fd.cpp)
#target_link_libraries(TEST_FD.x fwi2d sep)

add_executable(WE_MODELING.x ${MAIN_DIR}/we_modeling.cpp)
target_link_libraries(WE_MODELING.x fwi2d sep)

add_executable(FWI2D.x ${MAIN_DIR}/fwi2d.cpp)
target_link_libraries(FWI2D.x fwi2d sep)

#add_executable(FWI2D0.x ${MAIN_DIR}/fwi2d0.cpp)
#target_link_libraries(FWI2D0.x fwi2d nlsolvers sep)

add_executable(BSPLINES.x ${MAIN_DIR}/bsplines.cpp)
target_link_libraries(BSPLINES.x fwi2d sep)

#installing libraries and executables
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS fwi2d DESTINATION lib)
install(TARGETS WE_MODELING.x
                FWI2D.x
                #FWI2D0.x
                #BSPLINES.x
                #TEST_FD.x
                #TEST.x
        DESTINATION bin)