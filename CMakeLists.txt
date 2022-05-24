#Required CMake version
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

#Name of our cmake project
project(fwi2d LANGUAGES CXX)


#Require c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -ffast-math -fopenmp -w -no-pie")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -fopenmp -w -fPIC")

#optional openMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#optional MPI
if(MPI)
    message("MPI activated")
    find_package(MPI)
    if (MPI_FOUND)
        message("MPI enabled")
        include_directories( ${MPI_CXX_INCLUDE_PATH} )
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_MPI")
    endif(MPI_FOUND)
endif(MPI)

if(GPU)
    enable_language(CUDA)
    message("CUDA enabled")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++11 -O3 -w -Xcompiler -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCUDA -fPIC")
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

if(DOUBLE_PRECISION)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDOUBLE_PRECISION")
    message("Double precision activated")
else()
    message("Single precision activated")
endif()

if(PROFILING)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg") # add the flap -pg for profiling using gprof
    message("Code profiling activated")
endif()

if(DEBUG)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -DDEBUG") # add the flap -g for debugging using gdb
    message("Debug mode activated")
endif()

#Include SEPLIB directory
include_directories(./external/SEP/local/include)

#Adding path to link SEPLIB library
link_directories(./external/SEP/local/lib)

#add the rules for building the C++ code (execute the CMakelists file in the subdirectory)
add_subdirectory(src)
