#Required CMake version
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

#Name of our cmake project
project(fwi2d LANGUAGES CXX)


#Require c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -ffast-math -fopenmp -w")
#enable_language(ISPC)

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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g") # add the flap -g for debugging using gdb
    message("Debug mode activated")
endif()

if(GPU)
    enable_language(CUDA)
    message("CUDA enabled")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -w -O3 -ffast-math -std=c++11")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-pie")
endif()

#Required openMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#Include default directory
include_directories(/usr/local/include)
include_directories(/usr/include)

#Include SEPLIB directory
include_directories(/opt/SEP/SEP7.0/include)

#Include CUDA directory
include_directories(/usr/local/cuda/include)

#Adding path to link default libraries
link_directories(/usr/local/lib)

#Adding path to link SEPLIB library
link_directories(/opt/SEP/SEP7.0/lib)

#Adding path to link tcmalloc library
link_directories(/opt/gperftools/local/lib)

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

#add the rules for building the C++ code (execute the CMakelists file in the subdirectory)
add_subdirectory(src)
