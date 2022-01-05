cmake_minimum_required(VERSION 2.8)
project(SEPLIB-base-lib-corelibs-superset)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../include)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../sep)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../sep3d)


set (SOURCES 
axes.c	   drn.c	   header.c   list.c	    rite.c	   tagstream.c
conform.c  extra_keys.c    headval.c  mpi_thread.c  rite_file.c    transfer.c
convert.c  fortran_junk.c  init.c     ntraces.c     rite.inc	   type_usage.c
coord.c    	struct_sepf.f90    reed.c	   
debug.c    grid.c	   keys.c     struct_sepc.c)





add_library(superset ${SOURCES})
target_link_libraries(superset sep3d sep)
install(TARGETS superset DESTINATION lib)

if(SU_LIBRARIES)
include_directories(${SU_INCLUDES})
set(SOURCE2 reed_su.c rite_su.c init_su.c)
add_library(sepsu ${SOURCE2})
target_link_libraries(sepsu sep3d sep)
install(TARGETS sepsu DESTINATION lib)
endif(SU_LIBRARIES)
