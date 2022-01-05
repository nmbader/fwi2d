cmake_minimum_required(VERSION 2.8)
project(SEPLIB-base-lib-corelibs-sep3d)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../sep)



set (SOURCES 
copy_grid.c	     sep3d.f90 getkeyfirstbyte.c	getvalheaders.c		putnumberkeys.c
copykeys.c	     getkeyfmt.c	handle_file_pointers.c	putvalbindex.c
extractvalbyindex.c  getkeyindex.c	init_3d.c		putvalbname.c
fortran_junk.c	     getkeyname.c	insertvalbyindex.c	putvalheaders.c
  getkeytype.c	putaxis.c		reorder2.c
getaxis.c	     getnumberaxes.c	putgetgridwindow.c	reorder.c
getgridtag.c	     getnumberkeys.c	putgridtag.c		sep_3d_close.c
getheaderbytes.c     getvalbindex.c	putheadertag.c		tag_is_pipe.c
getheadertag.c	     getvalbname.c	putkey.c)



add_library(sep3d ${SOURCES})
TARGET_LINK_LIBRARIES(sep3d sep)
install(TARGETS sep3d DESTINATION lib EXPORT SEPlibExports)
