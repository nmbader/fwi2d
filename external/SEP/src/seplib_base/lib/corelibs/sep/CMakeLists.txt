cmake_minimum_required(VERSION 2.8)
project(SEPLIB-base-lib-corelibs-sep)

find_program(FLEX flex)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../include)


ADD_CUSTOM_COMMAND(OUTPUT sepgetpar_scan.c 
   COMMAND flex -l ${CMAKE_CURRENT_SOURCE_DIR}/sepgetpar_scan.l 
   COMMAND mv lex.yy.c sepgetpar_scan.c
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/sepgetpar_scan.l)



INCLUDE(TestBigEndian)

TEST_BIG_ENDIAN(BIGENDIAN)
IF(${BIGENDIAN})
     ADD_DEFINITIONS(-DWORDS_BIGENDIAN)
ENDIF(${BIGENDIAN})


set (SOURCES 	${CMAKE_CURRENT_BINARY_DIR}/sepgetpar_scan.c	 hcount.c	 puthead.c		   sepstrpar.c
auxclose.c	 head.c		 putlin.c		   sepstrput.c
auxin.c		 helixcart.c	 redin.c		   sepstrscrdata.c
auxpar.c	 hetch.c	 redout.c		   sepstrscrhead.c
auxputch.c	 init_io.c	 sepalloc.c		   sepstrsocketdata.c
auxputhead.c	 initpar.c	 separg.c		   sepstrsockethead.c
datapath.c	 input.c	 sepdoc.c		   sepstrsync.c
expandnm.c	 isapipe.c	 sep_err.c		   septagstream.c
fd_io.c		  sepfile.c		   sepwindow.c
fetch.c		 isclosed.c	 sepfsize.c		   slice.c
file_io.c	 isordinary.c	 sepgetch.c		   snap.c
findnm.c	 maketitle.c	 sepgetpar_decode.c	   sreed.c
	 make_unpipe.c	 	   srite.c
fortalloc.c	 math.c		 sepgetpar_string_store.c  sseek.c
fortran_junk.c	 multifd_io.c	 seploc.c		   ssize.c
	 noheader.c	 sepstrindata.c		   status.c
fsbrk.c		 noieeeint.c	 sepstrinhead.c		   streamlist.c
fullnm.c	 opensock.c	 sepstrinoutdata.c	   strformats.c
getnumberaxes.c  output.c	 sepstrinouthead.c	   tetch.c
getputaxis.c	 putch.c	 sepstroutdata.c	   warn.c
hclose.c	 putch_format.c  sepstrouthead.c)



set(SOURCES ${SOURCES} sep.f90       sep_func.f90  sep_mod.f90  sep_par.f90
 sep_io.f90    )


add_library(sep ${SOURCES})
install(TARGETS sep DESTINATION lib EXPORT SEPlibExports)
