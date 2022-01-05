cmake_minimum_required(VERSION 2.8)
project(SEPLIB-ratfor90)

set(DONOTHING
Prog.rules.C SEP.idoc.rules Prog.rules.C++ SEP.ifdef 
Dir.defs.bottom Doc.defs.bottom Prog.rules.F77 SEP.local.defs 
Prog.rules.F90 SEP.local.rules Doc.rules.action Prog.rules.bottom 
SEP.machine.type Doc.rules.book Prog.rules.depend SEP.obj.rules 
Doc.rules.clean Prog.rules.general SEP.overide.defaults Doc.rules.fig 
Prog.rules.java SEP.prog.rules Doc.rules.html Prog.rules.lib SEP.site.defs 
Doc.rules.idoc Prog.rules.shell Doc.rules.test README SEP.bottom 
SEP.distr.bottom SEP.distr.defs SEP.distr.dir Prg.rules.lib SEP.distr.rules 
sepweb.top SEP.distr.top Prog.defs.top septex.tar septex.pl sepbook.pl
not_found
    SEP.generic.defs
            SEP.lib.defs
   Doc.defs.top
  SEP.top
)


FIND_PROGRAM(PPMTOGIF ppmtogif)
if(PPMTOGIF) 
MESSAGE("IN HERE pptm")
   else(PPMTOGIF) 
  set(PPMTOGIF "do_nothing ppmtogif")
endif(PPMTOGIF)

find_program(PS2PDF epstopdf ps2pdf pstopdf)
if(PS2PDF)
else(PS2PDF)
  set(PS2PDF "do_nothing epstopdf")
endif(PS2PDF)



set(SEP ${CMAKE_INSTALL_PREFIX})
set(SEPBINDIR ${CMAKE_INSTALL_PREFIX}/bin)
set(SEPLIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
set(SEPINCDIR ${CMAKE_INSTALL_PREFIX}/include)
set(CC ${CMAKE_C_COMPILER})
set(CLD ${CMAKE_C_COMPILER})
set(F90C ${CMAKE_Fortran_COMPILER})
set(F90LD ${CMAKE_Fortran_COMPILER})
set(SYSCFLAGS ${CMAKE_C_FLAGS})
set(SYSF90FLAGS ${CMAKE_Fortran_FLAGS})

 
  set(min  "${CMAKE_CURRENT_SOURCE_DIR}/Prog.defs.compile.in")
  set(mout "${CMAKE_CURRENT_BINARY_DIR}/Prog.defs.compile")
  configure_file (${min} ${mout})

  set(xin  "${CMAKE_CURRENT_SOURCE_DIR}/local.defs.in")
  set(xout "${CMAKE_CURRENT_BINARY_DIR}/local.defs")
  configure_file (${xin} ${xout})

install(FILES ${DONOTHING} ${xout} ${mout}  DESTINATION include)
