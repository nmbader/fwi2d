#
# Prog.rules.F90
#
#
# Contains the rules for compiling and linking F90 programs

#============================================================================
#============================================================================
#                                     COMPILING
#============================================================================
#============================================================================


#The general idea is that you can specify within a single makefile
#what to compile in debug mode, single processor, and multi processor 
#(if available).  By setting DEBUG to yes or SINGLE to yes, you describe
#how ALL object files will be created.  By specifying the names of
#certain object files in UCDEBUG and UCSINGLE you can have a greater
#degree of control.
#

#---------------------------------------------------------------------------
#                                    Definitions
#

# User has the ability to specify different flags for programs to
# be compiled in debug and single mode.  If not specified default to
# the more general user defined flags

ifndef UF90DEBUGFLAGS
UF90DEBUGFLAGS = ${UF90FLAGS}
endif

ifndef UF90SINGLEFLAGS
UF90SINGLEFLAGS = ${UF90FLAGS}
endif

#
# We need create the flags that we are going to use for compiling.
# If the user has already set these flags in the makefile, don't bother.
# The order is important here.  
#


ifndef F90DEBUGFLAGS
F90DEBUGFLAGS =  ${UF90INCLUDES}  ${OF90FLAGS} ${OF90DEBUGFLAGS} ${UF90FLAGS}
endif
ifndef F90SINGLEFLAGS
F90SINGLEFLAGS = ${UF90INCLUDES} ${OF90FLAGS} ${F90OPTIMFLAGS} ${UF90SINGLEFLAGS}  
endif
ifndef F90FLAGS
F90FLAGS = ${UF90INCLUDES}  ${OF90FLAGS}  ${F90OPTIMFLAGS}  ${F90PARFLAGS} ${UF90FLAGS}  
endif

RATFLAGS90 = -SOURCE ${FULLSRC}/$(<F)
F90PRE = ${FPP} ${UF90DEFINES} ${OF90DEFINES} 
LOPTRAN = ${SEPBINDIR}/lop2f90
CLOPTRAN = /james/bin/all/clop2f90
DCLOPTRAN = /james/bin/all/dclop2f90

# No we need to define how we compile.  Default is a user override
# followed by DEBUG=yes, SINGLE=yes, object file in F90DEGUG, object file
# in F90SINGLE, generic compile.


#ifndef F90_COMPILE
#ifeq (${DEBUG}, yes)
#define F90_COMPILE
#	$(compile.init)
#	${F90C}   ${F90DEBUGFLAGS}   -o $@
#endef
#else #DEBUG FLAG SET TO NO
#ifeq (${SINGLE}, yes)
#H=$(strip $(findstring $*,${F90PAR}))
#ifneq "$H" ""
#define F90_COMPILE
#	$(compile.init)
#	${F90C}  ${F90FLAGS}   -o $@
#endef
#else
#define F90_COMPILE
#	$(compile.init)
#	${F90C}  ${F90SINGLEFLAGS}   -o $@
#endef
#endif
#else #SINGLE FLAG SET TO NO
#ifneq ($(strip $(findstring $*,${F90DEBUG})),)
#define F90_COMPILE
#	$(compile.init)
#	${F90C}  ${F90DEBUGFLAGS} -o $@
#endef
#else #FILE NOT FOUND IN F90DEBUG LIST
#ifneq ($(strip $(findstring $*,${F90SINGLE})),)
#define F90_COMPILE
#	$(compile.init)
#	${F90C}   ${F90SINGLEFLAGS}   -o $@
#endef
#else #FILE NOT IN F90SINGLE LIST
#define F90_COMPILE
#	$(compile.init)
#	${F90C} ${F90FLAGS}  -o $@ 
#endef
#endif #F90SINGLE LIST
#endif #F90DEBUG LIST
#endif #SINGLE = yes
#endif #DEBUG=y
#endif #F90_COMPLE NOT DEFINED
##
##
##---------------------------------------------------------------------------
#
ifeq (DEBUG, yes)  #IF IN DEBUG MODE
FLAGS = ${F90DEBUGFLAGS} bb
else
ifeq (SINGLE, no) 
FLAGS = ${F90FLAGS}
else
BASE=  ${UF90INCLUDES}  ${OF90FLAGS}  ${F90OPTIMFLAGS}
FLAGS=  ${BASE} $(patsubst %,${F90PARFLAGS},$(strip $(findstring $*,${F90PAR}))) ${UF90FLAGS} $(patsubst %,${UF90PARFLAGS},$(strip $(findstring $*,${F90PAR})))
endif
endif

ifndef F90_COMPILE
define F90_COMPILE
	$(compile.init)
	${F90C}   ${FLAGS}   -o $@
endef
endif





#---------------------------------------------------------------------------
#                                    Rules
#

#NOW THE RULES FOR COMPILING

#not sure if these should be SRCDIR

#If RATF90 = yes then we use the f90 compiler

ifeq (${RATF90},yes)

${OBJDIR}/%.o: ${SRCDIR}/%.rst
	${RATFOR90}  -sep  ${RATFLAGS90} <$< >$*.${F90EXT}
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT} 
endif

${OBJDIR}/%.o: ${SRCDIR}/%.rs
	${RATFOR90}  -sep  ${RATFLAGS90} <$< >$*.${F90EXT}
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT} 
endif

${OBJDIR}/%.o: ${SRCDIR}/%.rt
	${RATFOR90}   <$< >$*.${F90EXT}
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT} 
endif

${OBJDIR}/%.o: ${SRCDIR}/%.r
	${RATFOR90}   <$< >$*.${F90EXT}
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT}
endif

${OBJDIR}/%.o: ${SRCDIR}/%.RS
	${F90PRE} <$< | \
	${RATFOR90}   -sep  |\
	sed s/"DIR NOBOUNDS"// | \
	sed s/auxpar,auxputch/auxpar/  | ${F90PRE} > $*.${F90EXT}
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT} 
endif

endif #END OF RATF90=yes

${OBJDIR}/%.o: ${SRCDIR}/%.lop
	${RATFOR90}   <$< >$*.${F90EXT} -sep
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT} 
endif


${OBJDIR}/%.o: ${SRCDIR}/%.clop
	${RATFOR90}   <$< >$*.${F90EXT} -sep -complex
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT} 
endif


${OBJDIR}/%.o: ${SRCDIR}/%.dclop
	${RATFOR90}   <$< >$*.unrat -sep -dcomplex
	<$*.unrat  ${F90PRE}  >$*.${F90EXT}
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT} 
endif


%.${F90EXT}:	${SRCDIR}/%.lop
	${RATFOR90}   <$< >$*.${F90EXT} -sep

%.${F90EXT}:	${SRCDIR}/%.dclop
	${RATFOR90}   <$< >$*.unrat -sep -dcomplex
	<$*.unrat  ${F90PRE}  >$*.${F90EXT}

%.${F90EXT}:	${SRCDIR}/%.clop
	${RATFOR90}   <$< >$*.unrat -sep -complex
	<$*.unrat  ${F90PRE}  >$*.${F90EXT}
	@${RM} $*.unrat


%.${F90EXT}:	${SRCDIR}/%.r90
	${RATFOR90}   <$< >$*.${F90EXT} -sep


#${OBJDIR}/%.o: ${SRCDIR}/%.clop
#	${RATFOR90}   <$< >$*.unrat
#	<$*.unrat  ${F90PRE} | ${CLOPTRAN} >$*.${F90EXT}
#	${F90_COMPILE}  $*.${F90EXT}
#ifneq (${SAVE_FORTRAN}, yes)
#	@${RM} $*.${F90EXT}
#endif

#${OBJDIR}/%.o: ${SRCDIR}/%.dclop
#	${RATFOR90}   <$< >$*.unrat
#	<$*.unrat  ${F90PRE} | ${DCLOPTRAN} >$*.${F90EXT}
#	${F90_COMPILE}  $*.${F90EXT}
#ifneq (${SAVE_FORTRAN}, yes)
#	@${RM} $*.${F90EXT} 
#endif

${OBJDIR}/%.o: ${SRCDIR}/%.img
	${RATFOR90}   <$< >$*.unrat
	<$*.unrat  ${F90PRE} | ${LOPTRAN} >$*.${F90EXT}
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT} 
endif

${OBJDIR}/%.o: ${SRCDIR}/%.r90
	${RATFOR90}   <$< >$*.${F90EXT}
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT} 
endif

${OBJDIR}/%.o: ${SRCDIR}/%.rs90
	${RATFOR90}  -sep  ${RATFLAGS90} <$< >$*.${F90EXT}
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT} 
endif


${OBJDIR}/%.o: ${SRCDIR}/%.FS90
	${SAWF90} $< | sed s/C#/#/ | sed s/"DIR NOBOUNDS"// | \
	sed s/auxpar,auxputch/auxpar/  | ${F90PRE} > $*.%{FEXT}
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT}
endif






#aa=$(patsubst %.PPP,%,$(strip $(findstring $*,${F90PAR})))
#bb=$(patsubst %,F90P,$(strip $(findstring $*,${F90PAR})))

#ifeq (${bb},CCC)
#aa=yes
#else
#aa=no
#endif


${OBJDIR}/%.o:	${SRCDIR}/%.F
	${F90PRE} -DSOURCE='"'${FULLSRC}/$(<F)'"' < ${SRCDIR}/$*.F >$*.fix.${F90EXT}
	${F90_COMPILE} ${aa} ${bb}   $*.fix.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.fix.${F90EXT}
endif

${OBJDIR}/%.o:	${SRCDIR}/%.F90
	${F90PRE} -DSOURCE='"'${FULLSRC}/$(<F)'"' < ${SRCDIR}/$*.F90 >$*.fix.${F90EXT}
	${F90_COMPILE} ${aa} ${bb}   $*.fix.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.fix.${F90EXT}
endif

ifneq (${PREPROCESS}, no)
${OBJDIR}/%.o:	${SRCDIR}/%.${F90EXT}
	${F90PRE} -DSOURCE='"'${FULLSRC}/$(<F)'"' < ${SRCDIR}/$*.${F90EXT} >$*.fix.${F90EXT}
	${F90_COMPILE} ${aa} ${bb}   $*.fix.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.fix.${F90EXT}
endif
else #Don't run preprocessor
${OBJDIR}/%.o: ${SRCDIR}/%.${F90EXT}
	${F90_COMPILE}  ${SRCDIR}/$*.${F90EXT}
endif

${OBJDIR}/%.o: ${SRCDIR}/%.fs90
	${SAWF90} $<  | ${F90PRE} > $*.${F90EXT}
	${F90_COMPILE}  $*.${F90EXT}
ifneq (${SAVE_FORTRAN}, yes)
	@${RM} $*.${F90EXT}
endif



#
#
#---------------------------------------------------------------------------






#============================================================================
#============================================================================
#                                     LINKING
#============================================================================
#============================================================================




#---------------------------------------------------------------------------
#                                    Definitions
#



# We are building the link step in parts
# because the user might want to override it for certain programs
# and we want to make it as convenient as possible.

ifndef  F90LIBDIRS
F90LIBDIRS = $(addprefix -L,${UF90LIBDIRS} ${SITE_LIBDIRS} ${OF90LIBDIRS})
F9077LIBDIRS = $(addprefix -L,${UF77LIBDIRS} ${SITE_LIBDIRS} ${OF90LIBDIRS})
endif

ifndef F90LIBS
	F90LIBS = ${UF90LIBS} ${LOCALDEFAULTLIBS} ${OF90LIBS}
	F9077LIBS = ${UF77LIBS} ${LOCALDEFAULTLIBS} ${OF90LIBS}
endif

ifndef F90LDFLAGS
	F90LDFLAGS  = ${OF90LD_FLAGS_BEFORE} $(filter %.o,$^) ${OF90LD_FLAGS_AFTER}
endif

#NOW THE LINKING INFORMATION
ifndef F90_LN
define F90_LN
	${F90LD} ${F90LDFLAGS} ${F90LIBDIRS} ${F90LIBS}
	${INSTALL_AOUT}
endef
endif

#Special rules when running ratfor90 preprocessor over F77 code
ifndef F9077_LN
define F9077_LN
	${F90LD} ${F90LDFLAGS} ${F90IBDIRS} ${F9077LIBS}
	${INSTALL_AOUT}
endef
endif

#
#
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
#                                    Definitions
#
ifeq (${RATF90},yes)
${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.rst
	${F9077_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.rt
	${F9077_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.rs
	${F9077_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.r
	${F9077_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.RS
	${F9077_LN}

endif #END OF RATF90=yes

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.r90
	${F90_LN}

%.x ${BINDIR}/%.x:	${OBJDIR}/%.o ${SRCDIR}/%.rs90
	${F90_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.FS90
	${F90_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.r90
	${F90_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.F
	${F90_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.F90
	${F90_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.${F90EXT}
	${F90_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.fs90
	${F90_LN}

junk:
	echo ${MAKE_DEPEND}

ifeq (${MAKE_DEPEND},yes)
ifndef SOURCES
SOURCES= $(notdir $(wildcard ${SRCDIR}/*.r90)  $(wildcard ${SRCDIR}/*.F90) $(wildcard ${SRCDIR}/*.rs90) $(wildcard ${SRCDIR}/*.f90) $(wildcard ${SRCDIR}/*.lop) $(wildcard ${SRCDIR}/*.dclop) $(wildcard ${SRCDIR}/*.clop) $(wildcard ${SRCDIR}/*.F90))
endif
-include .make.dependencies.${SEP_ARCH}
.make.dependencies.${SEP_ARCH}: $(addprefix ${SRCDIR}/, ${SOURCES})
	@${SEPBINDIR}/Makedepend -s ${SRCDIR} -x ${BINDIR} -d ${OBJDIR} ${SOURCES} >$@
endif

