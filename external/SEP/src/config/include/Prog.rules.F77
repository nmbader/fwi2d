#
# Prog.rules.F77
#
#
# Contains the rules for compiling and linking F77 programs

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

ifndef UF77DEBUGFLAGS
UF77DEBUGFLAGS = ${UF77FLAGS}
endif

ifndef UF77SINGLEFLAGS
UF77SINGLEFLAGS = ${UF77FLAGS}
endif

#
# We need create the flags that we are going to use for compiling.
# If the user has already set these flags in the makefile, don't bother.
# The order is important here.  
#

ifndef F77DEBUGFLAGS
F77DEBUGFLAGS =  ${UF77INCLUDES} ${OF77FLAGS}  ${OF77DEBUGFLAGS} ${UF77FLAGS}
endif
ifndef F77SINGLEFLAGS
F77SINGLEFLAGS = ${UF77INCLUDES} ${OF77FLAGS} ${F77OPTIMFLAGS} ${UF77SINGLEFLAGS}  
endif
ifndef F77FLAGS
F77FLAGS = ${UF77INCLUDES} ${OF77FLAGS}   ${F77OPTIMFLAGS} ${F77PARFLAGS} ${UF77FLAGS}  
endif

F77PRE = ${FPP} ${UF77DEFINES} ${OF77DEFINES} -Dsource='"'${FULLSRC}/$(<F)'"'

# No we need to define how we compile.  Default is a user override
# followed by DEBUG=yes, SINGLE=yes, object file in F77DEGUG, object file
# in F77SINGLE, generic compile.


ifndef F77_COMPILE
ifeq (${DEBUG}, yes)
define F77_COMPILE
	$(compile.init)
	${F77C} ${F77DEBUGFLAGS}   -o $@
endef
else #DEBUG FLAG SET TO NO
ifeq (${SINGLE}, yes)
define F77_COMPILE
	$(compile.init)
	${F77C} ${F77SINGLEFLAGS}   -o $@
endef
else #SINGLE FLAG SET TO NO
ifeq ($(@F),  $(findstring $(@F),${F77DEBUG}))
define F77_COMPILE
	$(compile.init)
	${F77C} ${F77DEBUGFLAGS} ${SRCDIR}/%.fix.f  -o $@
endef
else #FILE NOT FOUND IN F77DEBUG LIST
ifeq ($(@F),  $(findstring $(@F),${F77SINGLE}))
define F77_COMPILE
	$(compile.init)
	${F77C} ${F77SINGLEFLAGS}   -o $@
endef
else #FILE NOT IN F77SINGLE LIST
define F77_COMPILE
	$(compile.init)
	${F77C} ${F77FLAGS}  -o $@ 
endef
endif #F77SINGLE LIST
endif #F77DEBUG LIST
endif #SINGLE = yes
endif #DEBUG=y
endif #F77_COMPLE NOT DEFINED
#
#
#---------------------------------------------------------------------------





#---------------------------------------------------------------------------
#                                    Rules
#

#NOW THE RULES FOR COMPILING

#not sure if these should be SRCDIR

#If RATF90 = yes then we use the f90 compiler


ifneq (${RATF90},yes)

${OBJDIR}/%.o: ${SRCDIR}/%.rst
	${SAW}  $<  > $*.rt
	${SAT}  $*.rt | ${RATFOR}   | ${F77PRE} -C > $*.${F77EXT}
	${F77_COMPILE}  $*.${F77EXT}
	@${RM} $*.rt $*.Rst
ifneq (${SAVE_FORTRAN} , yes)
	@${RM} $*.rt $*.${F77EXT}
endif
	
${OBJDIR}/%.o: ${SRCDIR}/%.rs
	${SAW}  $< | ${RATFOR}   | ${F77PRE} -C > $*.${F77EXT}
	${F77_COMPILE}  $*.${F77EXT}
ifneq (${SAVE_FORTRAN} , yes)
	@${RM} $*.${F77EXT}
endif

${OBJDIR}/%.o: ${SRCDIR}/%.rt
	${SAT}  $< |  ${RATFOR}  | ${F77PRE} -C  >$*.${F77EXT}
	${F77_COMPILE}  $*.${F77EXT}
	@${RM} $*.Rt
ifneq (${SAVE_FORTRAN} , yes)
	@${RM} $*.${F77EXT}
endif
	
${OBJDIR}/%.o: ${SRCDIR}/%.r
	${RATFOR}  $<   | ${F77PRE} -C >$*.${F77EXT}
	${F77_COMPILE}  $*.${F77EXT}
	@${RM} $*.R
ifneq (${SAVE_FORTRAN} , yes)
	@${RM} $*.${F77EXT}
endif

${OBJDIR}/%.o: ${SRCDIR}/%.RS
	${SAW} $< | ${RATFOR} -C | sed s/C#/#/ | sed s/"DIR NOBOUNDS"// | \
	sed s/auxpar,auxputch/auxpar/  | ${F77PRE} -C   > $*.${F77EXT}
	${F77_COMPILE}  $*.${F77EXT}
	@${RM} $*.Rs
ifneq (${SAVE_FORTRAN} , yes)
	@${RM} $*.${F77EXT}
endif

endif #END OF RATF90=yes


${OBJDIR}/%.o: ${SRCDIR}/%.FS
	${SAW} $< | sed s/C#/#/ | sed s/"DIR NOBOUNDS"// | \
	sed s/auxpar,auxputch/auxpar/   | ${F77PRE} -C > $*.%{FEXT}
	${F77_COMPILE}  $*.${F77EXT}
ifneq (${SAVE_FORTRAN} , yes)
	@${RM} $*.${F77EXT}
endif

ifeq (${PREPROCESS}, yes)
${OBJDIR}/%.o:	${SRCDIR}/%.${F77EXT}
	${F77PRE} < ${SRCDIR}/$*.${F77EXT} >$*.fix.${F77EXT}
	${F77_COMPILE}  $*.fix.${F77EXT}
ifneq (${SAVE_FORTRAN} , yes)
	@${RM} $*.fix.${F77EXT}
endif
else #Don't run preprocessor
${OBJDIR}/%.o: ${SRCDIR}/%.${F77EXT}
	${F77_COMPILE} ${SRCDIR}/$*.${F77EXT}
endif
	
${OBJDIR}/%.o: ${SRCDIR}/%.fs
	${SAWF} $<  | ${F77PRE} -C  > $*.${F77EXT}
	${F77_COMPILE}  $*.${F77EXT}
ifneq (${SAVE_FORTRAN} , yes)
	@${RM} $*.${F77EXT}
endif

	
${OBJDIR}/%.o: ${SRCDIR}/%.ft
	${SAT} $<   | ${F77PRE} -C > $*.${F77EXT}
	${F77_COMPILE}  $*.${F77EXT}
ifneq (${SAVE_FORTRAN} , yes)
	@${RM} $*.${F77EXT}
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

ifndef  F77LIBDIRS
F77LIBDIRS = $(addprefix -L,${UF77LIBDIRS} ${SITE_LIBDIRS} ${OF77LIBDIRS})
endif

ifndef F77LIBS
	F77LIBS = ${UF77LIBS}  ${LOCALDEFAULTLIBS} ${OF77LIBS} 
endif

ifndef F77LDFLAGS
	F77LDFLAGS  = ${OF77LD_FLAGS_BEFORE} $(filter %.o,$^) ${OF77LD_FLAGS_AFTER}
endif

#NOW THE LINKING INFORMATION
ifndef F77_LN
define F77_LN
	${F77LD} ${F77LDFLAGS} ${F77LIBDIRS} ${F77LIBS}    
	${INSTALL_AOUT}
endef
endif

#
#
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
#                                    Definitions
#

ifneq (${RATF90},yes)
${BINDIR}/%.x %.x:	 ${OBJDIR}/%.o  ${SRCDIR}/%.rst
	${F77_LN}
	
${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.rs
	${F77_LN}

${BINDIR}/%.x %.x:	 ${OBJDIR}/%.o ${SRCDIR}/%.rt
	${F77_LN}
	
${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.r
	${F77_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.RS
	${F77_LN}

endif #END OF RATF90=yes


${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.FS
	${F77_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.${F77EXT}
	${F77_LN}

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.fs
	${F77_LN}
	
${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.ft
	${F77_LN}
