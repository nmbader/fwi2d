#
# Prog.rules.C
#
#
# Contains the rules for compiling and linking C programs

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

ifndef UCDEBUGFLAGS
UCDEBUGFLAGS = ${UCFLAGS}
endif

ifndef UCSINGLEFLAGS
UCSINGLEFLAGS = ${UCFLAGS}
endif

#
# We need create the flags that we are going to use for compiling.
# If the user has already set these flags in the makefile, don't bother.
# The order is important here.  
#

ifndef CDEBUGFLAGS
CDEBUGFLAGS =   ${UCINCLUDES} ${OCFLAGS} ${OCDEFINES} ${UCDEFINES}  ${OCDEBUGFLAGS} ${UCFLAGS}
endif
ifndef CSINGLEFLAGS
CSINGLEFLAGS =  ${UCINCLUDES} ${OCFLAGS} ${OCDEFINES} ${UCDEFINES} ${COPTIMFLAGS} ${UCSINGLEFLAGS}  
endif
ifndef CFLAGS
CFLAGS =  ${UCINCLUDES} ${OCFLAGS}  ${OCDEFINES}  ${UCDEFINES} ${COPTIMFLAGS} ${CPARFLAGS} ${UCFLAGS}  
endif

# No we need to define how we compile.  Default is a user override
# followed by DEBUG=yes, SINGLE=yes, object file in CDEGUG, object file
# in CSINGLE, generic compile.



ifndef C_COMPILE
ifeq (${DEBUG}, yes)
define C_COMPILE
	$(compile.init)
	${CC} ${CDEBUGFLAGS} $< -DSOURCE='"'${FULLSRC}/$(<F)'"' -o $@
endef
else #DEBUG FLAG SET TO NO
ifeq (${SINGLE}, yes)
define C_COMPILE
	$(compile.init) 
	${CC}  ${CSINGLEFLAGS} $< -DSOURCE='"'${FULLSRC}/$(<F)'"' -o $@
endef
else #SINGLE FLAG SET TO NO
ifeq ($(@F),  $(findstring $(@F),${CDEBUG}))
define C_COMPILE
	$(compile.init)
	${CC} ${CDEBUGFLAGS} $< -DSOURCE='"'${FULLSRC}/$(<F)'"' -o $@
endef
else #FILE NOT FOUND IN CDEBUG LIST
ifeq ($(@F),  $(findstring $(@F),${CSINGLE}))
define C_COMPILE
	$(compile.init)
	${CC} ${CSINGLEFLAGS} $< -DSOURCE='"'${FULLSRC}/$(<F)'"' -o $@
endef
else #FILE NOT IN CSINGLE LIST
define C_COMPILE
	$(compile.init)
	${CC}   ${CFLAGS} $< -DSOURCE='"'${FULLSRC}/$(<F)'"' -o $@ 
endef
endif #CSINGLE LIST
endif #CDEBUG LIST
endif #SINGLE = yes
endif #DEBUG=y
endif #C_COMPLE NOT DEFINED
#
#
#---------------------------------------------------------------------------





#---------------------------------------------------------------------------
#                                    Rules
#

#NOW THE RULES FOR COMPILING
${OBJDIR}/%.o: ${SRCDIR}/%.${CEXT}
	${C_COMPILE}

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

ifndef  CLIBDIRS
CLIBDIRS = $(addprefix -L,${UCLIBDIRS} ${SITE_LIBDIRS} ${OCLIBDIRS})
endif

ifndef CLIBS
	CLIBS = ${UCLIBS} ${LOCALDEFAULTLIBS} ${OCLIBS}
endif



ifndef CFLGAGS
	CLDFLAGS  = ${OCLD_FLAGS_BEFORE} $(filter %.o,$^) ${OCLD_FLAGS_AFTER}
endif


#NOW THE LINKING INFORMATION
ifndef C_LN
define C_LN
	${CLD} ${CLDFLAGS} ${CLIBDIRS} ${CLIBS}
	${INSTALL_AOUT}
endef
endif


#
#
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
#                                    Definitions
#

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.${CEXT}
	${C_LN}
