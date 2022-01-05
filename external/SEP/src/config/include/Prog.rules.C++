#
# Prog.rules.C++
#
#
# Contains the rules for compiling and linking C++ programs

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

ifndef UC++DEBUGFLAGS
UC++DEBUGFLAGS = ${UC++FLAGS}
endif

ifndef UC++SINGLEFLAGS
UC++SINGLEFLAGS = ${UC++FLAGS}
endif

#
# We need create the flags that we are going to use for compiling.
# If the user has already set these flags in the makefile, don't bother.
# The order is important here.  
#

ifndef C++DEBUGFLAGS
C++DEBUGFLAGS =  ${UC++INCLUDES} ${OC++FLAGS} ${OC++DEFINES}  ${OC++DEBUGFLAGS} ${UC++FLAGS}
endif
ifndef C++SINGLEFLAGS
C++SINGLEFLAGS = ${UC++INCLUDES} ${OC++FLAGS} ${OC++DEFINES} ${C++OPTIMFLAGS} ${UC++SINGLEFLAGS}  
endif
ifndef C++FLAGS
C++FLAGS = ${UC++INCLUDES} ${OC++FLAGS}  ${OC++DEFINES} ${C++OPTIMFLAGS} ${C++PARFLAGS} ${UC++FLAGS}  
endif

# No we need to define how we compile.  Default is a user override
# followed by DEBUG=yes, SINGLE=yes, object file in C++DEGUG, object file
# in C++SINGLE, generic compile.


ifndef C++_COMPILE
ifeq (${DEBUG}, yada)
define C++_COMPILE
	$(compile.init)
	${C++C} ${C++DEBUGFLAGS} $< -DSOURCE='"'${FULLSRC}/$(<F)'"' -o $@
endef
else #DEBUG FLAG SET TO NO
ifeq (${SINGLE}, yes)
define C++_COMPILE
	$(compile.init)
	${C++C} ${C++SINGLEFLAGS} $< -DSOURCE='"'${FULLSRC}/$(<F)'"' -o $@
endef
else #SINGLE FLAG SET TO NO
ifeq ($(@F),  $(findstring $(@F),${C++DEBUG}))
define C++_COMPILE
	$(compile.init)
	${C++C} ${C++DEBUGFLAGS} $< -DSOURCE='"'${FULLSRC}/$(<F)'"' -o $@
endef
else #FILE NOT FOUND IN C++DEBUG LIST
ifeq ($(@F),  $(findstring $(@F),${C++SINGLE}))
define C++_COMPILE
	$(compile.init)
	${C++C} ${C++SINGLEFLAGS}	$< -DSOURCE='"'${FULLSRC}/$(<F)'"' -o $@
endef
else #FILE NOT IN C++SINGLE LIST
define C++_COMPILE
	$(compile.init)
	${C++C} ${C++FLAGS} $< -DSOURCE='"'${FULLSRC}/$(<F)'"' -o $@ 
endef
endif #C++SINGLE LIST
endif #C++DEBUG LIST
endif #SINGLE = yes
endif #DEBUG=y
endif #C++_COMPLE NOT DEFINED
#
#
#---------------------------------------------------------------------------





#---------------------------------------------------------------------------
#                                    Rules
#

#NOW THE RULES FOR COMPILING
${OBJDIR}/%.o: ${SRCDIR}/%.${C++EXT}
	${C++_COMPILE}

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

ifndef  C++LIBDIRS
C++LIBDIRS = $(addprefix -L,${UC++LIBDIRS} ${SITE_LIBDIRS} ${OC++LIBDIRS})
endif

ifndef C++LIBS
	C++LIBS = ${UC++LIBS} ${LOCALDEFAULTLIBS} ${OC++LIBS}
endif

ifndef C++FLGAGS
	C++LDFLAGS  = ${OC++LD_FLAGS_BEFORE} $(filter %.o,$^) ${OC++LD_FLAGS_AFTER}
endif


#NOW THE LINKING INFORMATION
ifndef C++_LN
define C++_LN
	${C++LD} ${C++LDFLAGS} ${C++LIBDIRS} ${C++LIBS}
	${INSTALL_AOUT}
endef
endif

#
#
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
#                                    Definitions
#

${BINDIR}/%.x %.x:	${OBJDIR}/%.o ${SRCDIR}/%.${C++EXT}
	${C++_LN}
