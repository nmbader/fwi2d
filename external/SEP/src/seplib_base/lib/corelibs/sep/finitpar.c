/* 
 * finitpar()
 *
 * initialize getpar global information.
 * used in Fortran programs. The fortran wrappers that are generated mean
 * that fortran programs call this routine as initpar().
 *
 * Author: Stewart A. Levin MRDC  4/14/89
 * Revised: Dave Nichols, SEP 10/27/90
 * 	use SEP style defines and added support for Dec-3100
 * W. Bauske IBM 04-21-91
 *	The code that uses FORTRAN getarg() allocated 1 too few bytes for
 *	each string it retrieved. 
 * Revised: Stewart A. Levin SEP 12-5-93
 *      switch to sepxargc, sepxargv to avoid shared library linkage
 *      problems
 * Revised: Martin Karrenbach, Karlsruhe Univ. 4-1-95
 *      added support for non-architecture specific  F90 compiler: NAG
 * Revised: Martin Karrenbach, Karlsruhe Univ. 9-1-95
 *      added support for non-architecture specific  HPF compiler: PGI
 *      uses a non-documented calling sequence of getarg_
 *      based on information from PGI directly (tested on SGI Powerchallenge)
 * Revised: Martin Karrenbach, Karlsruhe Univ. 9-7-95
 *      added support for non-architecture specific  GNU compiler: G77
 * Revised: Martin Karrenbach, Karlsruhe Univ. 11-17-95
 *      added support for INTEL Paragon
 * Revised: Robert Clapp 7-18-97 Moved stdlib include to sepcube
 */

#include <sepConfig.h>
#include "sep_fortran_internal.h"
#include <cfortran.h>

/* use Fortran library routine anyhow for shared library safety
 * #if defined SUN || defined sun
 * #define XARG
 * #endif 
 */

#if defined(CONVEX) 
#define XARG
#endif /* CONVEX */

#ifdef XARG

/* everything is predefined for some machines */

extern int xargc; extern char **xargv;
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void finitpar(void) 
_XFUNCPROTOEND
#else
void finitpar() 
#endif
{initpar(xargc,xargv);}

#else /*XARG*/

/* these machines don't have xargv and xargc predefined so we have to
  find out what they are */

/* some machines just need to copy some other external variable */

#if defined(DEC) || defined(DECALPHA)

/* DEC mips based workstation */
extern int f77argc; extern char** f77argv;

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void finitpar(void)
_XFUNCPROTOEND
#else
void finitpar()
#endif
{
 initpar(f77argc,f77argv);
}
#define FIXED
#endif

#if defined(SGI) & ( !defined(pgiFortran) )

/* SGI workstations without special compiler */ 
extern int f77argc; extern char** f77argv;

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void finitpar(void)
_XFUNCPROTOEND
#else
void finitpar()
#endif
{
 initpar(f77argc,f77argv);
}
#define FIXED
#endif

#ifdef ARDENT
extern int _UT_argc; extern char** _UT_argv;

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void finitpar(void)
_XFUNCPROTOEND
#else
void finitpar()
#endif
{
 initpar(_UT_argc,_UT_argv);
}

#define FIXED
#endif

#ifndef FIXED

/* The information was not available elsewhere try and pick it up 
 * using the fortran library routines.
 */

#include <stdio.h>
static char *cmdline;

#if defined(NAGf90Fortran) || defined(pgiFortran) || defined(gnu) || defined(IBMR2Fortran) || defined(gfortran)

#if defined(NAGf90Fortran)
#define iargc_ f90_unix_MP_iargc
#define getarg_ f90_unix_MP_getarg
#endif

#if defined(gfortran)
#define iargc_ _gfortran_iargc
#define getarg_ _gfortran_getarg_i4
#endif

#ifdef F90_TYPE_PGI
/* uses the normal getargc and iargc system routines */
#endif

#else

#ifdef AbsoftProFortran
#define iargc_ IARGC
#define getarg_ GETARG
#endif

#ifdef CRAY
#define iargc_ IARGC
#define getarg_ GETARG
#endif

#ifdef HP700
#define iargc_ FTN_300_IARGC
#define getarg_ Ftn_getarg
#endif



#endif

/* define how the getarg routine indexes the command line values */
/* the fist argument ( the program name can be index 1 or zero ) */
#if defined(NAGf90Fortran) 
#define INDEX_START 0
#else
#ifdef HP700 
/*#define INDEX_START 1 */
#define INDEX_START 0
#else
#define INDEX_START 0
#endif
#endif


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void finitpar(void)
_XFUNCPROTOEND
#else
void finitpar()
#endif
{
 extern int iargc_();  /* Fortran library function */
#ifdef CRAY 
 extern int getarg_(); /* Fortran library function */
#else
 extern void getarg_(); /* Fortran library function */
#endif
 char temp[256];
 int xargc; char **xargv;
 int i,j,idx,len;

#ifdef F90_TYPE_PGI
     int pghpf_0;
#endif
 xargc=iargc_();

#ifdef F90_TYPE_PGI
 xargc--;
#endif

 xargv = (char **) calloc(xargc+2,sizeof(char *));
 if(xargv == ((char **) NULL)) seperr ("finitpar: calloc failed\n");
 xargv[xargc+1] = (char *) NULL; /* standard null terminator */

 /* retrieve command line entries using Fortran library routine */
 for(i=0; i<=xargc; i++) {
     for(j=254; j>=0; j--) temp[j] = ' ';
     idx = i + INDEX_START;

#ifdef CRAY 
     len=32; /* number of 8byte integers in temp */
     j = getarg_(&idx,temp,&len);
#else
     len=256; /* number of bytes in temp */

#ifdef F90_TYPE_PGI
     pghpf_0=0;
     getarg_(&idx,temp,&pghpf_0,&pghpf_0,&len); 
#else
#if defined (HP700)
     getarg_(&idx,temp,&len); 
#else
     getarg_(&idx,temp,len);
#endif
#endif

#endif

     for(j=254; j>=0; j--) if(temp[j] != ' ' && temp[j] != '\0') break;
     xargv[i] = (char *) calloc(j+2,sizeof(char));

     if(xargv[i] == ((char *) NULL)) seperr ("finitpar: calloc failed\n");
     if(j>=0) strncpy(xargv[i],temp,j+1); /* copy and null terminate */
     xargv[i][j+1] = '\0';

     }
 xargc++;
 initpar(xargc,xargv);
}
#endif /* FIXED */

#endif /* XARG */

/* now make a fortran biding to the name initpar */
FCALLSCSUB0(finitpar,INITPAR,initpar)

