#ifndef _SEP_CUBE_H
#define _SEP_CUBE_H

/* Include file <sepcube.h>  This file contains function definitions for all
 * the routines in libcube. It should be includded in all the c-files in
 * that library.
 *
 * Author 10/23/91 D.Nichols - split up seplib.h
 *
 *  copyright (c) 1991 Stanford University
 *
 */

#include <sepConfig.h>
#include <stdio.h>
#include <prototypes.h>

/* Moved this to the general include file because almost everything uses
   malloc - questionable  -Bob*/

#if defined(HAVE_STDLIB_H)
#include <stdlib.h>
/* library name redefinitions, to prevent name conflicts  */
#if defined(CRAY)
#define _lib_sep(name) sep_##name
#define head _lib_sep(head)
#define input _lib_sep(input)
#define output _lib_sep(output)
#endif
#else
extern char *malloc ();
extern char *realloc ();
extern void free ();
#endif
#endif

#if defined(__STDC__) && defined(CRAY)
/* typedef int FILE_DESCR; */
#ifndef FILE_DESCR
#define FILE_DESCR int
#endif
#endif



/* typedef int FILE_DESCR; */
#ifndef FILE_DESCR
#define FILE_DESCR int
#endif



#include<sep_fortran.h>        
#include<sep_old_external.h>
#include<sep_main_external.h>
#include<sep_pars_external.h>
