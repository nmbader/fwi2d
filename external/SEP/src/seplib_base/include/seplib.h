#ifndef _SEPLIB_H
#define _SEPLIB_H

/* Include file <seplib.h>  This file sets up 
 *
 * <stdio.h>
 * <math.h>
 * pi
 * complex and FILE_DESCR typedefs
 *
 * modified 3/26/83  S. Levin - added check for NOHEADER to fake out
 *   input
 * modified 4/1/83 S. Levin & Shuki - added MIN and MAX definitions.
 *		      Added ifndef checks to make sure these weren't
 *		      predefined. Changed complex and FILE_DESCR to
 *		      reduce potential conflict.
 * adapted  4/4/83 S. Levin - separated declarations into separate
 *		      seplib.h file.
 * modified 7/26/83 S. Levin - declared pi static to avoid linkage conflict
 * modified 5/15/84 S. Levin - added and corrected some extern defs.
 * modified 9/14/85 S. Levin - for Convex - convex had cabs and scabs in math.h
 * modified 9/18/90 D.Nichols - Added prototyped definitions for ansi-C.
 * modified 11/13/90 D.Nichols - reed and rite have void* buffers in ansi-c.
 * modified 10/23/91 D.Nichols - split into three files
 *      sepcube.h - defines the functions in libcube
 *      sepmath.h - defines the functions in libsepmath
 *      seplib.h  - standard setup file that pulls in all the header files
 *		    required for a standard setup.
 *
 *   copyright (c) 1991 Stanford University
 *
 */

#include <stdio.h>
#include <math.h>
#include <signal.h>

#ifdef MAXFLOAT
#define SEP_MAXFLOAT       MAXFLOAT
#else
#define SEP_MAXFLOAT       ((float)3.40282346638528860e+38)
#endif


#include <sepcube.h>
#include <sepmath.h>

#if defined(HAVE_STDLIB_H)
#include <time.h>
#else
extern
#ifdef __cplusplus
  "C"
#endif
long time();
#endif


#endif
