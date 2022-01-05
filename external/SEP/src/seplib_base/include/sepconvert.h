#ifndef SEP_CONVERT_H
#define SEP_CONVERT_H 1

/* Include file for lib90aux  This file contains function definitions for all
 * the sep90 routines. It should be includded in all the c-files in
 * that library.
 *
 */

#include "prototypes.h"
#include<stdio.h>
#include<stdlib.h>
#ifdef __cplusplus
extern "C" {
#endif
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern int createfloat(int, float*, float*);
extern int getbackint(int, int*, int*);
extern int reed_char_real(char *, int, float*);
extern int rite_real_char(char *, int, float*);
extern int zasc(char*, char*,int);
_XFUNCPROTOEND
#else
extern int createfloat();
extern int getbackint();
extern int zasc();
extern int xdrbhdrsub();
extern int xdrhdrsub();
extern int reed_char_real();
extern int rite_real_char();
#endif
#ifdef __cplusplus
}
#endif

#endif
