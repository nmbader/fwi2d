#include <prototypes.h>
#include <stdio.h>


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern void finitpar(void);
extern int noieeeinterupt(void);
extern int fortalloc( char*,int);
extern int fortfree( char*);
extern char* fsbrk(unsigned int);
_XFUNCPROTOEND
#else /*END OF NO PROTO */
extern int noieeeinterupt();
extern void finitpar();
extern int fortalloc();
extern int fortfree();
extern char*fsbrk();
#endif /*END OF NO PROTO */

/*  $Id: sep_fortran_internal.h,v 1.1.1.1 2004/03/25 06:37:24 cvs Exp $ */
/*  $Id: sep_fortran_internal.h,v 1.1.1.1 2004/03/25 06:37:24 cvs Exp $ */
