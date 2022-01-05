#include <sepConfig.h>
/*
 interupts are disabled for improved speed 
11-12-90	Martin Karrenbach
*/

#ifdef DEC3100

/* include the DEC supplied fixup routines */ 
/* the included file needs the compile flag -DINEXACT
   see intr_handl_dec.c  for more details	*/
#define INEXACT
#include "intr_handl_dec.c"

#endif /* DEC3100 */

int noieeeinterupt()
{

#ifdef DEC3100

stp_flt_handl() ;
stp_int_handl() ;

#endif /* DEC3100 */

#ifdef SUN4

abrupt_underflow() ;

#endif /* SUN4 */

return 0 ;

}

/*  $Id: noieeeint.c,v 1.1.1.1 2004/03/25 06:37:24 cvs Exp $ */
