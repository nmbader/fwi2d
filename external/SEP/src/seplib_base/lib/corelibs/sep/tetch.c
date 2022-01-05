/*
 * tetch
 * history
 * kamal		7/18/1986	first version: conflicted with hetch
 * kamal		7/21/1986	implemented stew's suggestion of clearing
 *					up the conflict with hetch
 * kamal		7/23/1986	make it similar to fetch instead of
 *					hetch (i.e. do getch first).
 * stew			7/23/1986	simplified chopping logic using regex(3)
 * kamal		7/25/1986	make it return pointer to top of history
 *					instead of the end if no match is found.
 * stew			9/7/87		<varargs> for portability
 * dave			9/16/90		made ansi-c and posix compatible.
 * dave 		9/17/90  	Use stdarg for ANSI-C compilers
 */
#include <sepConfig.h>
#include <stdlib.h>
#include "fastpar.h"
#include <string.h>
#include <sep_pars_external.h>
#include "streamlist.h"



#if NeedFunctionPrototypes
#include <stdarg.h>
#else
#include <varargs.h>
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
/*VARARGS2*/
int tetch( char *tag, char *type, ... )
_XFUNCPROTOEND
#else
/*VARARGS0*/
int tetch(va_alist /*tag,type,ptr*/)
va_dcl
#endif
{
MIXED var;
#if NeedFunctionPrototypes
#else
char *tag, *type;
#endif
int rc;
va_list apdum;
static  streaminf *info;

#if NeedFunctionPrototypes
 va_start(apdum,type);
#else
va_start(apdum);
tag = va_arg(apdum,char *);
type = va_arg(apdum,char *);
#endif
switch(type[0]) {
case 'g': var.g = va_arg(apdum,double *);
	  rc = getch(tag,type,var.g); break;
case 's': var.s = va_arg(apdum,char *);
	  rc = getch(tag,type,var.s); break;
case 'f': case 'r': var.f = va_arg(apdum,float *);
	  rc = getch(tag,type,var.f); break;
default: var.i = va_arg(apdum,int *);
	  rc = getch(tag,type,var.i); break;
}
va_end(apdum);


if(0 < rc) return(rc);

if( info == 0 ) info = tag_info( "in", TAG_IN );

if( info->tqlen == 0 ) return 0;

return(getpar_decode(info->tetch_queue,info->tqlen,tag,type,var));

}
