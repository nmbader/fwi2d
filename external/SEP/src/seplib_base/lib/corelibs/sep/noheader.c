/*
 * noheader() - looks for noheader=y on command line
 *		returns 1 if found , 0 if not.
 *
 * Author S. Levin  3/26/83
 */
#include <sepConfig.h>
#include <sepcube.h>
#include "sep_main_internal.h"
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int noheader(void)
_XFUNCPROTOEND
#else
int noheader()
#endif
{
 char temp[100];

 if(getch("noheader","s",temp))
     if(temp[0] == 'y') return(1);
 return(0);
}
