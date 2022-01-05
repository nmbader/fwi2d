/*
 head

 USAGE
	FILE *sep_head()


DESCRIPTION
Returns a file descriptor for the output header pointed to by a head=
command line parameter or stdout if not defined. This is especially
useful to C programmers who don't need to go through the somewhat
clumsy interface of putch. 

*/
 
/* Revised   4/2/83  S. Levin -- Added option for head=stdout.
 * Revised   7/14/83 S. Levin -- Added err and perror calls for diagnostics.
 * Revised   5/15/84 S. Levin -- Reopen if header closed
 * Revised   6/25/86 S. Levin -- Added option for head=stderr for segy use.
 */
#include <sepConfig.h>
#include <stdio.h>
#include <sep_main_internal.h>
#include <sep_main_external.h>
#include "streamlist.h"


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
FILE *sep_head(void)
_XFUNCPROTOEND
#else
FILE *sep_head()
#endif
{
 static FILE *headfile = NULL;
 streaminf* info;

 if(headfile == NULL || isclosed(headfile) ) 
    {
    info = tag_info( "out", TAG_OUT );
    headfile = info->headfile;
    }
 return(headfile);
}
