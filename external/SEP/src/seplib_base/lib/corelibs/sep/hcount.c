/*
hcount

USAGE
int hcount(flag)

INPUT PARAMETERS
flag - int dummy argument

RETURN VALUES
 n = number of times hclose has been called

*/
/*
 * 
 * utility subroutine to check whether hclose has been called already
 * return value is number of previous hclose calls.
 * 
 */
#include <sepConfig.h>
#include <stdio.h>

#include <sep_main_external.h>
#include "streamlist.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int hcount(int flag)
_XFUNCPROTOEND
#else
int hcount(flag)
int flag;
#endif
{
streaminf *info;

 info = tag_info( "out", TAG_OUT );
 return( info->ready_out );

}
