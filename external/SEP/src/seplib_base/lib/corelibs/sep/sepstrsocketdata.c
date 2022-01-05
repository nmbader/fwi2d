/*
opensocketstream


USAGE
opensocketstream(streaminf *info)

DESCRIPTION
Opens up dataset of type STREAMSOCKOUT
*/
/*
*
 *
 * Author: Bob 12-17-97
 *
 */
#include <sepConfig.h>
#include <stdio.h>
#include <fcntl.h>

#include <assert.h>


#ifdef RS6000
#undef __STR__
#endif

#include <string.h>

#ifndef _NFILE
#if defined(vax) || defined(sun)
#include <sys/param.h>
#define _NFILE NOFILE
#endif
#if defined(CONVEX) 
#include <limits.h>
#define _NFILE OPEN_MAX
#endif
#endif

#include <rpc/types.h>

#include "streamlist.h"
#include "sep_main_internal.h"




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void open_socketstream( streaminf *info )
_XFUNCPROTOEND 
#else
void open_socketstream( info )
     streaminf* info;
#endif
{
    
    assert( info->entrytype == STREAMSOCKOUT );

	  info->dataname = (char *) malloc((int)strlen("follow_hdr")+1);
    strcpy(info->dataname,"follow_hdr");
   

    /* initialise the default I/O routines */
    init_io(info );
}
/*  $Id: sepstrsocketdata.c,v 1.2 2004/04/08 22:32:27 bob Exp $ */
