/*
 * LOGIC:
 *         Basic module for opening and closing sepstream objects.
 *         External entries defined here are:
 *
 *     sepstr_open( info )          initialize header and data stream.
 *     sepstr_copyh( infin, infout ) copy input header to output header.
 *     sepstr_hclose( info )         close header file.
 *
 * Author:  D. Nichols  1-8-93   SEP
 * Revised: S.A. Levin  5-6-95   Mobil/SEP don't use non-posix strdup
 * Revised: B. Biondi  29-7-95   commented error when does'nt find in=
 * Revised: R. Clapp   18-7-97   Prototypes
 *
 *
 */
#include <sepConfig.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>

#include <assert.h>

#include "streamlist.h"

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



char parambuf2[80123];

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void open_inoutstream( streaminf *info )
_XFUNCPROTOEND 
#else
void open_inoutstream( info )
streaminf* info;
#endif
{
    
    assert( info->entrytype == STREAMINOUT||  info->entrytype == STREAMSCR  );
/*    assert( info->headfile  != (FILE*)0 );*/

   
    if( info->headerbuf != 0 ){
	/* reading from existing header, go and find the input name */ 
/* commented by Biondo. Dave would probably disagree */
/*
fprintf(stderr,"open_inoutstream:  non existing Header file \n");
	if( sepstrpar( info, "in", "s", parambuf2 ) == 0  ){
            seperr( "unable to obtain in=... for tag \"%s\"\n",info->tagname);
        }
*/
	if( sepstrpar( info, "in", "s", parambuf2 ) == 0  ){
	        outname( info );
        }else if(0==strcmp(parambuf2,"-1"))   outname( info );
        else{ 
	  info->dataname = (char *) malloc((int)strlen(parambuf2)+1);
	  strcpy(info->dataname, parambuf2);
      }
	
    }else{
	/* creating a new output, go and make up the name */
        outname( info );
    }

    /* initialise the default I/O routines */
    init_io(info );
}
