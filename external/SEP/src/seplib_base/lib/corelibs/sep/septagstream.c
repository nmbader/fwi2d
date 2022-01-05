/*
 * LOGIC:
 *         Maintain a simple list of streaminf structures.
 *         Predefined tags are "in", 
 *	   "out", "head".  List is reordered on the
 *         fly to put most recently called tags at head of
 *         list.  New tags are added as necessary with a
 *         call to sepstr_open  used to initialize
 *         the entry.  No facility at present to close the
 *         files to make space.
 *

 * Author: S. Levin     11/16/86   SEP
 * Revised: S. Levin    3/10/87   SEP  corrected srite_ in/out flag
 *                                     addressed problem of a program
 *                                     that never calls input().
 * Revised: S. Levin    5/7/87    SEP  made _NFILE machine dependent
 * Revised: S. Levin    3/5/91    MRDC/DRL worked around Convex compiler
 *					problem
 * Revised: W. Bauske   4-18-91   IBM  Removed use of bcopy for RS/6000
 *
 * Revised: D.Nichols   1-8-93    SEP major revision to use new sepstream
 *                                logic. "head" is a pseudo-tag that actually
 *                                opens the "in" file.
 * Revised: Bob          8-97     Added different prototypes
 * Revised: Bob         12-97     Added STREAMSOCKOUT
 */
#include <sepConfig.h>
#include <stdio.h>

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


#include "streamlist.h"
#include <sep_main_external.h>
    
/* return streaminf pointer given input tag. 0 on error */
/* second argument is 1 for read, 0 for write */

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
streaminf * tag_info(const char *tag, enum tagtype type)
_XFUNCPROTOEND
#else
streaminf * tag_info(tag,type)
char *tag; enum tagtype type;
#endif
{
    streaminf *curr;

    curr =sepstr_head();

    while( curr != 0 ){
   	if( strcmp(tag,curr->tagname) == 0 ) return curr;
		curr = curr->next;
   }

   if(type==TAG_INQUIRE) return(0);
    /* fell through so we must create it and put it at the end of the list*/
    curr = sepstr_new( (char *) tag, type );


    sepstr_addend( curr ); 


    switch( curr->entrytype ){
     case STREAMIN:
        sepstr_in_head( curr );
        break;
     case STREAMOUT:
        sepstr_out_head( curr );
        break;
     case STREAMINOUT:
        sepstr_inout_head( curr );
        break;
     case STREAMSOCKOUT:
        sepstr_socket_head( curr );
        break;
     case STREAMSCR:
        sepstr_scr_head( curr );
        break;
    }

    return curr;
}

    
/* return streaminf pointer given file desriptor . 0 on error */

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
streaminf * fd_info( int fd )
_XFUNCPROTOEND
#else
streaminf * fd_info(fd)
int fd;
#endif
{
    streaminf *curr;

    curr =sepstr_head();

    while( curr != 0 ){
   	if( fd == fileno(curr->streamfile)	) {
	    return curr;
	}
	curr = curr->next;
    }

    /* fell through this is an error */
    seperr("fd_info(): no stream has this file descriptor-> %d \n",fd );
    return 0;
}
