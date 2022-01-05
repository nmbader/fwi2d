/*<

sep_tag_is_pipe

Usage:
ierr=sep_tag_is_pipe()

Return Values

-1= tag is a pipe, but not stdout
0 = not a pipe
1 = stdout is not a pipe



Input Parameters:


Output Parameters:



Description:
Check to see if stdout is a pipe


CATEGORY
Lib:Sep3d:Utility

COMPILE LEVEL
DISTR

>*/ 
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Sun Dec 14 13:24:39 PST 1997

Purpose: 

*/	 

#include <sep3d.h> 
#ifdef HAVE_UNISTD_H
#include<unistd.h>
#endif
#ifdef HAVE_STRING_H
#include<string.h>
#endif
#include "streamlist.h"
#include "sep_main_internal.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_tag_is_pipe(char *tag)
_XFUNCPROTOEND
#else
int sep_tag_is_pipe(tag)
char *tag;
#endif 

{ 
 streaminf *info;
 int haveheadfileno;

info = tag_info(tag, TAG_INQUIRE);  /* get info on this tag */
if(info == SEPPOINTNULL) return 1;   /* Not a an SEP History File */

if( info->headfile != SEPPOINTNULL ){
   haveheadfileno = isapipe(fileno(info->headfile));
} else {
  haveheadfileno = 0;
}


if ( info->isapipe || haveheadfileno ){
	if(strcmp(info->tagname, "out") == 0) return 1;
    /*data is going down a pipe (questionable if the stdout
     is a socket*/
	else return -1;
}
return  0;
} 


/*  $Id: tag_is_pipe.c,v 1.2 2004/04/08 22:32:27 bob Exp $ */
