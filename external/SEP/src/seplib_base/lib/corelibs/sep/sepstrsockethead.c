#include <sepConfig.h>
#ifdef RS6000
#undef __STR__
#endif

#include <stdio.h>
#include <string.h>

#ifdef HAVE_ERRNO_H
#include <errno.h>
#endif

#ifndef STDC_HEADERS
extern int errno;
#endif

#if defined(__APPLE__) || defined(LINUX)
#define USE_SOCKETS
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif

#include <assert.h>

#include <sep_main_external.h>
#include "sep_main_internal.h"
#include "streamlist.h"

/* 
  sepstr_socket_head

	USAGE
  sepstr_socket_head( streaminf *info )     

	DESCRIPTION
	Open a socket for type STREAMSOCKOUT. Used (to send hff and gff down 
  a pipe). Used instead of traditional socket approach because socket
  number is also passed down the history file.  Therefore syncing operation
  is done AFTER hclose (and not by this routine.)

	COMMENTS
  info is a sepstream info structure with the tag field , entrytype,
  and the headername filled in.
  

*/
/* 
Author:  Bob 12-97 
Modified: Bob 6-99 GNU typedef

*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void sepstr_socket_head( streaminf *info )
_XFUNCPROTOEND 
#else
void sepstr_socket_head( info )
     streaminf *info;
#endif
{
  char portstr[256];
	char temp_ch[256],hostnm[255];
	int sock;


    assert( info->entrytype == STREAMSOCKOUT ); 

	info->format_num = FMT_XDR_FLOAT; /* default output format */

	/* make a hetch queue to store saved parameters in */
	info->hetch_queue = new_queue( 251 );
	info->hqlen = 251;


	portstr[0] = '\0';
  info->sockfd=opensock1(portstr, 0);
  if( gethostname(hostnm,255)!=0 ) seperr("sepstr_socket_head(): getting hostname\n");

	free( info->headername);
  sprintf(temp_ch,"%s:%s\n",hostnm,portstr);
	info->headername = (char *) malloc((int)strlen(temp_ch));
	strcpy(info->headername,temp_ch);

	
   /* open the input/output stream */
   open_socketstream( info );


}
/*  $Id: sepstrsockethead.c,v 1.2 2004/04/08 22:32:27 bob Exp $ */
