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


#include <sys/types.h>
#include <sys/stat.h>

#include <assert.h>

#include <sep_main_external.h>
#include "sep_main_internal.h"
#include "streamlist.h"

/* 
  sepstr_scr_head

	USAGE
  sepstr_scr_head( streaminf *info )     

	DESCRIPTION

  
  This routine will open the stream for scratch usage (input/output)
  This can only be done if both header and data are (separate) ordinary files.

	COMMENTS
  info is a sepstream info structure with the tag field , entrytype,
  and the headername filled in.

  There are four possibilites for the headername specification,
  1) If it contains a "|" (pipe symbol) it is a command to use with popen.
  2) If it contains a colon it is a hostname:portnumber pair for opensock.
  3) It is "stdin".
  4) If none of the above are true it is a filename to be opened.

*/
/* 
Author:  Biondo??
Modified: Bob 7/97-Added different prototypes
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void sepstr_scr_head( streaminf *info )
_XFUNCPROTOEND 
#else
void sepstr_scr_head( info )
     streaminf *info;
#endif
{
  char portstr[256];
	char temp_ch[256],hostnm[255];
	int sock;

  struct stat statbuf;

    assert( info->entrytype == STREAMSCR );

    /* input header */
    if( strchr( info->headername, '|' ) != 0 ||   /*pipe */
  strchr( info->headername, ':' ) != 0 ||     /* socket  */
      strcmp( info->headername, "stdin" ) == 0  ||    /* stdin */
      strcmp( info->headername, "stdout" ) == 0   ){   /* stdout */
      seperr("Headername \"%s\" for tag \"%s\" cannot be opened as a scratch dataset because it is not a regular file \n",info->headername, info->tagname);
    }

    /* must be a file */

        /* the name is OK but it doesn't exist, create it for writing/read */
        info->headfile = fopen(info->headername,"w+");

  info->format_num = FMT_XDR_FLOAT; /* default output format */

  /* make a hetch queue to store saved parameters in */
  info->hetch_queue = new_queue( 251 );
  info->hqlen = 251;

    if( info-> headfile == (FILE*)0 ){
  info->valid = 0;
  return;
    }


    /* open the input/output stream */
    open_scrstream( info );


}
