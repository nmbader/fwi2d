#include <sepConfig.h>
#ifdef RS6000
#undef __STR__
#endif

#include <stdio.h>
#include <string.h>

#if defined(HAVE_ERRNO_H) || defined(__APPLE__)
#include <errno.h>
#else

#ifndef STDC_HEADERS
extern int errno;
#endif
#endif


#include <sys/types.h>
#include <sys/stat.h>

#include <assert.h>
#include <sep_main_external.h>
#include "streamlist.h"

/* 

  sepstr_inout_head( info )     
  
  info is a sepstream info structure with the tag field , entrytype,
  and the headername filled in.
  
  This routine will open the stream for input/output. 
  This can only be done if both header and data are (separate) ordinary files.
  
  There are four possibilites for the headername specification,
  1) If it contains a "|" (pipe symbol) it is a command to use with popen.
  2) If it contains a colon it is a hostname:portnumber pair for opensock.
  3) It is "stdin".
  4) If none of the above are true it is a filename to be opened. 

*/

/*
  Author:
		Dave???
  
  Modified:
   Bob  7-97 Added different style prototypes
   Bob  6-99 GNU style ifdef
*/

char parambuf3[1024];

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void sepstr_inout_head( streaminf *info )
_XFUNCPROTOEND 
#else
void sepstr_inout_head( info )
     streaminf *info;
#endif
{
    struct stat statbuf;

    assert( info->entrytype == STREAMINOUT ); 
    assert( info->headername != 0 );
    
    /* input header */
    if( strchr( info->headername, '|' ) != 0 ||   /*pipe */
	strchr( info->headername, ':' ) != 0 ||     /* socket  */
    	strcmp( info->headername, "stdin" ) == 0  ||    /* stdin */
    	strcmp( info->headername, "stdout" ) == 0   ){   /* stdout */
      seperr("Headername \"%s\" for tag \"%s\" cannot be opened as an in/out dataset because it is not a regular file \n",info->headername, info->tagname);
    }

    /* must be a file */

    /* stat the header file name to see if it already exists 
	and isn't zero length */
    if(-1 != stat(info->headername,&statbuf) && ((int)statbuf.st_size) != 0 ) {

        /* header exists  open for read/write  */
        info->headfile = fopen(info->headername,"r+");

        if( info->headfile == 0 ) {
	    seperr("Headername \"%s\" for tag \"%s\" cannot be opened read/write \n",info->headername, info->tagname);
	}

        /* make sure the input stream is fully buffered */
        setvbuf(info->headfile,0,_IOFBF,BUFSIZ); 

        /* read the header into the buffer for this stream */
        readhdr( info );
	fseek(info->headfile,0,2); /* make sure we are at the end of the file */

	/* figure out the existing data format */
        if( sepstrpar(info, "data_format","s,",parambuf3) ){
           info->format_num = get_format_num(parambuf3);
	}else{
	   info->format_num = FMT_XDR_FLOAT; /* default output format */
	}

    }else{

        /* header does not exist */
        if(errno != ENOENT && ((int)statbuf.st_size) != 0) {  /* invalid name */
            perror("sepstrinouthead, openfile()");
            seperr("Bad header name '%s' for tag %s \n",
                   info->headername,info->tagname);
        }

        /* the name is OK but it doesn't exist, create it for writing/read */
        info->headfile = fopen(info->headername,"w+");

	info->format_num = FMT_XDR_FLOAT; /* default output format */

	/* make a hetch queue to store saved parameters in */
	info->hetch_queue = new_queue( 251 );
	info->hqlen = 251;
    }

    if( info-> headfile == 0 ){
	info->valid = 0;
	return;
    }
    
    /* open the input/output stream */
    open_inoutstream( info );

}
