/*
 * LOGIC:
 *         module for preparing sepstream data file.
 *         External entries defined here are:
 *
 *     open_instream( info )          initialize input data stream.
 *
 * entries in info filled in by this routine are:
 *          valid, dataname, format_num and all the IO function pointers
 *
 *
 * Author:  D. Nichols  8-94   SEP
 * Revised: S. Levin    5/4/95 Mobil/SEP  Avoid non-posix strdup
 * Revised  B. Clapp   6/2/99  Start GNU conversion
 *
 */
#include <sepConfig.h>
#include <stdio.h>
#include <fcntl.h>


#ifdef RS6000
#undef __STR__
#endif

#include <string.h>


#include <assert.h>

#include "streamlist.h"
#include "strformats.h"
#include <sepcube.h>
#include "sep_main_internal.h"

/* 

 To open an input we open the
 input file specified by "in=.." in the header. 
 If the value is "in=stdin" or "in=follow_hdr" 
 we read data from the same stream as the header.

 The "format=..." specifier is checked in the header to decide how to
 open the input. Currently the valid values are "native..." and "xdr.." 
 "vplot", and "cm_...."

*/

char parambuf[80024];

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void open_instream(streaminf* info)
_XFUNCPROTOEND
#else
void open_instream( info )
streaminf* info;
#endif
{
   int count;
   char *ptr;


    assert( info->headerbuf != (char*) 0 );
    assert( info->entrytype == STREAMIN );

    /* for stdin we can override the input with in= on the command line */
    count = 0;
    if( !strcmp(info->tagname,"in") ) count = getch( "in", "s", parambuf );

    /* otherwise we look for in= in the header */
    if( count==0  && sepstrpar( info, "in", "s", parambuf ) == 0  ){
	/* invalid input, user will get an error when they try to
         * perform I/O on it */
        info->format_num = FMT_NOT_KNOWN;
	info->valid=0;
	return;
    }

    if( strcmp( parambuf, "stdin" )==0 ){
	/* stdin means follow the header unless this is "in" */
	if( strcmp(info->tagname,"in")==0 ){
	    ptr = (char*)malloc((int)strlen("stdin")+1);
	    strcpy(ptr,"stdin");
	    info->dataname = ptr;
	}else{
	    ptr = (char*)malloc((int)strlen("follow_hdr")+1);
	    strcpy(ptr,"follow_hdr");
	    info->dataname = ptr;
	}
    }else{
	ptr = (char*)malloc((int)strlen(parambuf)+1);
	strcpy(ptr,parambuf);
        info->dataname = ptr;
    }

    /* initialise the pointers to I/O routines */
    init_io( info ); 

    /* check the input data format */
    if( sepstrpar( info, "data_format", "s", parambuf )  == 0 ){
	info->format_num = FMT_XDR_FLOAT; /* default format */
    }else{
	info->format_num = get_format_num( parambuf ); 
    }

    /* check if the data format is one of the known ones */
    if( info->format_num == FMT_NOT_KNOWN ){
	seperr("unknown data format: %s \n %s \n",parambuf,
	       "known types are xdr..., native.., cm..., and vplot");
    }

}
