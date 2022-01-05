/* 
 * A routine to set the pointers to the appropriate I/O routines.
 * It checks the iotype files on the info structure and does the
 * "right thing". To add a new I/O type you must put a new case
 * in this routine and add a new type to the enumeration in 
 * sepstream.h
 *
 * Revised: Biondo 8/2/95 changed error condition for STREAMINOUT
 * Revised: Bob   7/20/97 changed to NeedFunctionProto
 * Revised: Bob   12/18/97 added STREAMSOCKOUT
 */
#include <sepConfig.h>
#include "streamlist.h"
#include "sep_main_internal.h"
#include <sep_main_external.h>

#include <string.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void init_io( streaminf* info )
_XFUNCPROTOEND 
#else
void init_io( info )
streaminf * info;
#endif
{


 
 if( !strcmp(info->dataname,"stdin") ||
      !strcmp(info->dataname,"follow_hdr" ) ||
      !strcmp(info->dataname,"stdout") ){

      /* use bufferd I/O if the data file is the same as the header
       * because we have already used buffered I/O to read/write the header */
      info->iotype = FILE_IO;

 }else if( strchr( info->dataname, ';' ) != 0 ) {
      info->iotype = MULTI_FD_IO ;

 }else{
        /* use raw I/O if the data file is separate from the header */
      info->iotype = FD_IO ;
 }
 if( (!strcmp(info->dataname,"follow_hdr")||!strcmp(info->dataname,info->headername)) && ((info->entrytype == STREAMINOUT) || (info->entrytype == STREAMSCR))){
     seperr("tag \"%s\" in/out dataset cannot be a pipe or the same file as the header\n",info->tagname);
 }


 switch( info->iotype ){
	case( FILE_IO ):
	    init_file_io( info );
	    break;
	case( FD_IO ):
	    init_fd_io( info );
	    break;
	case( MULTI_FD_IO ):
	    init_multifd_io( info );
	    break;
	default:
	    seperr("unknown IO type\n");
 }

}

