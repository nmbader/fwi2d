/*

NAME:
   sep_put_header_format_tag -- put ...

SYNOPSIS
   #include <sep3d.h>

   int sep_put_header_format_tag(char *tag_history, char *tag_header)

PARAMETER
   INPUT
   char *tag_history   :  tag of History File
   char *tag_header   :  tag to Header Format File

RETURN VALUE
    0 if successful
   -1 if fails for other reasons

DESCRIPTION
   Initialize the History File tag_history as a Sep3d History File.
   Initialize the Header Format File tag_header as a Sep3d Header 
   Format File. Link tag_header to tag_history.

KEYWORDS
   header 

SEE ALSO
   sep3d

AUTHOR
   Hector Urdaneta , July ... 1995

*/
#include <sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include <sep3d.h>
#include <unistd.h>
#include "streamlist.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_put_header_format_tag(const char* tag_history, char* tag_header)
_XFUNCPROTOEND
#else
int sep_put_header_format_tag(tag_history, tag_header)
char *tag_history; 
char *tag_header;
#endif
{
    streaminf *info;

    /* Initialize history file as Sep3d history file */
    info = tag_info(tag_history, TAG_INQUIRE);


    info->headerformatfile = (char *) malloc(1+(int)strlen(tag_header));
    strcpy(info->headerformatfile,tag_header);
    auxputch("hff","s", tag_header, tag_history);


    /* Open HFF file */
 
    switch( info->entrytype ){
      case( STREAMIN ):
  			seperr("Can not call sep_put_grid_format_tag for STREAMIN tag \n");
				break;
      case( STREAMOUT ):
/* Notice that is auxscr to allow routines to get information
on header keys also if it is a Header Format File linked with an output 
History File */
	      auxout(tag_header);
        if(1==sep_tag_is_pipe((char *) tag_history)){
          /*if we are dealing with stdout tag we want this file to
            disapear at the end of program execution */
        }
	      break;
      case( STREAMINOUT ):
	auxscr(tag_header);
	break;
      case( STREAMSOCKOUT ):
	auxscr(tag_header);
	break;
    }    
    return 0;
}
