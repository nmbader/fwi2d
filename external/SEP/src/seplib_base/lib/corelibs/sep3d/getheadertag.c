/*<
sep_get_header_format_tag
USAGE
   int sep_get_header_format_tag(tag_history,tag_header)
INPUT PARAMETER
   char* - tag_history     tag of History File
OUTPUT PARAMETER
   char** - tag_header     tag to Header Format File

RETURN VALUE
   -1=if tag_history is not a SEP History File
    0=if successful
   +2=if tag_history is an Sep77 History File

DESCRIPTION
   Get tag for Header Format file pointed by tag_history.
CATEGORY
Lib:Sep3d:Header access

COMPILE LEVEL
DISTR
>*/
/*
KEYWORDS
   header 

SEE ALSO
   sep3d

AUTHORS
   Hector Urdaneta , July ... 1995
   Biondo Biondi , July ... 1995

*/

/* 
*
* LOGIC: 
* Determine if history tag points to an SEP history file
* check if header tag is already specified in sepinfo
* case +: write tag in sepinfo to tag_header
* case -:
*        if file is used for IN or INOUT 
*            check for header_tag in history file
*                  case +: copy it 
*                  case -: is a SEP77 History file
*        else
*            check if header tag is written on command line
*            if not create a default header tag 
*                  tag_header = tag_history + "^^"
*
*/

#include <sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include<time.h>
#include "sep3d.h"
#include "streamlist.h"
#include "sep_main_internal.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_get_header_format_tag(const char* tag_history, char** tag_header)
_XFUNCPROTOEND
#else
int sep_get_header_format_tag(tag_history, tag_header)
char *tag_history; 
char **tag_header;
#endif
{
    streaminf *info;
    char tag_header_buf[1024],scratch_file[1024],*junk_tag;
    int test=0,put_end,found;
     time_t tempt;

   info = tag_info(tag_history, TAG_INQUIRE);  /* get info on this tag */
   if(info == SEPPOINTNULL) return 1;   /* Not a an SEP History File */





   if((info->headerformatfile) == SEPSTRNULL) {
      found=0;
     if((info->entrytype == STREAMINOUT)){
       if(auxpar("hff","s", tag_header_buf, tag_history) == 1){
         found=1;
         *tag_header= (char *) malloc(1+(int)strlen(tag_header_buf));
		     strcpy(*tag_header, tag_header_buf);
		     info->headerformatfile = (char *)malloc(1+(int)strlen(tag_header_buf));
		     strcpy(info->headerformatfile, *tag_header);
        }
     }
	   if((info->entrytype == STREAMIN) && found==0){ /*Input file */
	     if(auxpar("hff","s", tag_header_buf, tag_history) == 1) test=1;
	     if(getch("hff_in","s", tag_header_buf) == 1) test=1; 
	     if(test == 0 || !strncmp("-1",tag_header_buf, 2)) 
		      return 2;     /* is a simple SEP77 history file */
	     else { /* tag is not in sepinfo but it is in history file,
		        write it into sepinfo.                       */
				/*This is a new error condition to force people to call init_3d at
    			the begining of any sep3d program */
   			if(0==strcmp("in",tag_history))
     			seperr("sep_get_header_tag():must call init_3d at the begining of your sep3d program \n");
         *tag_header=expandnm(tag_header_buf,(char*)NULL);
          /*
         *tag_header= (char *) malloc(1+(int)strlen(*tag_header));
		     strcpy(*tag_header, tag_header_buf);
          */
		     info->headerformatfile = (char *)malloc(1+(int)strlen(*tag_header));
		     strcpy(info->headerformatfile, *tag_header);
		     /* open HFF appropriately */
       	 auxin(*tag_header);
		     return 0;
     	 }
	    }    
      else if(found==0){
        /*check for header_format_file being overwritten on command line */
	      if(!getch("hff","s", tag_header_buf)) {
		      /* for the special case of stdout info->headername = stdout
		      so we must getch the name from the command line */
		      if(strcmp(info->tagname, "out") == 0){
            /*data is going down a pipe (questionable if the stdout
              is a socket*/
            if ( info->isapipe || isapipe(fileno(info->headfile)) ){

              /*pt1 of 2 of lame trick.  Basically we are going to create
              a temporary hff file which we unlink from.  When sep3d_close
              is called the contents of the file are copied down a socket.
              Far from elegant, but it should work. Logic is otherwise 
              we will have to create temp buffers for the hff file because
              we can't open socket until the history file is sent with
              the correct socket # info */ 
 
              /*first create the temporary file name */
              datapath(scratch_file);
              time(&tempt);
              srand(tempt);
              put_end=rand();
              sprintf(tag_header_buf,"%s%s%d",scratch_file,"hff_",put_end);
						}
            else if(findnm(fileno(info->headfile), tag_header_buf, 
              sizeof(tag_header_buf)) ==0){
			        seperr("When stdout is not in the same directory you need to specify the Headerr Format file with hff on the command line") ;
						}
		      }
          else  strcpy(tag_header_buf, info->headername);
		      strcat(tag_header_buf,"@@");
	      }
        *tag_header=expandnm(tag_header_buf,(char*)NULL);
         /*
        *tag_header= (char *) malloc(1+(int)strlen(*tag_header));
	      strcpy(*tag_header, tag_header_buf);
        */
	      sep_put_header_format_tag(tag_history, tag_header_buf);
	      return 0;
	   }
    }
    else { /* sepinfo contains tag_header, so just copy it */
    	*tag_header= (char *) malloc(1+(int)strlen(info->headerformatfile));
      strcpy(*tag_header,info->headerformatfile);
	    return 0;
    }

return(0);
}
/*<
fget_header_format_tag

USAGE
 int fget_header_format_tag (tag_history,tag_header)

INPUT PARAMETERS
  tag_history - char*   tag of history file

 OUTPUT PARAMETERS
  tag_header - char*   tag of header file


DESCRIPTION
  Allows fortran to access tag of header file (could be written in a better way
  that is totally transparent, but don't want to bother for now)

CATEGORY
Lib:Sep3d:Header access

COMPILE LEVEL
DISTR
>*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int fget_header_format_tag(char * tag_history, char *tag_header)
_XFUNCPROTOEND
#else
int fget_header_format_tag(tag_history, tag_header)
char * tag_history;
char *tag_header;
#endif
{
char *temp[1];
int ierr;

ierr=sep_get_header_format_tag(tag_history,temp);
if(ierr==0) {
  strcpy(tag_header,*temp);
  free(*temp);
}

return ierr;
}


