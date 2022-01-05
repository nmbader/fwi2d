/*<
sep_get_grid_format_tag 

   
USAGE 
int sep_get_grid_format_tag(tag_history,tag_grid)

INPUT PARAMETER
   char-*tag_history    tag of History File

OUTPUT PARAMETER
   char-**tag_grid     tag to Header Format File

RETURN VALUE
   -1= if tag_history is not a SEP History File
    0= if successful
   +2= if tag_history has no Grid Format File associated with it.

DESCRIPTION
   Get tag for Header Format file pointed by tag_history.

CATEGORY
Lib:Sep3d:Grid access

COMPILE LEVEL
DISTR
>*/
/*
KEYWORDS
   grid 

SEE ALSO
   sep3d

AUTHOR
   Biondo Biondi , July ... 1995

*/

/* 
*
* LOGIC: 
* Determine if history tag points to an SEP history file
* check if grid tag is already specified in sepinfo
* case +: write tag in sepinfo to tag_grid
* case -:
*        if file is used for IN or INOUT 
*            check for grid_tag in history file
*                  case +: copy it 
*                  case -: is a SEP77 History file
*        else
*            check if grid tag is written on command line
*            if not create a default grid tag 
*                  tag_grid = tag_history + "^^"
*
*/

#include <sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include "sep3d.h"
#include "streamlist.h"
#include "sep_main_internal.h"
#include<time.h>


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_get_grid_format_tag(const char* tag_history, char** tag_grid)
_XFUNCPROTOEND
#else
int sep_get_grid_format_tag(tag_history, tag_grid)
char *tag_history; 
char **tag_grid;
#endif
{
    streaminf *info;
    char tag_grid_buf[1024],scratch_file[1024];
	  int test=0,put_end;
    int found;
     time_t tempt;

		

		/*changed so the first call can be SEP3d, might not be correct ...*/
    info = tag_info(tag_history, TAG_INQUIRE);  /* get info on this tag */
    if(info == SEPPOINTNULL) return 1;   /* Not a an SEP History File */




 if((info->gridformatfile) == SEPSTRNULL) {
     found=0;
     if((info->entrytype == STREAMINOUT)){
       if(auxpar("gff","s", tag_grid_buf, tag_history) == 1){
          found=1;
          *tag_grid=expandnm(tag_grid_buf,(char*)NULL);
           /*
          *tag_grid= (char *) malloc(1+(int)strlen(tag_grid_buf));
          strcpy(*tag_grid, tag_grid_buf);
          */
         info->gridformatfile = (char *)malloc(1+(int)strlen(*tag_grid));
         strcpy(info->gridformatfile, *tag_grid);
       } 
     }
     if((info->entrytype == STREAMIN) && found==0){ /*Input file */
       if(auxpar("gff","s", tag_grid_buf, tag_history) == 1) test=1;
       if(getch("gff_in","s", tag_grid_buf) == 1) test=1;
       if(test == 0 || 0==strncmp("-1",tag_grid_buf, 2)){
          return 2;     /* is a simple SEP77 history file */
				}
       else { /* tag is not in sepinfo but it is in history file,
            write it into sepinfo.                       */
  			/*This is a new error condition to force people to call init_3d at
				   the begining of any sep3d program */
				  if(0==strcmp("in",tag_history))
			    seperr("sep_get_grid_tag():must call init_3d at the begining of your sep3d program \n");
          *tag_grid=expandnm(tag_grid_buf,(char*)NULL);

         /*
         *tag_grid= (char *) malloc(1+(int)strlen(*tag_grid));
         strcpy(*tag_grid, tag_grid_buf);
         */
         info->gridformatfile = (char *)malloc(1+(int)strlen(*tag_grid));
         strcpy(info->gridformatfile, *tag_grid);
         /* open GFF appropriately */
         auxin(*tag_grid);
         return 0;
       }
      }
      else if(found==0){
        /*check for grid_format_file being overwritten on command line */
        if(!getch("gff","s", tag_grid_buf)) {
          /* for the special case of stdout info->gridname = stdout
          so we must getch the name from the command line */
          if(strcmp(info->tagname, "out") == 0){
            if ( info->isapipe || isapipe(fileno(info->headfile)) ){

              /*pt1 of 2 of lame trick.  Basically we are going to create
              a temporary gff file which we unlink from.  When sep3d_close
              is called the contents of the file are copied down a socket.
              Far from elegant, but it should work. Logic is otherwise
              we will have to create temp buffers for the gff file because
              we can't open socket until the history file is sent with
              the correct socket # info */

              /*first create the temporary file name */
              datapath(scratch_file);
              time(&tempt);
              srand(tempt);
              put_end=rand();
              sprintf(tag_grid_buf,"%s%s%d",scratch_file,"gff_",put_end);
            }

            else if(findnm(fileno(info->headfile), tag_grid_buf, sizeof(tag_grid_buf)) ==0){
              seperr("When stdout is not in the same directory you need to specify the Grid Format file with gff on the command line") ;
            }
          }
         /* Doesn't account for auxilary stuff directed along pipes  */
          else  strcpy(tag_grid_buf, info->headername);
          strcat(tag_grid_buf,"@@@@");
        }
        *tag_grid= (char *) malloc(1+(int)strlen(tag_grid_buf));
        strcpy(*tag_grid, tag_grid_buf);
        sep_put_grid_format_tag(tag_history, tag_grid_buf);
        return 0;
     }
    }
    else { /* sepinfo contains tag_grid, so just copy it */
      *tag_grid= (char *) malloc(1+(int)strlen(info->gridformatfile));
      strcpy(*tag_grid,info->gridformatfile);
       if(0==strncmp("-1",*tag_grid, 2)){ return 2;}
      	else{ return 0;}
    }

  return(0);
}


/*<
fget_grid_format_tag

USAGE
 int fget_grid_format_tag (tag_history,tag_grid)

INPUT PARAMETERS
  tag_history - char*   tag of history file

 OUTPUT PARAMETERS
  tag_grid - char*   tag of grid file


DESCRIPTION
  Allows fortran to access tag of grid file (could be written in a better way
  that is totally transparent, but don't want to bother for now)

CATEGORY
Lib:Sep3d:Grid access

COMPILE LEVEL
DISTR
>*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int fget_grid_format_tag(char * tag_history, char *tag_grid)
_XFUNCPROTOEND
#else
int fget_grid_format_tag(tag_history, tag_grid)
char * tag_history;
char *tag_grid;
#endif
{
char *temp[1];
int ierr;

ierr=sep_get_grid_format_tag(tag_history,temp);
if(ierr==0) {
  strcpy(tag_grid,*temp);
  free(*temp);
}

return ierr;
}




