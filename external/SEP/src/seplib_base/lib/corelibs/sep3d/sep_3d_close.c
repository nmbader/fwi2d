/*$

=head1 NAME

sep_3d_close - close sep3d format files


=head1 SYNOPSIS

sep_3d_close()


=head1 DESCRIPTION

Closes sep_3d program, allows piping.

=head1 SEE ALSO

L<init_3d>



=head1 LIBRARY

B<sep3d>

=cut


>*/ 
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Sun Dec 14 11:44:56 PST 1997

Purpose: 

*/	 

#include <sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include<stdio.h>
#include <sep3d.h>
#if defined(__APPLE__) || defined(LINUX)
#define USE_SOCKETS
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "streamlist.h"
#include "sep_main_internal.h"
#include <sep3d.h>
#define EOL 014



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void sep_3d_close(void) 
_XFUNCPROTOEND
#else
void sep_3d_close() 
#endif 
{ 
char *tag_header[1];
char *tag_grid[1];
char temp_ch[1024];
streaminf *info_history,*info_hff_in,*info_gff_in;
streaminf *info_hff_out,*info_gff_out;
int hff=0,gff=0,sock,timeout,ierr,next,whence;
off_t seekpos;




 /*CHECK TO SEE IF STDOUT IS BEING SENT DOWN PIPE */
if(1==sep_tag_is_pipe("out")){
	 info_history = tag_info("out", TAG_INQUIRE); 
   /*FIRST DEAL WITH HEADER */
	 if((info_history->headerformatfile) != SEPSTRNULL &&
    0==getch("hff","s",temp_ch)){
		 hff=1;
		/*WE ARE SENDING THE OUTPUT DOWN A PIPE AND DID NOT SPECIFY A HFF FILE
      ON THE COMMAND LINE. THIS MEANS WE CREATED A TEMP FILE IN 
      sep_get_header_format_tag. */
			if(0!=sep_get_header_format_tag("out",tag_header))
				seperr("sep_3d_close(): trouble getting tag_header \n");
			info_hff_in = tag_info(*tag_header, TAG_INQUIRE); 
    	if( info_hff_in->headfile == (FILE*)0 ){
				/*AUXCLOSE FOR THE HEADER FILE HAS ALREADY BEEN CALLED. THIS
          PROBABLY MEANS THAT SOMEONE HAS STARTED WRITING OUT HEADERS.
          LATER WE SHOULD SET THIS UP SO WE CAN COPY EVERYTHING THAT
          HAS ALREADY BEEN WRITTEN OUT THE SOCKET.  FOR NOW JUST GIVE
          AN ERROR */
        seperr("Can not call auxclose on the headers (most likely by writing out the headers) before sep_3d_close if you wish to pipe\n");
      }
      /*NOW WE NEED TO DO SOME MAGIC, FIRST CREATE*/
			auxsockout("hff_out");
			info_hff_out = tag_info("hff_out", TAG_INQUIRE); 
			putch("hff","s",info_hff_out->headername);
    }
   /*NOW THE GRID*/
	 if((info_history->gridformatfile) != SEPSTRNULL &&
    0==getch("gff","s",temp_ch)){
		 	gff=1;
		/*WE ARE SENDING THE OUTPUT DOWN A PIPE AND DID NOT SPECIFY A GFF FILE
      ON THE COMMAND LINE. THIS MEANS WE CREATED A TEMP FILE IN 
      sep_get_grid_format_tag. */
			if(0!=sep_get_grid_format_tag("out",tag_grid))
				seperr("sep_3d_close(): trouble getting tag_grid \n");
			info_gff_in = tag_info(*tag_grid, TAG_INQUIRE); 
    	if( info_gff_in->headfile == (FILE*)0 ){
				/*AUXCLOSE FOR THE GRID FILE HAS ALREADY BEEN CALLED. THIS
          PROBABLY MEANS THAT SOMEONE HAS STARTED WRITING OUT GRID.
          LATER WE SHOULD SET THIS UP SO WE CAN COPY EVERYTHING THAT
          HAS ALREADY BEEN WRITTEN OUT THE SOCKET.  FOR NOW JUST GIVE
          AN ERROR */
        seperr("Can not call auxclose on the grid (most likely by writing out the headers) before sep_3d_close if you wish to pipe\n");
      }
      /*NOW WE NEED TO DO SOME MAGIC, FIRST CREATE*/
			auxsockout("gff_out");
			info_gff_out = tag_info("gff_out", TAG_INQUIRE); 
			putch("gff","s",info_gff_out->headername);
    }
}


hclose();

/*WITH THAT HCLOSE WE SHOULD HAVE STARTED THE LINKING ROUTINE
  SO NOW WE NEED TO DO THE STANDARD SYNCHRONIZATION */
  timeout=5; getch("timeout","i",&timeout);

	if(hff==1){
		if( (info_hff_out->sockfd=socklisten(info_hff_out->sockfd, timeout ))==-1) 
      seperr("sep_3d_close(): pipe synch failed for hff\n");
	
			info_hff_out->headfile = fdopen( info_hff_out->sockfd, "w" );
			if( info_hff_out->valid && info_hff_out->ioinf == 0 ){
    	  (*info_hff_out->open_func)( info_hff_out, &(info_hff_out->ioinf) );
			}



		/*AT THIS POINT WITH ANY LUCK WE HAVE ESTABLISHED COMMUNICATION. THE
      NEXT TRICK IS TO COPY THE CONTENTS OF THE OLD HFF FILE TO THE NEW
      HFF FILE */
			seekpos=(off_t) 0;whence=0;
	 		if( 0!=fseek(info_hff_in->headfile,(long)seekpos,whence) ) perror("file_seek");


    	while( (next = getc(info_hff_in->headfile)) != EOT && next != EOF ){
        /* read until we get an EOT or EOF */
				ierr=putc(next,info_hff_out->headfile);
    }
		fflush(info_hff_out->headfile);

			auxhclose("hff_out");
			auxclose(*tag_header);
			unlink(*tag_header);
			auxin(*tag_header);
			info_hff_in = tag_info(*tag_header, TAG_INQUIRE); 
			free(*tag_header);
			free(info_history->headerformatfile);
			info_history->headerformatfile=(char *)malloc(1+(int)strlen("hff_out"));	
			strcpy(info_history->headerformatfile,"hff_out");
   }

	if(gff==1){
		if( (info_gff_out->sockfd=socklisten(info_gff_out->sockfd, timeout )) ==-1)
      seperr("sep_3d_close(): pipe synch failed for gff\n");

			info_gff_out->headfile = fdopen( info_gff_out->sockfd, "w" );
			if( info_gff_out->valid && info_gff_out->ioinf == 0 ){
    	  (*info_gff_out->open_func)( info_gff_out, &(info_gff_out->ioinf) );
			}
		/*AT THIS POINT WITH ANY LUCK WE HAVE ESTABLISHED COMMUNICATION. THE
      NEXT TRICK IS TO COPY THE CONTENTS OF THE OLD HFF FILE TO THE NEW
      HFF FILE */
			seekpos=(off_t) 0;whence=0;
	 		if( 0!=fseek(info_gff_in->headfile,(long)seekpos,whence) ) perror("file_seek");


    	while( (next = getc(info_gff_in->headfile)) != EOT && next != EOF ){
        /* read until we get an EOT or EOF */
				ierr=putc(next,info_gff_out->headfile);
    }

			auxhclose("gff_out");
			auxclose(*tag_grid);
			unlink(*tag_grid);
			auxin(*tag_grid);
			info_gff_in = tag_info(*tag_grid, TAG_INQUIRE); 

		/*AT THIS POINT WITH ANY LUCK WE HAVE ESTABLISHED COMMUNICATION. THE
      NEXT TRICK IS TO COPY THE CONTENTS OF THE OLD HFF FILE TO THE NEW
      HFF FILE */
			free(*tag_grid);
			free(info_history->gridformatfile);
			info_history->gridformatfile = (char *) malloc(1+(int)strlen("gff_out"));	
			strcpy(info_history->gridformatfile,"gff_out");
   }
	/*IN THEORY, AT LEAST, WE ARE DONE */


 

return  ;
} 


/*  $Id: sep_3d_close.c,v 1.2 2004/04/08 22:32:27 bob Exp $ */
