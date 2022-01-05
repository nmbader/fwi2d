/*$

=head1 NAME

make_unpipe -unpipe a seplib file (therefore back seakable)

=head1 SYNOPSIS

make_unpipe(tag)

=head1 INPUT PARAMETERS

=over 4

=item tag - char* 

      tag to unpipe

=back


=head1 DESCRIPTION

Make sure that a seplib input dataset is not a pipe.
If it isn't a pipe we do nothing.
If it is we make a temporary dataset and then
copy the data to a new file which is made temporary by
unlinking it after it has been created. This ensures that it
will disappear when the job ends.


=head1 COMMENTS

WARNING DO NOT CLOSE THE FILE AFTER make_unpipe() has been called 

=head1 LIBRARY
B<sep>

=cut


>*/
/*
 * Modified  5/8/96  Stew Levin, expanded strdup() to malloc + strcpy
 *                   so LINUX gcc would shut up.
*/
#include <sepConfig.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "streamlist.h"
#include <assert.h>
#include "sep_main_internal.h"
#include <sep_main_external.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int make_unpipe( char *tag )
_XFUNCPROTOEND 
#else
int make_unpipe( tag )
char *tag;
#endif
{
char copybuf[BUFSIZ];
char scratch_file[4096];
streaminf *info;
FILE* tfile;
int num;

    assert( tag != 0 );

    if( (info = tag_info( tag, TAG_IN ) ) == 0 ){
        seperr("make_unpipe: unable to open tag \"%s\" \n", tag );
    }

    assert( info->entrytype == STREAMIN );

    /* leave it alone if it isn't a pipe */
    if( info->isapipe == 0 ) return (0);

    if( info->ioinf == 0 ) (*info->open_func)( info, &(info->ioinf) );
    if( !info->valid ){seperr("make_unpipe(): invalid input tag %s\n",tag);}

    /* make up a temporary dataset name */
    if(-1 == mkstemp(strcat(strtok(datapath(scratch_file),";"),"Unpipe_XXXXXX"))) {
      perror("make_unpipe(): mkstemp");
    }

    /* scratch_file now points to a null terminated scratch file name 
       guaranteed to be unique and identify whom it was made for */ 
    tfile = fopen(scratch_file,"w");
    if(  tfile == 0 ){
	perror("Opening scratch file in make_unpipe()");
	seperr("Unable to create scratch file \"%s\" to buffer a pipe",scratch_file);
    }


    /* copy the whole input to the file */
    while( ( num = (int)(info->read_func)(info,info->ioinf,copybuf,BUFSIZ )) != 0 ){
        if( num < 0 ){
	    perror("make_unpipe");
	    seperr("make unpipe read failed for tag \"%s\" \n",info->tagname);
        }
	if( fwrite( copybuf, 1, num, tfile ) != num ){
	    perror("make_unpipe");
	    seperr("make_unpipe: write failed for tag \"%s\" \n",info->tagname);
	}
    }
    /* close scratch file */
    fclose( tfile );

    /* close piped input file */
    (*info->close_func)(info, info->ioinf );
    free( info->dataname);
    info->ioinf = 0;
     
    /* setup info to read from scratch file now */
	  info->dataname = (char *) malloc( 1 + strlen(scratch_file)  );
     if(info->dataname != ((char *) NULL)) strcpy(info->dataname ,
                                               scratch_file  );
    init_fd_io( info );

    /* reopen temporary file */
    (*info->open_func)( info, &(info->ioinf) );

    /* Now we have our scratch file. */
    unlink(scratch_file);
    /* We have just turned this file into a ghost which is doomed
       to die completely when this program terminates */

    /* update the info structure */
    info->isapipe  = 0;
    info->backseekable = 1;
		return(0);
}
