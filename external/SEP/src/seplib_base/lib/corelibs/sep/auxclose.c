/*$

=head1 NAME

auxclose - Close a SEPlib history file
 
=head1 SYNOPSIS

void auxclose(tag)

=head1 INPUT PARAMETERS 

=over 4 

=item	char* - tag  

      name of history file

=back

=head1 DESCRIPTION

	Auxclose closes, if necessary, an auxiliary input 
	or output previously opened with auxin or auxout.

=head1 COMMENTS

	Parameter `tag' is the name by which the auxiliary
	file is known.  All internal buffers associated 
	with that file will be freed as well.

=head1 LIBRARY

B<sep>

=cut

>*/

/*

SEE ALSO
	L<auxin>, L<auxout>, L<auxinout>, L<auxpar>, seplib

KEYWORDS   file close auxillary

WHERE
	./cube/libcube/auxclose.c

*/

/*
 * Revised: 7-22-87 Stewart A. Levin  force truncation of output
 * Revised: 9-15-90 Dave Nichols, tidy up for ANSI compilers.
 * Revised: 9/94 Dave Nichols,  rewrite using new seplib IO architecture.
 * Revised: 10/21/94 Dave Nichols, check for valid close fn. before using it.
 * Revised:  7/19/97 Robert Clapp, prototyping, returns
 */
#include <stdio.h>
#include "streamlist.h"
#include "sep_main_internal.h"
#include <sep_main_external.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int auxclose( const char *name)
_XFUNCPROTOEND 
#else
int auxclose(name)
char *name;
#endif
{
    streaminf *info;

    if( (info =tag_info(name,TAG_INQUIRE)) != 0 ){

     if( info->headerbuf != 0 ) free(info->headerbuf);
     if( info->entrytype == STREAMIN ){
	if( info->headfile != 0 ) fclose(info->headfile);
     }else{
	if( info->headfile != 0 ) { 
	    fflush( info->headfile);fclose(info->headfile);
        }
     }
     if( info->close_func != 0 ) (*info->close_func)(info, info->ioinf );
     if( info->headername != 0 ) free( info->headername);
     if( info->dataname != 0 ) free( info->dataname);

     sepstr_del( info );
     free(info);
   }
	return(0);
}
