/*$

putlin - put a line into history file

=head1 SYNOPSIS 
putlin(string)

=head1 INPUT PARAMETERS

=over 4

=item putlin - char* 

      line to put in history file

=back

=item DESCRIPTION

Put text line into history file 

=head1 SEE ALSO

L<putch>, L<puthead>

=head1 LIBRARY

B<sep>

=cut

>*/

/*
 * Modified 2/8/83  S. Levin   line buffered output, write to head= file
 *			       instead of stdout if found. Wrote head()
 *			       function to return a stream pointer to
 *			       the header output stream. Changed putlin,
 *			       puthead, etc. to call head().
 *
 * Modified 3/4/83  S. Levin   added putch() as an alias for putch_    
 *
 * Modified 7/14/83 S. Levin   added err() and perror() diagnostics
 */
#include <sepConfig.h>
#include    <stdio.h>
#include "sepstream.h"
#include <assert.h>

#include <sep_main_external.h>
#include <sep_pars_external.h>
#include "../include/streamlist.h"
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int putlin(char *line)
_XFUNCPROTOEND
#else
int putlin(line)
char *line;
#endif
{
	FILE *file; extern FILE *head();
	streaminf* info;

	assert( line != 0 );

	info = tag_info( "out", TAG_OUT );

	assert( info != 0 );

	if(  info->headfile == 0 ){
	     seperr("putlin() called after header is closed");
	};

	if( !info->header_title ) write_title( info );
	 	
	file = info->headfile;
	fprintf (file,"\t%s\n",line);
	if(ferror(file))
	   {
	    perror("putlin()");
	    seperr("putlin() I/O error on output header\n");
	   }
	fflush(file);
	return 0;
}
