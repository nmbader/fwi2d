/*$

=head1  NAME

hclose - close the seplib output history file

=head1 SYNOPSIS

	void hclose ()

=head1 DESCRIPTION 

Closes the stdout SEPLIB history file

=head1 COMMENTS

The end of a history file is delimited by a ctrl D (octal 004) so that
seismic data can be placed in back of it e.g. when passing data down a
pipe.  Hclose does this and must be called after all additions to the
output header and before writing output data.

=head1 DIAGNOSTICS

One of the functions of the include file <sep.wrapup> is to check
that the header has been closed before the program finishes.

=head1 BUGS
The seplib I/O routines, such as rite, srite, do not check to make
sure hclose has been called before writing to output.

=head1 KEYWORDS 

seplib close

=head1  LIBRARY

B<sep>


=cut 

*/


/*
 * hclose() function to clean up and close output header. This will 
 * automatically append a ^D and flush to permit next process in a
 * pipe to start running.  This should be called BEFORE writing data
 * to output file and AFTER all relevant parameters have been appended
 * to the output header.
 *
 * Author:  S. Levin   2/14/83
 * Revised: S. Levin   2/18/83    Added checking against header update after
 *				  hclose() call.
 * Revised: S. Levin   4/23/84    Write ^D trailer only if data to follow
 *				  This eases manual editing of headers.
 * Revised: S. Levin   1/31/86    New pipe synchronization for fast hetch
 *				  which avoids one char reads to avoid
 *				  passing ^D:  write pid as 9 ascii digit
 *				  number preceding ^D.  Sleep for a few
 *				  seconds waiting for SIGALARM signal
 *				  sent by other end of pipe.  hetch on
 *				  other end of pipe will send back the
 *				  signal to signify it's read the ^D.  If
 *				  no synch signal arrives the sleep expires
 *				  and processing continues assuming the
 *				  other end is older style and won't skip
 *				  ^D or is lucky.
 * Revised: S. Levin   1/25/89    Downgrade multiple hclose() calls to warning.
 * Revised D.Nichols   7/17/90    New pipe synchronization that works across
 *                                hosts. Write hostname, port number, & pid down
 *                                pipe and waits for a connection or 5 secs
 *				  or a sigalrm for an old version program.
 * Revised D.Nichols   9/11/92    Add data format tooheader before closing.
 *                                
 *                                
 */

#include <sepConfig.h>
#include <stdio.h>
#include "sep_main_internal.h"
#include <sep_main_external.h>
#define EOT 004
#define EOP 006

#include <assert.h>

#include "streamlist.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int hclose(void)
_XFUNCPROTOEND 
#else
int hclose()
#endif
{

static streaminf *outinfo;

if( outinfo == 0 )  outinfo = tag_info( "out", TAG_OUT );
if( outinfo == 0 ) return 0;

sepstr_hclose( outinfo );

return 0;

}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int auxhclose(const char *tag)
_XFUNCPROTOEND
#else
int auxhclose(tag)
const char * tag;
#endif
{

streaminf *info;

assert( tag != 0 );

info = tag_info( tag, TAG_INQUIRE );

if( info == 0 ) return 0;

sepstr_hclose( info );

return 0;
}
