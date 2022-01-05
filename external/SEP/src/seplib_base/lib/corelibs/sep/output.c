/*<
output


USAGE
FILE *output() -C;  int output() -Fortran

DESCRIPTION
Output returns a stream pointer (or a file descriptor for Fortran
users) that can be used to write output data from a program.  The
location of this data is determined according to the following
priorities:

COMMENTS
1) Look for out= on the command line.

2) Use stdout if it is a pipe.

3) Construct an output file name by attaching the prefix returned from
 datapath() to the name of the program executing at the time.

If the final path starts with an environment variable the current value of
that variable will be written in the output header.

C programmers rarely need to call these routines directly.  The
include file <sep.startup> invokes them to define the global variables
`outstream' and `outfd'. Most prograns in C and on FOrtran use the
routines "sreed" and "srite" to read and write data. These programs do
not use file pointers explicitly but output() is called by the
routines to determine the output dataset location.

SEE ALSO
input, rite, srite, datapath, seplib

DIAGNOSTICS

Output terminates program execution with an appropriate error message
if the output file could not be opened.  

BUGS 
It is risky to mix output stream I/O and file descriptor I/O because
of stream buffering. flush should be called whenever you switch
out of stream mode.  When the output is a pipe, hclose (9s) must
be called before the first write to output.

KEYWORDS: write output outstream

CATEGORY
Lib:Sep:I/O

COMPILE LEVEL
DISTR
>*/
/* wfile = output()   -- Find a file descriptor for data output

	if (getch finds out=filename) then open it
	if (redirected output and oldstyle prog) use stdout
	else / * make up a good name * / use default datapath
	    / * and open that new file * /
            open a file for writing at datapath_ProcessName
		example
			makefile:   datapath=/scr/jon/D$(DATA)
			implies     out=/scr/jon/Dwz25window

	puthed out=outfile in=outfile

 Revised: S. Levin   2/11/83
	Setup scheme to default foot via .datapath files
	or DATAPATH environment variable.
 Revised: S. Levin   2/15/83
	Added output() to return stream pointer  check datapath for stdout
	default.  Fixed null user id bug when running in background.
 Revised: S. Levin   2/21/83
	Changed default to /scr/username/_ because # generally signals
	the start of a comment in a makefile.
 Revised: S. Levin   3/1/83
	Moved datapath code to separate subroutine datapath()
 Revised: S. Levin   3/2/83
	Deleted diagnostic printout
 Revised: S. Levin   3/15/83
	Use new isapipe subroutine to automatically send data down pipe
	unless specifically overridden
 Revised: S. Levin   3/26/83
	Call noheader() .
 Revised: S. Levin   7/14/83
	Added perror() calls to improve diagnostics.
 Revised: S. Levin   8/12/83
	Removed hidden 'o=s' out=stdout abbrev.  out=s still works.
 Revised: S. Levin   4/10/84
	Removed unnecessary putch of out=
 Revised: S. Levin   4/29/84
	Added dummy use of hetch() to force copying input to output header
 Revised: S. Levin   5/15/84
	Changed default output datafile naming procedure.  Use findnm()
	to look for output header name in current directory and construct
	datapath suffix of the form {header name minus intial H or final
	.[CAPITAL LETTER(S)] tail}.{Program name}   Otherwise just use
	original program name scheme.
 Revised: S. Levin   5/18/84
	Changed final .{program name} simply to .DATA - in most cases
	the program name is already embedded in the header name (which
	is what we really wanted to know anyhow).
 Revised: Joe Dellinger 7/26/84
	noheader=y only implies out=stdout if you have also redirected the
	header away from stdout.
 Revised: S. Levin   2/7/85
	Changed final .DATA - to suffix of output header+# char
	This avoids clobbering of input for <a.H prog >a.M
 Revised: S. Levin   10/1/85
	Changed final # to @.  Same reasons as earlier change
 Revised: D. Nichols 9/16/90  update to posix & ansi-c 
 Revised: D. Nichols 7/29/91  Use filename expansion routine expandnm()
 Revised: R. clapp  7/18/97   Added lib_internal include
 */
#include <sepConfig.h>
#include <stdio.h>
#include <sep_main_external.h>
#include "sep_main_internal.h"
#include "streamlist.h"


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
FILE * output (void)
_XFUNCPROTOEND
#else
FILE * output ()
#endif
{
    static  FILE * wfile = NULL;/* will hold output file descriptor */
    streaminf *info;

    if (wfile == NULL || isclosed (wfile))/* need to find an output */
    {
	info = tag_info( "out", TAG_OUT );
	if( info->ioinf == 0 ) (*info->open_func)(info, &(info->ioinf) );

	if( !info->valid ){
	    wfile =  0 ;
	}else{
	   wfile = info->streamfile;
        }
    }
    return (wfile);
}
