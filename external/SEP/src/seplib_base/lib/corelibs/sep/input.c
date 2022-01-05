/*<
input

USAGE
FILE *input()-C ; int  input()

DESCRIPTION
Input returns a stream pointer (or a file descriptor for Fortran
users) that can be used to read input data into a program.  The
location of this data is determined according to the following
priorities:

COMMENTS

1) Look for in= on the command line.

2) Look for in= on the input history file.

3) Use standard input if it has been redirected.

The string found on the command line or in the history file may start
with an environment variable. e.g. in="$SCR/mydata.H@". If the environment
variable is set in the current environment then that value will be used.
If it is not set then the value found in the header file or on the command
line will be used. The "output()" routine writes the current value at the time
of creation in the output header.

C programmers rarely need to call these routines directly.
The include file <sep.startup> invokes them to define the
global variables `instream' and `infd'. Most people use the
routines "sreed" and "srite" to read an write data, these do not
require a file pointer as an input parameter. 

SEE ALSO
output, seplib, sreed, srite

DIAGNOSTICS
Input terminates program execution with an appropriate error message
if the input file could not be found or opened.

BUGS
It is risky to mix input stream I/O and file descriptor I/O
because of look ahead buffering.  <sep.startup> goes to some
pain to handle the case of input from stdin intelligently.

KEYWORDS: input stream stdin

CATEGORY
Lib:Sep:I/O

COMPILE LEVEL
DISTR
>*/

/* input -- fetch pointer to input file and open it. */
/* Processing history:
	?	jon	wrote the program
	2-3-83	ron	change 2nd fetpar from "out" to "out in"
	2-7-83	jon	change 2nd fetpar to hetch
	2-9-83  stew	removed fetpar call. clarified logic.
			return stdin fd when not otherwise defined.
			recognize out=stdout or in=stdin keywords
			as well to recognize data on pipe or old
			style programs with noheader=y and data on
			stdin.
	2-14-83	stew	evaluate once and save. Added input() to return
			stream pointer.
	2-17-83 stew	Force unbuffered input if from stdin.
	3-2-83	stew	Deleted diagnostic printout
	3-15-83 stew	Do not look for out= on input header!
	7-14-83 stew	Added perror() calls to improve diagnostics
	8-12-83 stew	Removed hidden 'i=s' in=stdin parameter. in=s still O.K.
	7-11-89 stew    Don't putch if hclose() has already been called.
	9-16-90  dave   made ansi-c and posix compatible 
 Revised: D. Nichols 7/29/91  Use filename expansion routine expandnm()
 Revissed: R. Clapp  7/18/97  Added include lib_internal
*/
#include <sepConfig.h>
#include <stdio.h>
#include <string.h>
#include "streamlist.h"
#include <sep_main_external.h>
#include "sep_main_internal.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
FILE *input(void)
_XFUNCPROTOEND
#else
FILE *input()
#endif
{
	static FILE *rfile = NULL ;
	streaminf* info;

   /* logic:
       get the stream info structure for tag "in" and return the
	streamfile member.
   */
   if(rfile == NULL || isclosed(rfile)) {
	info = tag_info( "in", TAG_IN );

	if( ! info->valid) return 0;

        if( info->ioinf == 0 ) (*info->open_func)(info, &(info->ioinf) );

	if( ! info->valid) return 0;

        rfile = info->streamfile ;
   }
   return (rfile);
}
