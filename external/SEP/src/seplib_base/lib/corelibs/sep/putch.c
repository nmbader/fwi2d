/*$

=head1 NAME

putch - put argument in output history file

=head1 SYNOPSIS

  int putch(name, format, val)

=head1 INPUT PARAMETERS

=over 4

=item name  -char*

      parameter name

=item format - char

      format of variable (d,f,s,l)

=item val- void*

      variable value of variable

=back


=head1 RETURN VALUE

=over 4

0 = if sucessful

=back


=head1 DESCRIPTION

This function is a seplib companion to the input routines fetch,
getch, and hetch that appends values to the output history file with
expressions of the form `name=value'.  



=head1 COMMENTS

Possible format conversions are
"d" (or "i"), "f" (or "r"),"m", "g" and "s" for integer, floating point,long long,
double precision, and strings respectively.

These values may be diverted from the default history file on standard
output with the command line parameter `head=filename'.

=head1 EXAMPLES

From C:

putch ( "intpar", "i", &int);

putch ( "intpar2", "d", &int);

putch ( "float", "f", &float);

putch ( "floatpar", "r", &float);

putch ( "somename", "s", "name");



From Fortran:

call putch ( 'intpar', 'i', int)

call putch ( 'intpar2', 'd', int)

call putch ( 'floatpar', 'f', float)

call putch ( 'float', 'r', float)

call putch ( 'somename', 's', 'name')


 
=head1 SEE ALSO

L<getch>, L<fetch>, L<hetch>, L<puthead>

=head1 DIAGNOSTICS

Program execution is terminated when an invalid format is passed.

=head1 BUGS

It is not presently possible to specify multiple output names in one
call as can presently be done on input with fetch, etc.

=head1 KEYWORDS 

parameter output 

=head1 LIBRARY

B<sep>

=cut
*/

 /*
 * Modified 2/8/83  S. Levin   line buffered output, write to head= file
 *			       instead of stdout if found. Wrote head()
 *			       function to return a stream pointer to
 *			       the header output stream. Changed putlin,
 *			       puthead, etc. to call head().
 *
 * Modified 3/4/83  S. Levin   added putch() as an alias for putch_    
 * Modified 7/14/83 S. Levin   added seperr () and perror() diagnostics
 * Modified 12/28/85 S. Levin  separated formatting from I/O for general use
 * Modified 9/7/87  S. Levin   <varargs> for portability
 * Revised: dave 9/17/90  Use stdarg for ANSI-C compilers
 * Revised: dave 8/94  Use new seplib internals
 */
#include <sepConfig.h>
#include    <stdio.h>
#include <sep_pars_external.h>
#include "streamlist.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int putch ( char *name, char *type, void* val )
_XFUNCPROTOEND
#else
int putch ( name, type, val )
char *name, *type;
char* val ;
#endif
{
static streaminf* outinfo = 0;

if( outinfo == 0 ) outinfo = tag_info( "out", TAG_OUT );

return sepstrput( outinfo, name, type, val );
}
