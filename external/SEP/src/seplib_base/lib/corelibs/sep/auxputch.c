/*$
=head1 NAME

auxputch - put parameter into auxilary file

=head1 SYNOPSIS

C<int auxputch(name, type, val, tag)>

=head1 INPUT PARAMETERS

=over 4

=item name - 	char*

      name variable name

=item type- 	char* 

      format ((d our i),(f or r),(g or s),l,m)

=item	val - void* 

      variable value of variable

=item	tag - char* 

      tag of history file

=back

=head1 RETURN VALUE

0 = if successful


=head1 DESCRIPTION

Writes Parameters to an auxilary history file

=head1 COMMENTS

This function is a seplib companion to the input routine auxpar.
Auxputch writes or appends values to the auxiliary history file
Possible format conversions are "d" (or "i"), "f" (or "r"), "g" and "s"
for integer, floating point, double precision, and strings respectively.

The default auxiliary output header name given in string `tag'
for the auxiliary output header may be
overriden by coding `tag=header_file' on the command line.
(This is the same action taken by snap.)
The output header file is created if necessary.

Examples:
from C:
auxputch ( "int", "i", &int, "tag")
auxputch ( "int", "d", &int, "tag")
auxputch ( "float", "f", &float, "tag")
auxputch ( "float", "r", &float, "tag")
auxputch ( "name", "s", "name", "tag")

from Fortran:
call auxputch ( 'int', 'i', int, 'tag')
call auxputch ( 'int', 'd', int, 'tag')
call auxputch ( 'float', 'f', float, 'tag')
call auxputch ( 'float', 'r', float, 'tag')
call auxputch ( 'name', 's', 'name', 'tag')

=head1 SEE ALSO

L<auxpar>, L<slice>, L<auxclose>,L<putch>

=head1 DIAGNOSTICS

Program execution is terminated when an invalid format is passed.
Not all sequences of auxout) and auxputch have
been tested.

=head1 KEYWORDS 

auxillary header parameter output

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
 * Modified 3/4/83  S. Levin   added putch() as an alias for putch_    
 * Modified 7/14/83 S. Levin   added err() and perror() diagnostics
 * Modified 12/28/85 S. Levin  separated formatting from I/O for general use
 * Author  12/28/85  S. Levin  auxputch
 * Modified 9/7/87  S. Levin   <varargs> for portability
 * revised  9-16-90  dave      made ansi-c and posix compatible
 * Revised: 9/17/90  dave      Use stdarg for ANSI-C compilers
 *
 */
#include <sepConfig.h>
#include    <stdio.h>
#include    "streamlist.h"
#include <sep_pars_external.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int auxputch( const char *name, const char *type, const void* val, const char* tag )
_XFUNCPROTOEND
#else
int auxputch( name, type, val, tag )
const char *name;
const char  *type, *tag ;
const void * val;
#endif
{
streaminf *info;


info = tag_info( tag, TAG_OUT );


return( sepstrput( info, (char *) name, (char *) type, (void *) val ) );

}
