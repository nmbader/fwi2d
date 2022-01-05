/*$

=head1 NAME

hetch - grab parameter from input history file

=head1 SYNOPSIS

  int hetch(name, format, variable)

=head1 INPUT PARAMETERS

=over 4

=item name -   char* 

      parameter name

=item format -  char* 

      format of variable (d,f,s)

=back 

=head1 OUTPUT PARAMETERS

=over 4

=item variable  -void * 

      value of variable

=back


=head1 RETURN VALUE

=over 4

=item x = number of matches found

=back


=head1 DESCRIPTION

This function is a seplib extension of the older getpar() that
extracts values from input history file expressions of the form
`name=value'.  

=head1 COMMENTS

More than one name may be searched by separating the
alternatives by blanks or commas.  Possible format conversions are "d"
(or "i"), "f" (or "r"), "g" and "s" for integer, floating point,
double precision, and strings respectively.  If no match is found, the
variable is not altered.

Hetch returns a count of the number of matches it found.  Also, the
input history file is copied to the output history file the first time
hetch (or fetch) is invoked.

Search and copy of the input history file can be supressed using the command
line parameter `noheader=y' in which case hetch returns 0.

=head1 DIAGNOSTICS

Program execution is terminated when
the command argument list cannot be found.

=head1 BUGS

It should be possible to limit the length of retrieved strings and
to extract individual characters.
The length of the input header is limited to 30,000 characters when
piped into a program.

=head1 KEYWORDS 

parameter input header

=head1 SEE ALSO

L<fetch>, L<getch>

=head1 LIBRARY

B<sep>

=cut

*/
/* 
 * hetch
 *
 * Revised  1-31-86  stew  new pipe syncronization:  look for 9 digit
 *                         pid preceding ^D and send SIGALRM to it to
 *			   let the other end of the pipe get on with
 *			   business as usual without further delay.
 * Revised  6-23-86  stew  always call puttitle()
 * Revised  7-21-86  kamal make it take care of tetch too
 * Revised  7-24-86  stew  corrected getpar interaction problems.
 * Revised  9-7-87   stew  <varargs> for portability
 * Revised: dave 9/17/90  Use stdarg for ANSI-C compilers
 * Revised  dave 7/17/90 Use internet socket to pass back a message saying
 *                       is OK to send data. Look for "PIPE" followed by
 *                       a hostname followed by a port number to connect to.
 *                       define USE_SOCKETS to get the new method.
 * Revised  dave 12/17/91 Added call to force input buffering for pipes.
 * Revised  dave     1994 Use new sepstrpar function.
 */
#include <sepConfig.h>
#include <stdio.h>

#include <sep_pars_external.h>
#include "streamlist.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int hetch( char *tag, char *type, void *val )
_XFUNCPROTOEND
#else
int hetch( tag, type, val )
char *tag, *type;
void *val ;
#endif
{
 static streaminf *ininfo = 0;

 if( ininfo == 0 )  ininfo= tag_info( "in", TAG_IN );

 return sepstrpar( ininfo, tag, type, val );

}
