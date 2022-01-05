/*$
=head1 NAME

fetch - grab a parameter from the command line or history file

=head1 SYNOPSIS

C<int fetch(name, format, variable)>

=head1 INPUT PARAMETERS

=over 4

=item char* - name   

      parameter name

=item char* - format 

      format of variable (d,f,s)

=back

=head1 OUTPUT PARAMETERS

=over 4

=item variable - void*

       variable value of variable

=back

=head1 RETURN VALUE

=over 4

=item x = number of matches found

=back

=head1 DESCRIPTION

This function extracts values from input history file and command line
expressions of the form `name=value'.  


=head1 COMMENTS

More than one name may be searched by separating the alternatives by blanks 
or commas.  Possible format conversions are "d" (or "i"), "f" (or "r"), "g", 
"s" and "1" for integer, floating point, double precision, strings, 
and boolean respectively.  (For boolean, y, Y, or 1 are returned 
as integer 1 and n, N, or 0 are returned as integer 0.)  
If no match is found, the variable is not altered.

When a keyword `par=filename' is encountered on the command line, that
file (and any that it might in turn point to) is scanned as well.
Fetch returns a count of the number of matches it found.  It also
copies the input history file to the output history file if this has
not already been done.

Search and copy of the input history file can be suppressed using the command
line parameter `noheader=y'.

=head1 EXAMPLES 

Some examples of how this routine works:
Command			Input		Result

C<fetch("nx","d",&nx)>nx=100		nx is set to 100, 1 argument found

C<fetch("nx","d",&nx)>nx=		0 arguments found

C<fetch("nx nxx","d",&nx)>nxx=30		nx is set to 30, 1 argument found

C<fetch("nx","d",&nx)>nx=50 nx=60	nx is set to 60, 2 arguments found

C<fetch("nt","d",&nx)>	nx=80		0 arguments found

C<fetch("nt ntt","d", &nt)> ntt=30 nt=40	nt is set to 40, 2 arguments found

C<fetch("nt ntt","d", &nt)> nt=30 ntt=40	nt is set to 40, 2 arguments found


=head1 SEE ALSO

L<getch>, L<hetch>, L<putch>, L<auxpar>

=head1 DIAGNOSTICS

Program execution is terminated with extreme prejudice if a par file
cannot be opened or the command argument list cannot be found.

=head1 BUGS

It should be possible to limit the length of retrieved strings and
to extract individual characters.  It is also not presently possible
to redirect search to an internal argument list.

=head1 KEYWORDS 

parameter input fetch 

=head1  LIBRARY

B<sep>

=cut
*/

/* edit history
clayton		1981	wrote getpar
hale		1982	?  INPAR?
JFC		1-22-83	stdin copied to stdout and parsed first.
ron		1-28-83 allow for multiple names of a variable
ron		2-1-83  does not count entries where there is no
			argument as being found
ron		2-2-83  changed the stdin read from character by character
			to line by line
ron		2-4-83	altered the stdin read to not keep lines starting with '#'
JFC		2-7-83  hetch introduced
stew            2-8-83  changed fetch_'s copy operation to write to
			FILE returned by head()
stew		2-9-83  hetch looks for noheader=y keyword and returns if there.
			don't copy past CTL D in hetch.
ron		2-15-83 added count_found so that the routines in fetch and
			getpar would use the same found.
stew		2-17-83 forced hetch_ to unbuffer stdin before copying
ron		2-23-83 added getm_getsav so that fetch call will not affect
			value until the last parameter is found
stew		3-1-83  defined time() in terms of time_t typedef in sys/types.h
stew		3-26-83 call noheader()
stew		7-14-83 added err and perror calls for improved diagnostics
stew		8-8-83  split fortran versions to separate files for linkage
stew		1-8-86  compatibility with fast getpar
stew		9-6-87	<varargs> for portability
dave 		9/17/90  Use stdarg for ANSI-C compilers

	 copyright (c) 1991 Stanford University

*/

#include <sepConfig.h>

#include <sep_pars_external.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int fetch( char *name, char *type, void* val )
_XFUNCPROTOEND
#else
fetch(name,type,val)
char* name, *type;
char *val;
#endif
{
    int rc;


    switch(type[0]) {
    case 'g': rc=getch(name,type,val); break;
    case 's': rc=getch(name,type,val); break;
    case 'f': case 'r': rc=getch(name,type,val); break;
    default: rc=getch(name,type,val); break;
    }
    if(0 < rc) return(rc);

    switch(type[0]) {
    case 'g': rc=hetch(name,type,val); break;
    case 's': rc=hetch(name,type,val); break;
    case 'f': case 'r': rc=hetch(name,type,val); break;
    default: rc=hetch(name,type,val); break;
    }
    return(rc);
}
