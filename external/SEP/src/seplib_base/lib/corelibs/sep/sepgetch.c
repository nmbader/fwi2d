/*$

=head1 NAME

getch - grab a parameter from command line

=head1 SYNOPSIS

C<int getch(name, format, value)>


=head1 INPUT PARAMETERS

=over 4

=item name-  char* 

      parameter name

=item format -char* 

      format of variable (d,f,s,m)

=back

=head1 OUTPUT PARAMETERS

=over 4

=item value - void*

       variable value of variable

=back

=head1 RETURN VALUE

=over 4

=item x = number of matches found

=back


=head1 DESCRIPTION

This function is the seplib equivalent of the older getpar() for
extracting values from command line expressions of the form
`name=value'.  


=head1 COMMENTS

More than one name may be searched by separating the
alternatives by blanks or commas.  Possible format conversions are "d"
(or "i"), "f" (or "r"), "g" and "s" for integer, floating point,
double precision, and strings respectively.  If no match is found, the
variable is not altered.

When a keyword `par=filename' is encountered on the command line, that
file (and any that it might in turn point to) is scanned as well.
Getch returns a count of the number of matches it found.

 
=head1 DIAGNOSTICS

Program execution is terminated with extreme prejudice if a par file
cannot be opened or the command argument list cannot be found.
It is also an error for parfiles to be nested more then 10 deep.

=head1 BUGS

It should be possible to limit the length of retrieved strings and
to extract individual characters.  It is also not presently possible
to redirect search to an internal argument list.

=head1 KEYWORDS 

parameter command line argument

=head1 SEE ALSO

L<fetch>, L<hetch>, L<putch>, L<auxpar>, L<getch_add_string>

=head1 LIBRARY

B<sep>

=cut

*/

/*
 *
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Modified 8/4/89  Dave Nichols 
 *	   made the routine varargs, copied form Stew's stuff.
 * Revised: dave 9/17/90  Use stdarg for ANSI-C compilers
 * Revised: dave    1994  Make a copy of the buffer getch_add_string.
 * Revided: bob    7/97   Added prototypes, include lib_internal.h
 */

#include <sepConfig.h>
#include <stdlib.h>
#include "fastpar.h"
#include <string.h>
#include <sep_pars_external.h>

#define GETCH_QUEUE_SIZE 1023

hash_item *getch_queue[GETCH_QUEUE_SIZE];
int getch_queue_size = GETCH_QUEUE_SIZE;

static int first_invoke = 1;
extern int sepxargc; extern char **sepxargv;

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int getch( const char *tag, char *type, void *ptr )
_XFUNCPROTOEND
#else
int
getch(tag,type,ptr)
char *tag, *type;
char* ptr;
#endif
{
 register int iargc;
 MIXED var;


 if(first_invoke) {
    for(iargc=1; iargc<sepxargc; iargc++) {
	getpar_string_store(getch_queue,getch_queue_size,sepxargv[iargc]);
  }
    first_invoke = 0;
    }

 switch(type[0]) {
 case 'g': var.g = (double *)ptr; break;
 case 's': var.s = (char *)ptr; break;
 case 'f': case 'r': var.f = (float *)ptr; break;
 case 'm': var.m = (long long *)ptr; break;
 default: var.i = (int *)ptr; break;
 }
	
 return (getpar_decode(getch_queue,getch_queue_size,tag,type,var));
}
/*$

=head1 NAME

getch_add_string - add parameters to the command line

=head1 SYNOPSIS

getch_add_string(par)


=head1 INPUT PARAMETERS

=over 4

=item char* - par   

      string to add to the command line

=back


=head1 DESCRIPTION

Add a string to the database of command line arguments.  Useful
if you want to tell SEPlib behave in a certain manner (e.g. don't
output a history file.)

getpar, seplib, fetch, hetch, putch, auxpar getch_add_string

=head1 LIBRARY

B<sep>

=cut

*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int getch_add_string(char *string)
_XFUNCPROTOEND
#else
int getch_add_string(string)
char *string;
#endif
{
     char* copy;
     int iargc;

    if(first_invoke) {
        for(iargc=1; iargc<sepxargc; iargc++) 
    	    getpar_string_store(getch_queue,getch_queue_size,sepxargv[iargc]);
        first_invoke = 0;
    }
    else sepxargc++;

     copy = (char*)malloc( (int)strlen( string ) + 1 );
     strcpy( copy, string );
     getpar_string_store(getch_queue,getch_queue_size,copy);
		return(0);
}
