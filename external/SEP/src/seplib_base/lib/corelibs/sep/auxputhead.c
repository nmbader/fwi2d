/*$
=head1 NAME

auxputhead - put a formated string to seplib history file

=head1 SYNOPSIS

C<auxputhead(tag,format, value1, value2, ...)>

=head1 INPUT PARAMETERS

=over 4

=item tag    - char*

      Tag of file tow rite out

=item format - char* 

      C format statement

=item value1..n  - 

      void  values

=back

=head1 DESCRIPTION

This function permits the C programmer to add information to the
output history file using printf (3) format control.  

=head1 COMMENTS
These values may
be diverted from the default history file on standard output with the
command line parameter `head=filename'.  Thus one might code

C<puthead("out","\\tn1=%d n2=%d esize=%d\\n",nsamples,ntraces,4);>

to reflect a change in size and format of the input data.

=head1 SEE ALSO

printf(3), seplib, L<putch>,L<puthead>

=head1 DIAGNOSTICS

See printf(3) for error handling.

=head1 KEYWORDS 

header print printf

=head1 LIBRARY

B<sep>

=cut

>*/
/*	
 *
 *	auxputhead(" printf format with %s, %d etc","VAX-11-",780);
 *
 * Modified 2/8/83  S. Levin   line buffered output, write to head= file
 *			       instead of stdout if found. Wrote head()
 *			       function to return a stream pointer to
 *			       the header output stream. Changed putlin,
 *			       puthead, etc. to call head().
 * Modified 3/4/83  S. Levin   added putch() as an alias for putch_    
 * Modified 7/14/83 S. Levin   added err() and perror() diagnostics
 * Modified 9/7/87  S. Levin   <varargs> for portability
 * Revised  9/19/88 J.C Dulac  vprintf remplaced by vfprintf
 * Revised  9/16/90 D. Nichols Standardize on vfprintf
 * Revised: dave 9/17/90  Use stdarg for ANSI-C compilers
 * Revised: Bob 10/99  Introduce auxputhead
 */
#include <sepConfig.h>
#include <stdio.h>
#include    "streamlist.h"
#include <sep_pars_external.h>
#include <sep_main_external.h>


#if NeedFunctionPrototypes
#include <stdarg.h>
#else
#include <varargs.h>
#endif

/*	put a line on the header.
 *	the arguments are identical to printf conventions
 */

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
/*VARARGS1*/
int auxputhead(char *tag, char *fmt, ... )
_XFUNCPROTOEND
#else
/*VARARGS0*/
auxputhead(va_alist)	/*not Fortran callable*/
va_dcl
#endif
{
#if NeedFunctionPrototypes
#else
	char *fmt,*tag;
#endif
	va_list apdum;

streaminf *info;
info = tag_info( tag, TAG_OUT );



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
 	va_start(apdum,fmt);
_XFUNCPROTOEND
#else
	va_start(apdum);
	/* get format */
	fmt = va_arg(apdum,char *);
#endif
  if( !info->header_title ) write_title( info );
	vfprintf(info->headfile,fmt,apdum);
	fprintf(info->headfile,"\n");
	if(ferror(info->headfile))
	   {
	    perror("putch()");
	    seperr("putch() I/O error on output header\n");
	   }
	fflush(info->headfile);
	va_end(apdum);
	return 0;
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sfauxputlin(char *tag, char *line)
_XFUNCPROTOEND
#else
int sfauxputlin(tag,line)
char *line,*tag;
#endif
{

return(auxputhead(tag,"%s\n",line));
}
