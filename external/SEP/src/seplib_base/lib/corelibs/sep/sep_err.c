/*$

=head1 NAME

seperr - print line and exit with a failure code

=head1 SYNOPSIS

C<void seperr(string)>

=head1 INPUT PARAMETERS

=over 4

=item	char* -  string 

      C print statement, see comments

=back


=head1 DESCRIPTION

This function permits the programmer to terminate program execution
and print an error message using printf (3) format control.  The
message is prefixed with the program name to prevent confusion.

From fortran you can only specify a string to print, no format control
is available.

=head1 COMMENTS
From C:

C<seperr(format, value1, value2, ...)>
char *format

From Fortran:

C<call seperr('error message')>

=head1 SEE ALSO

printf(3), seplib, L<sepwarn>

=head1 DIAGNOSTICS

The system exit code used is -1.

=head1 BUGS

Fortran programmers can only pass string values.

=head1 KEYWORDS 

error exit quit

=head1 LIBRARY

B<sep>

=cut



*/
/*	
 *	
 * Revised  4/23/84    stew     echo to stdout (header) as well
 * Revised  9/7/87     stew	<varargs> for portability
 * Revised  9/16/88    jclau    fvprintf instead to vprintf.
 * Modified 9/16/90  Dave Nichols  
 *	   all compilers now have vfprintf or it is in compatability library.
 * Revised: dave 9/17/90  Use stdarg for ANSI-C compilers		
 * Revised: dave 10/9/91  renamed seperr
 */
#include <sepConfig.h>
#include <stdio.h>
#include <sep_main_external.h>
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
#include <stdarg.h>
_XFUNCPROTOEND
#else
#include <varargs.h>
#endif


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
/*VARARGS1*/
int seperr(const char *format, ... )
_XFUNCPROTOEND
#else
/*VARARGS0*/
int seperr(va_alist)
va_dcl
#endif
{
	va_list apdum;
#if NeedFunctionPrototypes
#else
	char *format;
#endif
	extern char **sepxargv;

#if NeedFunctionPrototypes
 	va_start(apdum,format);
#else
	va_start(apdum);
	format = va_arg(apdum,char *);
#endif
      sep_err_prog();

/*----- stderr -----*/

        /* print out name of program causing error */
	if( sepxargv ) fprintf( stderr, "%s: ", sepxargv[0] );

        /* print out remainder of message */
	vfprintf( stderr, format, apdum );
	fflush(stderr);

	va_end(apdum);

/*----- stdout -----*/

#if NeedFunctionPrototypes
 	va_start(apdum,format);
#else
	va_start(apdum);
	format = va_arg(apdum,char *);
#endif
        /* print out name of program causing error */
	if( sepxargv ) fprintf( stdout, "%s: ", sepxargv[0] );

        /* print out remainder of message */
	vfprintf( stdout, format, apdum );
	fflush(stdout);

	va_end(apdum);
	exit(23);
}
