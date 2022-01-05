/*$

=head1 NAME

sepwarn - print a string and return a specified value

=head1 SYNOPSIS

C<int sepwarn(warn,format,value1,value2,...)>

=head1 INPUT PARAMETERS

=over 4

=item	warn -   int    

      return value to issue

=item	format - string  

      format statement

=item	value* - void    

      additional return values

=back


=head1 DESCRIPTION

This function issues a warning message  and returns (warn) to the calling 
program

From fortran you can only specify a string to print, no format control
is available.

=head1 COMMENTS
From C:

C<sepwarn(format, value1, value2, ...)>
or 
C<return(sepwarn,format, value1, value2)>
char *format

From Fortran:

C<call sepwarn('warning message')>
or
C<stat=sepwarn('warning message')>

=head1 SEE ALSO

L<seperr>

=head1 LIBRARY

B<sep>

=cut



>*/

/*
SEE ALSO

DIAGNOSTICS


KEYWORDS: error exit quit

*/
/*	
 * Create  6/99       Bob Written (hacked seperr) 
*/

/*  seperr edit history
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
int sepwarn(const int warn,const  char *format, ... )
_XFUNCPROTOEND
#else
/*VARARGS0*/
int sepwarn(va_alist)
va_dcl
#endif
{
	va_list apdum;
#if NeedFunctionPrototypes
#else
	int  warn;
	char *format;
#endif
	extern char **sepxargv;

#if NeedFunctionPrototypes
 	va_start(apdum,format);
#else
	va_start(apdum);
	format = va_arg(apdum,char *);
#endif


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

  /* print out remainder of message */
	vfprintf( stdout, format, apdum );
	fflush(stdout);

	va_end(apdum);
  return (warn);
}
