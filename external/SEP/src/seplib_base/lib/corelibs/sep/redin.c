/*
redin

USAGE
int redin()

DESCRIPTION
Redin checks whether standard input is a terminal.  If not, stdin is
assumed redirected and a 1 is returned.  Otherwise a 0 is returned.

SEE ALSO
isatty(3), doc, seplib

KEYWORDS: stdin redirect
*/
/*
 *	redin() returns a 1 if the input is redirected and a 0 if it is not
 *	redout() returns a 1 if the output is redirected and a 0 if it is not.
 *
 * By redirected, it means that the standard input is coming from
 * a pipe ("|"), or from a file ("<").
 * %
 *
 * modified  1/26/83  S. Levin   rewrite
 * modified  3/2/83   S. Levin   converted to generic fileno(xxx) rather than 0,1...
 */
#include <sepConfig.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "sep_main_internal.h"


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int redin(void)
_XFUNCPROTOEND
#else
int redin()
#endif
{
	return (!isatty(fileno(stdin)));
}
