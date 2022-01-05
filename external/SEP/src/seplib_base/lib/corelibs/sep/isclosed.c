/* 
 isclosed

USAGE
isclosed(stream)

INPUT PARAMETERS
stream -FILE file to check

DESCRIPTION
Check if stream is closed and available for reuse
 
*/
/*
 * Author: Stewart A. Levin  5/15/84	Grabbed code from auxin
 * Revised: Stewart A. Levin 2/26/95    SYS V string functions
 * Revised: Robert G. Clapp 7/18/97    Added prototyping
 & Revised: Robert G. Clapp 6/2/99  GNUed it
 */
#include <sepConfig.h>
#include <stdio.h>

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#include <sys/stat.h>




#if defined(HAVE_ERRNO_H) || defined(__APPLE__)
#include <errno.h>
#endif



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int isclosed(FILE* stream)
_XFUNCPROTOEND
#else
int isclosed(stream)
FILE *stream;
#endif
{
    struct stat stbuf; int ok = 1;
    if( -1 == fstat(fileno(stream),&stbuf)) {
			if (errno == EBADF) ok = 0;
			else{
				/*we will assume that this also means a bad file for now */
        ok=0;
			}
	}
  return(!ok);
}
