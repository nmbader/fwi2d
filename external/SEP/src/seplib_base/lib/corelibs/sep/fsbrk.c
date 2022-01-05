/*
 * sbrk() with error checking for use with saw, etc.
 *
 * use malloc on POSIX systems because sbrk is not guaranteed.
 */
/*
Edited Bob 6/1/99 Removed chech sbrk
*/

#include <sepConfig.h>
#include <sys/types.h>
#include <stdlib.h>
#include <sepcube.h>
#include "sep_fortran_internal.h"


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
char *fsbrk( unsigned int nbytes)
_XFUNCPROTOEND 
#else
char *fsbrk(nbytes)
unsigned int nbytes;
#endif
{

 char *rc;

 if( nbytes == 0 ) return (char*)0;

#if 1
/* CRAY gets void pointer fomr malloc */
 rc = (char *) malloc( nbytes );
 if(rc == ((char *) 0 )) {
#else
 extern caddr_t sbrk();

 rc=(char *) sbrk(nbytes);
 if(rc == ((char *) -1)) {
#endif
	perror("fsbrk");
	seperr("unable to allocate %d bytes\n",nbytes);
	}
 return(rc);
}
