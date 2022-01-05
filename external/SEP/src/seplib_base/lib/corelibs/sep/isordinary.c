/*<
isordinary

USAGE
isordinary(path)

INPUT PARAMETERS
path - char path to file
 
DESCRIPTION
check if path is ordinary file

*/

/*
 *
 * Author: S. Levin  5/16/84
 */
#include <sys/types.h>
#include <sys/stat.h>

#include "sep_main_internal.h"
/*<
fdordinary

USAGE
fdordinary(path)

INPUT PARAMETERS
path - char path to file
 
DESCRIPTION
check if file descriptor is attached to an ordinary file

*/

#include <sepConfig.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int fdordinary(int fd )
_XFUNCPROTOEND
#else
int fdordinary(fd)
int fd;
#endif
{
 struct stat stbuf;

 if(0 != fstat(fd,&stbuf)) return(0);  /* extraordinary if unable to stat!!! */

#ifdef _POSIX_SOURCE
 if( S_ISREG(stbuf.st_mode) ) return(1);
#else
 if((stbuf.st_mode & S_IFMT) == S_IFREG) return(1);
#endif

 return(0);
}

/*<
isordinary

USAGE
isordinary(path)

INPUT PARAMETERS
path - char path to file
 
DESCRIPTION
check if path is ordinary file

*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int isordinary(char *path)
_XFUNCPROTOEND
#else
int isordinary(path)
char *path;
#endif
{

 struct stat stbuf;

 if(0 != stat(path,&stbuf)) return(0);  /* extraordinary if unable to stat!!! */

#ifdef _POSIX_SOURCE
 if( S_ISREG(stbuf.st_mode) ) return(1);
#else
 if((stbuf.st_mode & S_IFMT) == S_IFREG) return(1);
#endif

 return(0);
}
