/*
isatape

USAGE
int isatape (filedes)

filedes - int file descriptor

DESCRIPTION
Isatape checks whether the specified file descriptor is associated
with a tape and returns 1 if it is a raw tape, 2 if it's a cooked
tape, 0 otherwise.  A 1 implies that seeks cannot be used for that
tape.  (ioctl statements would have to be used instead).  Isatape also
can be used to recover from I/O errors.

SEE ALSO
lseek(2), ioctl(2), mt(4), seplib

BUGS
Isatape does not handle fortran units - use the inquire statement
to check if direct access (i.e. seeks) are possible.
I have never seen cooked tape used successfully.

KEYWORDS: tape seek
*/
/*
 * isatape(fd)  -  check if file descriptor corresponds to a tape drive
 *		   returns 1 if tape, 0 otherwise.
 */
/*
 * Author - Stewart A. Levin  8/29/83
 * Revised - Stewart A. Levin  11/22/85  for Convex
 * Revised - Stewart A. Levin  2/7/86  system independence
 * Revised - Dave Nichols      9/19/90 Posix compliance
 * Revised - Bob Clapp         10/6/98 for Linux
 * Revised - Bob Clapp          6/1/99 GNU ifdef
 *
 */
#include <sepConfig.h>

#if defined(HAVE_SYS_TYPES_H)
#include <sys/types.h>
#endif
#ifdef __APPLE__
#include <unistd.h>
#endif

#include "sep_main_internal.h"
#include "sep_main_external.h"


#ifdef HAVE_ERRNO_H
#include <errno.h>
#else
extern int errno;
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int isatape(int fd)
_XFUNCPROTOEND
#else
int isatape(fd)
int fd;
#endif
{

/* 
 * strict posix does not have the MTIOCGET ioctl 
 * but we can do it with convex, sun4, dec3100, rs6000 and linux
 *
 * Add other machine deopendent code here,
 *
 */

#if defined(HAVE_SYS_IOCTL_H) &&  !defined(SGI) && !defined(__APPLE__)
#include <sys/ioctl.h>


#if defined(HAVE_SYS_TAPE_H)
#include <sys/tape.h>
#else

#include <sys/mtio.h>
#if defined(HAVE_ERRNO_H)
#include <errno.h>
#endif

 struct mtget buf;
 if(-1 == ioctl(fd,(int) MTIOCGET,(char *) &buf))
	if(errno == EBADF)
		seperr("isatape: %d is not a valid file descriptor\n",fd);
	else return(0);
 if(buf.mt_type == ((short) 0)) return(0);
 return(1);

#endif

#else

 /* 
  * attempt to use only posix calls
  *
  * assume anything that is a block or character special file and 
  * isn't a tty, a directory, a pipe or a regular file is a tape 
  */

#include <sys/stat.h>

struct stat buf;

if( fstat(fd,&buf) != 0 ) return 0;
 
 if(  isatty(fd) || S_ISREG(buf.st_mode) ||  S_ISFIFO(buf.st_mode) || 
	S_ISDIR(buf.st_mode) ) return(0);

 if( S_ISCHR(buf.st_mode) || S_ISBLK(buf.st_mode) ) return(1);
 return(0);
#endif

}
