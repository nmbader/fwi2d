/*<
isapipe

USASGE
int isapipe (filedes)

INPUT PARAMETERS
filedes -int File descriptor

DESCRIPTION
Isapipe checks whether the specified file descriptor is associated
with a pipe and returns 1 if so, 0 if not.  This tells whether or not
seeks can be used for that file.

COMMENTS
SEE ALSO
lseek(2),  fseek(3), seplib, sseek

BUGS
Raw tape I/O drivers ignore seeks silently.  This situation must be
checked separately using stat(2).  Isapipe does not handle fortran
units - use the Iinquire statement to check if direct access
(i.e. seeks) are possible.


CATEGORY
Lib:Sep:I/O

COMPILE LEVEL
DISTR


>*/
/*
 * isapipe(fd)
 * 
 * checks whether a given file descriptor is associated with a pipe
 * returns 1 if so, 0 if no.
 * author - Stewart A. Levin  3/15/83
 * modified - S. Levin 3/20/83  Old algorithm failed a significant fraction
 *		of the time.  New one, much slower, uses /etc/pstat (8)
 *		to print inode and file tables which are then cross indexed
 *		with fstat() output to find the open file entry and check
 *		for the pipe flag 'P'.  One alternative could be to use
 *		ps output to look for other processes with the same parent
 *		but this could fail due to stopped or background jobs as
 *		well as a sh (1) bug that implements long pipelines as
 *		subpipes (and hence changes parents!).
 * modifed - S. Levin 3/21/83  Adapted pstat code directly in line to speed
 *		up access to the inode and file list.  Suggested by Jeff
 *		Thorson.  Also access system info once on the assumption
 *		that preallocated stdin, stdout, stderr files are the only
 *		pipes the user wouldn't know about.
 *
 * modified - S. Levin  5/6/83   mimic system utility tee by declaring
 *		errno extern and trying to seek on the file descriptor.
 *		We then check for error (-1 return value) and compare
 *		the errno value to " seek on a pipe " error number (29).
 *
 * modified - S. Levin 9/16/85   Convex returns errno 22 (EINVAL) when
 *		it should return 29 (ESPIPE) for seek on pipe.
 *
 * modified - D. Nichols 9/19/90   Convex OS8.1 now retruns the correct 
 *		 error so we need only check for (ESPIPE).
 * modified - R. Clapp  7/18/97    Removed external errno because should
 *                                 be in errno.h added include unistd.h
 *                                 so lseek explicitly defined
 *
 */
#include <sepConfig.h>
#include <stdio.h>

#include <sep_main_external.h>
#include <sys/types.h>
#include <sys/file.h>

#if defined (HAVE_ERRNO_H) || defined(__APPLE__)
#include <errno.h>
#endif

#include <unistd.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int isapipe(int fd)
_XFUNCPROTOEND
#else
int isapipe(fd)
int fd;
#endif
{

register long rc;

rc = lseek(fd,0L,1);
/*
fprintf(stderr,"isapipe: fd=%d, lseek rc=%d, errno=%d\n",fd,rc,errno);
*/
if(-1 == rc)
	 if( errno == ESPIPE ) return(1);
return(0);

}
