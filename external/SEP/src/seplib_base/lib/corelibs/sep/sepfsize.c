/*
fsize
Usage
sep_file_size_t fsize(filedes)

INPUT PARAMETERS
filedes - int filescriptor

DESCRIPTION
 This routine finds the length of a file and returns its size in number
of bytes. filedes is the file descriptor that is returned by a open()
command or a creat() command. 

*/
#include <sepConfig.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "streamlist.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
sep_file_size_t fsize (int filedes)
_XFUNCPROTOEND 
#else
sep_file_size_t fsize (filedes)
int filedes;
#endif
/*
 *	returns size of file in bytes
 *
 * modified 1/26/83  S. Levin corrected length and simplfied
 * modified 3/27/83  return -1 if fstat fails.
 * modified 12/22/95 D. Nichols Made ssize a separate file.
 */
{
	struct stat buf;
	sep_file_size_t size;

	if(0 != fstat (filedes,&buf)) return(-1);
	size=(sep_file_size_t) buf.st_size;
	return size;
}
