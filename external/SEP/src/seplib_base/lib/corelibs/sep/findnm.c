/*
findnm

DESCRIPTION
findnm \- return name of file in current directory given a file descriptor

USAGE
int findnm(fd,string,maxlen)

INPUT PARAMETERS
fd - int file id
maxlen - int maximum length for filename 
string- char[]  name of stdout file

COMMENTS

Findnm is a utility function that looks in the current directory for
the name of a file attached to the file descriptor, fd.  The resulting
filename, if found, is placed in, string, and its length is returned.
Up to maxlen characters will be copied.  A count of zero is returned
if no name could be found.

SEE ALSO
csh(1), pwd, sh(1), find(1), basename(1), cd(1), chdir(2), getwd(3)

DIAGNOSTICS
Findnm returns the value -1 in case of string overflow.

BUGS
If the file has more than one link in the current directory,
only one name will be returned.

KEYWORDS: filename find
*/
/*
 * findnm(fd,string,maxlen) - generate unqualified name associated with fd
 *			   char string[maxlen]; int fd;
 *			   returns actual number of characters
 *
 * Author: stew     - only looks in current directory
 * Revised: stew  10-18-84  check if "." is in same file system before scanning.
 * Revised: stew   2-7-85   return -1 on string overflow.
 * Revised: Dave Nichols   9/16/90  add posix compatibility
 * Revised: stew  10-14-92 force _SUN for RS6000 for non-posix compatibility
 * Revised: stew  2-25-95  Solaris mods
 * Revised: martin  9-1-95  Linux mods
 * Revised: Bob     6-1-99  Switched to GNU ifdef
 */
#include <sepConfig.h>
#include "sep_main_internal.h"
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>


#if defined(HAVE_DIRENT_H) || defined(__APPLE__)
# include <dirent.h>
#define NAMLEN(dirent) strlen((dirent)->d_name)
#else
# define dirent direct
#define NAMLEN(dirent) (dirent)->d_namlen
# ifdef HAVE_SYS_NDIR_H
#  include <sys/ndir.h>
# endif
# ifdef HAVE_SYS_DIR_H
#  include <sys/dir.h>
# endif
# ifdef HAVE_NDIR_H
#  include <ndir.h>
# endif
#endif

/*#define DIRCAST (struct direct *)*/
#define DIRCAST (struct dirent *)

#ifndef D_INO_IN_DIRENT
#define d_fileno d_ino
#endif




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int findnm(int fd, char* string, int maxlen)
_XFUNCPROTOEND
#else
int findnm(fd,string,maxlen)
int fd; char *string; int maxlen;
#endif
{
 DIR *dirp;
 struct dirent *thisdir;
 struct stat stbuf,ptbuf;
 int nchar = 0;

 if((-1 == fstat(fd,&stbuf)) || (-1 == stat(".",&ptbuf)))
 	{ perror("findnm"); return(-1); }
 string[0] = '\0';

 /* see if file descriptor refers to same file system */
 if(stbuf.st_dev != ptbuf.st_dev) return(0);

  if(((DIR *) NULL) == (dirp = opendir("."))) return(-1);
  for(thisdir = readdir(dirp); thisdir != DIRCAST NULL ;
  thisdir = readdir(dirp)) if(((int) (stbuf.st_ino)) == 
#ifdef D_INO_IN_DIRENT
	((int) thisdir->d_ino)
#else
	((int) thisdir->d_fileno)
#endif
	) break;
      
  
  if(thisdir != DIRCAST NULL) {
		nchar= (int)NAMLEN(thisdir) ;
	  strncpy(string,thisdir->d_name,maxlen);
	  }
  closedir(dirp);

  if(nchar > maxlen) return(-1);
  return(nchar);

}
