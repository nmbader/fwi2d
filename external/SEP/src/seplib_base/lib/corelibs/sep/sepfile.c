/*<
file

DESCRIPTION
file \- open, optionally creating if necessary, a file

USAGE
int file(name, mode)

INPUT PARAMETERS
name -char* name of file;
mode -mode  mode to open file


COMMENTS

File opens a file called name, given as the address of a
null-terminated string.  If mode=0 it will try to open it for reading.
If mode=1 it will attempt to creat() the file and open it for reading
and writing.  If mode=2 it will first try to open it for reading and
writing and will then try to create a new file if this fails.  

In every case its file descriptor is returned.

SEE ALSO 
seplib(9), reed(9s), rite(9s), close(2) 

DIAGNOSTICS 

Program execution will be terminated via seperr() if: a needed directory
is not searchable; the file does not exist and the directory in which
it is to be created is not writable; the file does exist and is
unwritable; the file is a directory; there are already too many files
open.

BUGS 

Cannot override 0664 creation mode without separate chmod() call.  No
option to return error flag instead of terminating job.

KEYWORDS: file create open.

*/

/*
 *	opens and creates file
 *	arguments:
 *		name	character string containing filename
 *		mode	0 open only
 *			1 create only
 *			2 try open, then create
 *	returns file descriptor
 *
 *  Modified Dave Nichols, 10/22/91 : Added file name expansion calls
 *  Modified Dave Nichols, 12/22/95 : Removed file name expansion calls
 *       this was a design defect and badly impacted vplot linking outside
 *       seplib.
 *   Modeified Robert Clapp, 7/15/97 : Added include file unistd, not sure
 *                                     what the if conditions hsould be
 */

#include <sepConfig.h>
#include <sepcube.h>
#include <fcntl.h>
#include <unistd.h>
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int file (char *name, int mode)
_XFUNCPROTOEND
#else
int file (name,mode)
char *name;
int mode;
#endif
{
	int filedes;
	switch (mode)
	{
	case 0:
		if ((filedes = open(name,0)) ==-1)
			{
			 perror("file()");
			 seperr("file() can't open file %s\n",name);
			}
		break;
	case 2:
		if ((filedes = open (name,2)) ==-1)
		case 1:
			if ((filedes = creat (name,0664)) ==-1)
			{
			 perror("file()");
			 seperr("file() can't create file %s\n",name);
			}
			else
			{
				if (-1 == close (filedes))
				 {
				  perror("file()");
				  seperr("file() unable to close file %s\n",name);
				 }
				if(-1 == (filedes=open(name,2)))
				 {
				  perror("file()");
				  seperr("file() unable to open file %s\n",name);
				 }
			}
	}
	return (filedes);
}
