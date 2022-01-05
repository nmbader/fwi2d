/*$
=head1 NAME

datapath  - The datapath to put a seplib binaries

=head1 SYNOPSIS

#include <seplib.h>
char *datapath(prefix)
char *prefix

=head1 DESCRIPTION

Datapath outputs an intelligent prefix character string that can be
used to generate output file names.  This prefix is generated
according to the following priorities:

=head1 COMMENTS

1) Look for datapath= on the command line.

2) Look for the environmental variable DATAPATH

3) Look for datapath= in a file named `.datapath' in the current directory.

4) Look for datapath= in a file named `.datapath' in the user's home directory.

5) Use the default datapath `DEFAULT_DATA_DIR/login_name/' .
DEFAULT_DATA_DIR dir is defined when seplib is compiled. The default
value is often "/scr".

A ".datapath" file has the following format. An entry consists of an
(optional) host name followed by datapath=path, where path is the
datapath you wish to use. If you wish to place files in a directory you
should terminate the path component with a /.  The first entry should
contain no hostname, it is the default datapath for all hosts. If an
entry is founds with a hostname that matches the initial part of the
current hostname (without the domainname) then that datapath will be
used instead.

e.g.
datapath=/scr/dave/
oas datapath=/scr4/dave/
robson datapath=/scrx/dave

This will result in a datapath of "/scr/dave/" except on machines oas
and robson.

If a datapath starts with a $ (e.g. datapath=$SCR/dave/ ) then the
value of the environment variable will be used to start the datapath.
See "output" for a description of how environment variables are
handled in the information written to the header file.

The value returned by datapath is the address of the input string.

=head1 SEE ALSO

output, input, L<slice>, seplib

=head1 DIAGNOSTICS 

Datapath terminates program execution with an appropriate
error message if it has trouble reading a .datapath file.  

=head1 BUGS

The long list of rules is not easily remembered.  It is an attempt,
however, to mimic standard UNIX conventions such as  csh  or 
ex use to establish defaults.

=head1 KEYWORDS 

datapath output file binary

=head1 LIBRARY

B<sep>

=cut
*/
/*
  datapath(path)   -- Generate a default data routing prefix

	    if (fetch finds datapath=path)
	    else if environment variable DATAPATH exists use its contents
	    else if file .datapath exists use its contents for datapath
	    else if file ~/.datapath exists use its contents for datapath
	    else use default datapath /scr/login_name/_
 Author: S. Levin   2/21/83
	adapted code from output() subroutine
 Revised: S. Levin  3/29/83
	use getenv to find user's home directory
 Revised: S. Levin  4/4/84
	use geteuid rather than getuid
 Revised: D. Nichols  2/12/91
	Added ansi C and posix compatibility
	Added hostname dependent datapath files.
 Wes B. IBM 03-23-91 messed with temporary print statements but removed them
     copyright (c) 1991 Stanford University
 Revised: D. Nichols  1994
	 Stop looking in the header file. ( fetch() -> getch() )
	Revised R. Clapp 1997 Moved stlib.h include to sepcube.h, changed STDC to 
                       NeedFunctionPrototypes
	Revised R. Clapp 6/1/99 Switched to GNU protyping
 */

#include <sepConfig.h>

#if defined(HAVE_SYS_UTSNAME_H)
#include <sys/utsname.h>
#else
#include <sys/param.h>
#endif

#include <sepcube.h>
#include <stdio.h>
#include <pwd.h>
#include <sys/types.h>
#include <string.h>

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#else
extern char *getlogin() ;
extern struct passwd *getpwuid();
extern int geteuid();
#endif

#if defined(HAVE_LIMITS_H)
#include <limits.h>
#endif

#if !defined (HAVE_STDLIB_H)
extern char *getenv ();
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static void readpath(FILE*, char[]);
_XFUNCPROTOEND
#else
static void readpath();
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
char *datapath( char *datapth)
_XFUNCPROTOEND
#else
char *datapath(datapth)
	char *datapth;
#endif
{
#ifndef DEFAULT_DATA_DIR
	static char default_dir[] =  "/tmp/" ;
#else
	static char default_dir[] = DEFAULT_DATA_DIR ;
#endif

	FILE *pathfile; 
	char *envptr, *username, *homedir; char temp[100];
             memset(temp,'\0',sizeof(temp));

	    /* retrieve or construct a default path name */

  	       if(getch("datapath","s",datapth))
			; 
	       else if((envptr = getenv("DATAPATH")) != NULL) 
						strcpy(datapth,envptr);
	       else if((pathfile = fopen(".datapath","r")) != NULL)
						readpath(pathfile,datapth);
	       else if (
		  (NULL != (homedir = getenv("HOME")))
		     		 && 
		       (NULL != (pathfile = fopen(
		       strcat(strcpy(temp,homedir),"/.datapath"),"r")) )
						) readpath(pathfile,datapth);
	       else if ((username = getlogin()) != NULL && *username != '\0')
						sprintf(datapth,"%s/%s/_",
							default_dir,username);
	       else if ( (username = getenv("USER")) != NULL && *username != '\0')
	 	     				sprintf(datapth,"%s/%s/_",
							default_dir,username);
	       else if ((username = getpwuid(geteuid())->pw_name) != NULL
						&& *username != '\0')
		     				sprintf(datapth,"%s/%s/_",
							default_dir,username);
	return(datapth); /* return address of char array */
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static void readpath(FILE *pathfile,char datapth[])
_XFUNCPROTOEND
#else
static void readpath(pathfile,datapth)
FILE *pathfile; char datapth[];
#endif
	{


#ifndef MAXHOSTNAMELEN
#define MAXHOSTNAMELEN 255
#endif

#if defined(_POSIX_SOURCE)
#include<limits.h>

#ifdef PATH_MAX

#ifndef MAXPATH /* define MAXPATH only if not set elsewhere */
#define MAXPATH PATH_MAX
#endif

#else 

#ifndef MAXPATH /* define MAXPATH only if not set elsewhere */
#define MAXPATH _POSIX_PATH_MAX
#endif

#endif

struct utsname name;

#else /* not POSIX */

#ifndef MAXPATH /* define MAXPATH only if not set elsewhere */
#define MAXPATH 255
#endif

#endif

 	/* dimension these to maximum pathlen and max hostnamelen */
	char tmppath[MAXPATH],host[MAXHOSTNAMELEN],thishost[MAXHOSTNAMELEN];


	int hostlen;
 
	/* read the default datapath */ 
	if( 0 >=  fscanf(pathfile,"datapath=%s",datapth) )
		seperr(" Error reading .datapath file");

	/* get this hostname */

#if defined(_POSIX_SOURCE)
	if( uname(&name) != 0 ){
		perror("trying to get hostname");
                seperr(" Error calling uname ");
        }
	strcpy(thishost, name.nodename);

#else /* not POSIX */
	if( gethostname( thishost, sizeof(thishost) ) != 0 ){
		perror("trying to get hostname");
		seperr(" Error calling gethostname ");
	}
#endif

 	/* strip off domain part, use up to first period or space*/
	hostlen = (int)strcspn( thishost, ". ");
        thishost[hostlen]='\0';

	/* read any more lines in the file and se if the hostname matches */
	while( 0 <  fscanf(pathfile,"%s datapath=%s",host,tmppath) ){
		if( strcmp( thishost, host ) == 0  ){
			strcpy( datapth, tmppath );
		}
	}
	
	if(fclose(pathfile));
	return;
	}
