/*
maketitle

USAGE
maketitle(string)

Output Parameters
string -char* Contains user, program name, and time

Description
Print program name, user, time and system in header file 

*/
/* Modified 9/21/85 stew: handle null getlogin() so rsh stuff works */
/* maketitle 12/28/85 stew: format into string */
/* maketitle 2/13/91 dave: changed for ansi and POSIX compatibility */
/* maketitle 7/18/97  bob: prototypes*/
 
#include <sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <stdio.h>
#include <time.h>

#include "sep_main_internal.h"

#if defined(SYS_UTSNAME_H)
#include <sys/utsname.h>
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
char * maketitle( char *string)
_XFUNCPROTOEND
#else
char * maketitle(string)
char *string;
#endif
{

#if defined(SYS_UTSNAME_H)
	struct utsname name;
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
_XFUNCPROTOEND
#else
	extern time_t time();
#endif
	extern char **sepxargv;
	time_t tloc; char pwnam[32]; struct passwd *uident;
	char *statstr; char hostnam[256];

	time(&tloc);
	pwnam[0] = '\0';
	statstr = getlogin();
	if(statstr != (char *) NULL) strcpy(pwnam,statstr);
	if(pwnam[0] == '\0') {
		uident = getpwuid(geteuid());
		if(uident != (struct passwd *) NULL) strcpy(pwnam,uident->pw_name);
	}

	hostnam[0]='\0';  

#if defined(SYS_UTSNAME_H)
        if( uname(&name) == 0 ){
        	strcpy(hostnam, name.nodename);
        }

#else /* not POSIX */
 	gethostname( hostnam, sizeof(hostnam) );
#endif


	sprintf(string,"\n\n%s:   %s@%s   %s",*sepxargv,pwnam,hostnam,ctime(&tloc));
	return(string);
}
