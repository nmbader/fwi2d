/*
 * 	File name expansion.
 * 	--------------------
 *
 *      char* expandnm( char* name) 
 *
 *	Attempt to expand environment variables in the filename. 
 *      This routine looks for a leading component starting with '$'
 *      If found it takes the string up to the next '/'.
 *	
 *      This is assumed to be an environment variable. It is replaced by
 *	1) The value of the environment variable
 *	2) The value returned by getch
 *	3) The value returned by  sepstrpar using the headername.
 *
 *
 *  Note this routine allocates fresh memory for the expanded name, this may 
 *  be a memory leak, however it should only be called a few times
 *  per program so the problem should be minor. If this assumption is not
 *  valid we may have to rethink the routine.
 *
 *
 * Author   7/12/91  D. Nichols : Original routine, please log any changes.
 * Revised  5/4/95   S. Levin: avoid non-posic strdup() function
 * Revised  7/18/97  R. Clapp: Prototypes, explicit casts
 * Revised  6/1/99   R. Clapp Changed to GNU Prototyping
 *
 */

#include <sepConfig.h>
#include <stdio.h>
#include <string.h>
#include "streamlist.h"
#include <sep_main_external.h>
#include <sep_pars_external.h>


/* FROM GNU Bash */
#if !defined (MAXPATHLEN) && defined (HAVE_LIMITS_H)
#include <limits.h>
#endif /* !MAXPATHLEN && HAVE_LIMITS_H */

#if !defined (MAXPATHLEN) && defined (HAVE_SYS_PARAM_H)
#include <sys/param.h>
#endif /* !MAXPATHLEN && HAVE_SYS_PARAM_H */

#if !defined (MAXPATHLEN) && defined (PATH_MAX)
#define MAXPATHLEN PATH_MAX
#endif /* !MAXPATHLEN && PATH_MAX */

/* Yecch!  Who cares about this gross concept in the first place? */
#if !defined (MAXPATHLEN)
#  define MAXPATHLEN 1024
#endif /* MAXPATHLEN */


#if !defined (HAVE_STDLIB_H)
extern char *getenv ();
#endif /* HAVE_STDLIB_H  */


/* utility function to decode leading environment variable
 * also used in output() */

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
char* envhead( char* name )
_XFUNCPROTOEND
#else
char* envhead( name )
char * name;
#endif
{
    char *envname;
    int len;

    if( name[0] != '$' ) return (char*)NULL;
    if( name[1]=='{') {
      len=(int) strcspn(name,"}");
      if( (envname = (char*)malloc( len-1)) == (char*)NULL )
	seperr("envhead() unable to allocate work buffer");
    strncpy( envname, &name[2], len-2 );
        envname[len-2] = '\0';

    }
    else{
    /* get the environment variable */
    len = (int)strcspn( name, "/" );
    if( (envname = (char*)malloc( len )) == (char*)NULL )
	seperr("envhead() unable to allocate work buffer");
    strncpy( envname, &name[1], len-1 );
    envname[len-1] = '\0';
    
    }

    return envname;
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
char* expandnm( char* name , char* hdrname )
_XFUNCPROTOEND
#else
char* expandnm( name , hdrname )
char *name;
char *hdrname;
#endif
{
streaminf *info=0;

if( hdrname != 0 ) info = tag_info( hdrname, TAG_IN );

return expand_info(  name, info );

}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
char* expand_info( char* name , streaminf* info )
_XFUNCPROTOEND
#else
char* expand_info(  name , info )
char *name;
streaminf * info;
#endif
{
    char *newname;
    char *envname;
    char enval[MAXPATHLEN];
    char *ptr;
    int len,newlen,i;

    if( ( envname = envhead(name) ) == (char*)NULL ) {
	ptr = (char*)malloc((int)strlen(name)+1);
	strcpy(ptr,name);
	return(ptr);
	}

    /* find a value for the environment variable */
    if( ( ptr = getenv(envname)) != (char*)NULL)
	strcpy( enval, ptr );
    else
	if( !getch( envname, "s", enval ) ){
	    i=0;
	    if( info != 0  ) i = sepstrpar(  info, envname, "s", enval ); 
	    if( !i ) seperr("Unable to obtain pathname variable %s\n",envname );
        }
	

   len = (int)strlen(envname);
   newlen = (int)strlen(name)+(int)strlen(enval)-len;
   if( (newname = (char*) malloc(newlen+1)) == (char*)NULL )
	seperr("expandnm() unable to allocate work buffer");

    strcpy( newname, enval ); 
    strcat(newname, &name[len+1] );


    return(newname);
}
