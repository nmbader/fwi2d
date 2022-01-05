/*$

=head1 NAME

fullnm  - expand filename to a fully qualified name with all path prefixes

=head1 SYNOPSIS

C<int fullnm(filename,maxlen)>


=head1 INPUT PARAMETERS

=over 4

=item filename - char* 

      filename

=item maxlen   -  int 

      maximum length

=back

=head1 DESCRIPTION

Fullnm is a utility function that filename expands a string.  This
includes generating the pwd (1) prefix for unqualified names as well
as metacharacters a la csh (1).  Up to the first 'maxlen' characters
of the output will replace the input and a count of any remaining
characters will be returned. It also attempts to strip autoumounter
paths from the directory name.

=head1 SEE ALSO

csh(1), pwd(1), sh(1), find(1), basename(1)

=head1 DIAGNOSTICS

Fullnm returns the value -1 in case of error.

=head1 BUGS
Only csh metacharacters are recognized.

=head1 KEYWORDS 

expand filename

=head1 LIBRARY

B<sep>

=cut

*/
/*
 * Revised  5/16/84 stew:  forget csh metachars ... simply use getwd()
 * revised  9-16-90  dave   made ansi-c and posix compatible
 * revised  2-19-91  carlos  remove /tmp_mnt fron the front of full pathnames.
 *			     This is necessary when using the automounter.
 * revised  7-17-92  martin  changed to posix function  getcwd()
 * revised  7-19-97  Bob   changed NeedProto
 * revised  6-1-99   Bob  Added GNU proto
 */


#include <sepConfig.h>

#include <prototypes.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

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



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static void strshft(char *, int);
_XFUNCPROTOEND 
#else
extern char *getcwd();
static void strshft();
#endif /* NeedFunctionPrototypes*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static void strshft(char *, int);
_XFUNCPROTOEND 
#else
extern char *getcwd();
static void strshft();
#endif /* NeedFunctionPrototypes*/


#include "sep_main_internal.h"
#include <sep_main_external.h>

/* need stdio.h to get definition of NULL */
#include <stdio.h>
#include <string.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int fullnm(char *string, int maxlen)
_XFUNCPROTOEND 
#else
int fullnm(string,maxlen)
char *string; int maxlen;
#endif
{

  char path[MAXPATHLEN];
  int i, mcopy, ncopy, inlen, pathlen;


  if(0 == getcwd(path,(size_t)MAXPATHLEN)) {
			 perror("fullnm:");
			 seperr("fullnm(): %d, %s\n",MAXPATHLEN,path);
  }

  /* strip /tmp_mnt off the path */
  if( !strncmp( path,"/tmp_mnt",8) ) strshft(path,8);

  /* return if already fully qualified or starts with an environment variable */
  if(string[0] == '/' || string[0] == '$' ) return(0); 

  inlen = (int)strlen(string);
  pathlen = (int) strlen(strcat(path,"/"));
  string[--maxlen] = '\0'; /* insure null padding */
  mcopy = MIN(inlen,maxlen-pathlen);
  if(mcopy > 0) {  /* null terminated string */
	for(i=mcopy-1; i>=0; --i) string[pathlen+i] = string[i];
	string[pathlen+mcopy] = '\0';
	}
  ncopy = MIN(maxlen,pathlen);
  if(ncopy > 0) memcpy(string,path,ncopy); /* prepend current directory */
  return(inlen-mcopy);
}

/* shifts string p by n characters */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static void strshft(char *p, int n)
_XFUNCPROTOEND 
#else
static void strshft(p,n)
char *p; int n;
#endif
{
char *q = p+n;
	while ((*p++ = *q++) != '\0') 
	;
}
