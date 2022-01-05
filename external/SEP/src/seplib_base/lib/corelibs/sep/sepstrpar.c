/*
 * int sepstrpar( streaminf *info, char* tag, char *type, void *val )
 *
 * derived from hetch and auxpar this routine gets a parameter from
 * a seplib file. The seplib file is defined by the info argument
 * This function is now called by hetch and auxpar
 */

/*
 *  source file:   ./libcube/sepstrpar.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Modified 8/4/89  Dave Nichols 
 *	   made the routine varargs, copied from Stew's stuff.
 * Revised: dave 9/17/90  Use stdarg for ANSI-C compilers
 * Revised: dave 8/94  New seplib IO functions.
 * Revised: stew 2/25/95 use SYS V instead of BSD string functions
 *                       use regecmp/regex for Solaris
 * Revised: hector 7/21/95 deleted extra "{" in ANSI_DECL of
 *                         static void make_queues( streaminf* info )
 * Revised: hector 8/3/95 changed again "{" in ANSI_DECL of
 *                         static void make_queues( streaminf* info )
 * Revised: Biondo 8/2/95 changed again position of "{" in ANSI_DECL of 
 *                         static void make_queues( streaminf* info )
 * Revised: Bob    7/18/97 prototypes
 */

#include <sepConfig.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <string.h>
#include <assert.h>


#include <sepcube.h>
#include "sepstream.h"
#include "fastpar.h"


extern char **sepxargv;

#undef REGC 
#ifdef SOLARIS 
#define REGC 1
#endif

#if defined(__APPLE__) || defined(CYGWIN)
#define re_comp regcomp
#define re_exec regexec
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static char * tetch_chop(char *command,char *history);
_XFUNCPROTOEND 
#else
static char* tetch_chop();
#endif

#define HETCH_QUEUE_SIZE 4095
#define TETCH_QUEUE_SIZE 4095

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static void make_queues( streaminf* info )
_XFUNCPROTOEND 
#else
static void make_queues( info )
streaminf* info;
#endif
{

char* headx;
char saveit;
int lhx;

/* First, determine where the chopping will be for tetch.
 * start - headx is the part that will be used for hetch and tetch
 * headx - end is the part will be used by hetch only
 */
    if(sepxargv==0) seperr("Must call initpar() at the start of your program\n");
    
    
   /*Commented for now because reg_ex slowness on LINUX and tetch not used */
   /* headx=tetch_chop(*sepxargv,info->headerbuf); */
    headx=info->headerbuf;
    /* determine size of headx */
    lhx=headx - info->headerbuf;
    /* first process headx */
    if(lhx > 0) {

	/* allocate space for the tetch queue */
	info->tetch_queue = new_queue( TETCH_QUEUE_SIZE);
	info->tqlen = TETCH_QUEUE_SIZE;
	saveit = *headx; *headx = '\0';
        getpar_push_input(info->headerbuf,0);
        getpar_scan(info->tetch_queue,info->tqlen);
	*headx = saveit;
    }

    if(*headx) {
	/* allocate space for the hetch queue */
	info->hetch_queue = new_queue( HETCH_QUEUE_SIZE);
	info->hqlen = HETCH_QUEUE_SIZE;
        if(info->headfile) if(isatty(fileno(info->headfile))) return;
	getpar_push_input(headx,0);
	getpar_scan(info->hetch_queue,info->hqlen);
    }
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int sepstrpar( streaminf *info, char* tag, char *type, void *val )
_XFUNCPROTOEND 
#else
int sepstrpar( info, tag, type, val )
streaminf *info;
char *tag, *type;
void *val;
#endif
{

 MIXED var;
 int retval;


 assert( info->entrytype == STREAMIN || info->entrytype == STREAMINOUT || info->entrytype == STREAMSCR );

 if( info->headerbuf == 0 && info->entrytype == STREAMIN ) {
   seperr("Attempt to get parameter %s from tag %s which is not a valid description file \n",tag,info->tagname);
 }

 retval = 0 ;

 if( info->hetch_queue == (hash_item**)0 ) make_queues(info);

 switch(type[0]){
   case 'g': var.g = (double*)val; break;
   case 'm': var.m = (long long*)val; break;
   case 's': var.s = (char *) val; break;
   case 'f': case 'r': var.f = (float *) val; break;
   default: var.i = (int *)val; break;
 }
 
 if( info->hetch_queue != (hash_item**)0 ) 
   retval = getpar_decode(info->hetch_queue,info->hqlen,tag,type,var);

 if(retval == 0 && info->tetch_queue!=(hash_item**)0 )
   retval = getpar_decode(info->tetch_queue,info->tqlen,tag,type,var);

 return(retval);

}

#ifdef HP700
#define  _INCLUDE_XOPEN_SOURCE
#include <regex.h>
	regex_t re;
	int rcode;
#else /* ! HP700 */
#ifdef REGC
#ifdef SOLARIS
#include <libgen.h>
#else
#include<regex.h>
#endif
extern char *__loc1;
#else /* ! SOLARIS */
#ifdef SGI
#include "re_comp.h"
#else
	extern char *re_comp();
	extern int  re_exec() ;
#endif /*SGI*/
#endif /* SOLARIS */
#endif /* HP700 */

/*
returns pointer into history at point where it first encounters the
command name in the history.  If it doesn't find it, it returns the
address of the end of the history string (i.e. a pointer to the terminating
null.  A command name is matched at the beginning of a line according to
an exact match to command:\ \ \   or the regular expression

	^\([^ \t\n]*\)/command:[ ][ ][ ]

where command is filled in appropriately.  This matches things like

	command:  or  /usr/local/command:

i.e. possible outputs of puttitle().
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static char *tetch_chop(char *command,char *history)
_XFUNCPROTOEND
#else
static char *tetch_chop(command,history)
char *command, *history;
#endif
{
	char expression[128];
	char newcommand[128];
	char *ptr, *eptr;
#ifdef REGC
	char *regptr, *regmatch;
#endif
	int match, comlen;


/*
 * first build regular expression
 */
#define WHITE " \t\n"
#define SPECIAL "\\[.^*$"

        if( command ==0 ) return history;
	
	ptr = strrchr(command,'/');
	if(ptr == ((char *) NULL)) ptr = command;
	(void) sprintf(newcommand,"%s:   ",command);
	comlen = (int)strlen(newcommand);
	(void) sprintf(expression,"^[^%s]*%c",WHITE,'/');  /*initial pathname*/
	eptr = expression+strlen(expression);
	while(*ptr) {				     /* then command name */
	    if(((char *) NULL)  != strrchr(SPECIAL,*ptr)) *(eptr++) = '\\';
	    *(eptr++) = *(ptr++);
	    }
	strcpy(eptr,":   ");			    /* then : and three blanks */

#ifndef HP700


#ifdef REGC 
	regptr = regcmp(expression, (char *) NULL);
	if(regptr == (char *) NULL) {
		fprintf(stderr,"tetch regex botch: %s\n",regptr);
		return(history+strlen(history));
		}
#else
	ptr = re_comp(expression);
	if(ptr != (char *) NULL) {
		fprintf(stderr,"tetch regex botch: %s\n",ptr);
		return(history+strlen(history));
		}
#endif
/*
 * now scan history for a match
 */
	ptr = history;
	while(*ptr) {
		if(strncmp(ptr,newcommand,comlen) == 0) {
#ifdef REGC
		   free(regptr);
#endif
		   return(ptr);
		}else{
#ifdef REGC
		regmatch = regex(regptr, ptr);
		if (regmatch == (char *) NULL) {
			while(*ptr) if ('\n' == (*(ptr++))) break;
		} else {
			free(regptr);
			return(__loc1);
		}
#else /* !SOLARIS */
		match = re_exec(ptr);
		switch(match) {
		case 0: /* move to next line of history */
		    while(*ptr) if('\n' == (*(ptr++))) break;
		    break;
		case 1: /* our match */
		    return(ptr);
		case -1: /* internal error */
		    return(history);
		}
#endif /* SOLARIS */
		}
	}
#ifdef REGC
	free(regptr);
#endif
#else
	if ((rcode=regcomp(&re, expression, REG_NEWLINE)) != 0) {
		fprintf(stderr,"tetch regex botch: %d\n",rcode);
		return(history+strlen(history));
		}
/*
 * now scan history for a match
 */
	ptr = history;
	while(*ptr) {
		if(strncmp(ptr,newcommand,comlen) == 0) {
		   return(ptr);
		}else{

		match = regexec(&re, ptr, (size_t) 0, NULL, 0);
		switch(match) {
		case REG_NOMATCH: /* move to next line of history */
		    while(*ptr) if('\n' == (*(ptr++))) break;
		    break;
		case 1: /* our match */
		    return(ptr);
		default: /* internal error */
		    return(history);
		}
		}

	}
	regfree(&re);

#endif

	return(history);
}
