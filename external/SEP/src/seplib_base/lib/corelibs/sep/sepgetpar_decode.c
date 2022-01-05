/*
 *
 *  source file:   ./filters/loclib/getpar_decode.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

/* Revised 3-8-86 stew  Added timestamp to recover older functionality
 *			when presented with multiple tags
 * Revised 9-18-86 joe  Added '1' type ... y or 1 for yes, n or 0 for no
 * Revised 9-3-92 dave  Cleaned up null entry "tag=" to return no match.
 * Revised 2-25-95 stew Use SYS V string function
 * Revised 7-18-97 bob  Prototypes==good
 * Revised 12-18-97 biondo  made type='l' equivalent to type='1' 
   for saw compatibility
   Revised  6-2-99   Bob Begun translation to GNU standards
 */
#include <sepConfig.h>
#include<string.h>
#include<strings.h>
#include <stdio.h>

#if defined (HAVE_STDLIB_H)
#include<stdlib.h>
#else
extern int atoi();
#endif /* HAVE_STDLIB  */

#include <math.h>
#include <sepcube.h>
#include "fastpar.h"

/*prototypes for functions just used in this subroutine */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static int getpar_neq(register int,register char*);
static int tag_split(register char*,register int,char**,int*);
static int getpar_getval(hash_item*,char*,MIXED);
_XFUNCPROTOEND 
#else
static int getpar_neq();
static int tag_split();
static int getpar_getval();
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
hash_item ** new_queue( int len )
_XFUNCPROTOEND 
#else
hash_item ** new_queue( len )
int len;
#endif
{
hash_item **ret;

 ret = (hash_item**) malloc(sizeof(hash_item*)*len);
        if( ret == 0 )
           seperr("unable to allocate space for hetch_queue \n");

 memset( ret, '\0',  sizeof(hash_item*)*len);
 return ret;
}


/*
 * return 1 if char c is not in string s; else return 0
 */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static int getpar_neq(register int c,register char *s)
_XFUNCPROTOEND 
#else
static int getpar_neq(c,s)
register int c;
register char *s;
#endif
{
	do {
		if(*s == c) {
			return(0); 
		}
	} 
	while(*s++);
	return(1);
}

/*
 * split off first tag "n1" from input tags "n1 nt ne"
 */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static int tag_split(register char *tag,register int tlen,char **subtag,int *sublen)
_XFUNCPROTOEND 
#else
static int tag_split(tag,tlen,subtag,sublen)
register char *tag;
char **subtag;
register int tlen;
int *sublen;
#endif
{
 register int i,j;

 for(i=0; i<tlen; i++) if(tag[i] != ' ') break;
 if(i == tlen) return(0); /* all bytes consumed */
 *subtag = tag+i;
 for(j=i+1; j<tlen; j++) if(tag[j] == ' ') break;
 *sublen = (j-i);
 return(1);
}


/*
 * take stored string value and convert it according to "type" format
 * result stored at "ptr"
 *
 * type formats:
 *               "i" or "d"  to convert to integer
 *               "r" or "f"  to convert to real
 *               "g"         to convert to double precision
 *               "s"         to keep as a string value
 *
 * stored string values may specify a vector of numerics:
 *
 *               3.0,5x3.5,-1.0,3*2.2
 * 
 * yields the result:
 *
 *               3.0,3.5,3.5,3.5,3.5,3.5,-1.0,2.2,2.2,2.2
 *
 * the function's return value will be the count (10) of items converted
 */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static int getpar_getval(hash_item *foundit,char *type,MIXED ptr)
_XFUNCPROTOEND 
#else
static int getpar_getval(foundit,type,ptr)
hash_item *foundit;
char *type;
MIXED ptr;
#endif
{
    register char *sptr, *str;
    register int ival, jval;
    register int index, endindex;
    float flt;
    long long lg;
    double dubble;
    int integer;
    index=0;
    str = foundit->val;
    ival = foundit->vlen;

    if( str == NULL ) return 0;

    if( *type == 's' ){
	memcpy(ptr.s,str,ival);
	ptr.s[ival]='\0';
	return(1);
    }

    /* all other types may have more than one value encoded */

    while(ival > 0) {
	endindex= index+1;
	    sptr = str;
	    jval = ival;
	    while(jval && getpar_neq((int) (*sptr),"*x,"))
	    	{ sptr++; jval--;}
	    if(jval > 0)
	        if(*sptr=='*' || *sptr=='x') {
	    	endindex= index+atoi(str);
	       	str= sptr+1;
	    	ival = jval-1;
	    }
	switch(*type) {
	case 'm':
		lg= atoll(str);
		while(index<endindex) ptr.i[index++]= lg;;
		break;
	case 'd':
	case 'i':
		integer= atoi(str);
		while(index<endindex) ptr.i[index++]= integer;
		break;
	case 'l':
	case '1':
		if (str[0] == 'y' || str[0] == '1' || str[0] == 'Y')
			integer = 1;
		else
			integer = 0;
		while(index<endindex) ptr.i[index++]= integer;
		break;
	case 'f':
	case 'r':
		flt= atof(str);
		while(index<endindex) ptr.f[index++]= flt;
		break;
	case 'g':
		dubble= atof(str);
		while(index<endindex) ptr.g[index++]= dubble;
		break;
	default:
		seperr("getpar() unknown conversion type %c\n",*type);
	}
	while((--ival) && ((*(str++)) != ',')); /* skip past next comma */
    }
    return(endindex);
}

/*
 * take string "n1 nt ne" and look up stored parameters with those
 * names, returning value with most recent timestamp
 */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int getpar_decode(hash_item **q,int qlen,const char *tag,char *type,MIXED val)
_XFUNCPROTOEND 
#else
int getpar_decode(q,qlen,tag,type,val)
hash_item **q;
int qlen;
char *tag, *type;
MIXED val;
#endif
{
 char *subtag;
 int sublen, count = 0 ;
 char *next, *end;
 hash_item *foundit, *saveit;
 int bigtime = -1;


 if( qlen == 0 ) return 0;

 next=(char *) tag; end = ((char *) tag)+strlen(tag);
 while(tag_split(next,(int)(end-next),&subtag,&sublen)) {
    foundit = getpar_hash_lookup(q,qlen,subtag,sublen);
    if(foundit != ((hash_item *) NULL))
	if(bigtime < foundit->timestamp) {
	    bigtime = foundit->timestamp;
	    saveit = foundit;
	    }
    next = subtag+sublen;
    }
 if(bigtime >= 0) count = getpar_getval(saveit,type,val);
 return(count);
 }
