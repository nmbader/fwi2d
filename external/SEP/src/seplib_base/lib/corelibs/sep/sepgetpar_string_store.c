/*
 *
 *  source file:   ./filters/loclib/getpar_string_store.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

/* Rewrite 3-8-86 stew  Don't invoke parser unless we see par=.  This
 *			leaves the shell's parsed output of the command
 *			line intact, fixing joe's complaint.  Also a
 *			bit faster.
 *	   2-15-88 joe  Should be hash_item "**q", not "*q". Oops.
 *         11-2-90 stew Commented out "if(vlen<=0)return" to allow
 *			removal or erasure of preexisting definitions.
 *	   4-24-92 dave if it has quotes in it let getpar_scan figure out
 *			what to do.
 *     7-18-97 bob prototypes, declare voids, etc
 */
#include <sepConfig.h>
#include "fastpar.h"
#include <stdio.h>

#include <sepcube.h>
#include <string.h>
#include <strings.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void getpar_string_store(hash_item **q,int qsize,register char *string)
_XFUNCPROTOEND 
#else
void getpar_string_store(q,qsize,string)
hash_item **q;
int qsize;
register char *string;
#endif
{
 register int tlen, vlen;
 register char *val, *val2;

 val = strchr(string,'\"');
 val2 = strchr(string,'\'');

 if( val != 0 || val2 != 0 ) {
    /* if it has quotes in it, let getpar_scan figure it out */
    getpar_push_input( string, 0 );
    getpar_scan(q,qsize);
    return;
 }

 /* first check if this is a par=val string */
 val = strchr(string,'=');
 if(val == ((char *) NULL)) return;
 tlen = (int)(val - string);
 if(tlen <= 0) return;
 val++;
 vlen=(int)strlen(val);
	
 if(vlen <= 0) {
	/* empty entry gets a value of NULL */
	getpar_hash_store(q,qsize,string,NULL,tlen,0);
 }else{
     /* process quotes */
     getpar_hash_store(q,qsize,string,val,tlen,vlen);
     if( tlen==3 && 0==bcmp(string,"par",3) ) {
 	getpar_stack_par(val);
	getpar_scan(q,qsize);
	}
 }
}
