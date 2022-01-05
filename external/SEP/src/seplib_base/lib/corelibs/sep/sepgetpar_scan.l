NONWHITE	[^ \t\n]
WHITE	[ \t\n]
ALPHA	[A-Za-z]
ALPHANUM	[A-Za-z0-9_-]
TAG	{ALPHA}{ALPHANUM}*=
SQ	\'([^'\n]*\'\')*[^'\n]*['\n]
DQ	\"([^"\n]*\"\")*[^"\n]*["\n]
%{
/* lexical scanning for fast getpar */
/* Revised 3-8-86 stew  Added time stamp to enable older method of handling
 *			multiple tags.  Moved par= intiialization to
 * 			separate routine to avoid code duplication. 
 * Revised 5-15-89 farrell  Size in suballoc forced to WORDSIZE boundary.
 * Revised 9-3-92 dave  deal with null tags properly.
 * Revised 2-25-95 stew changed bcopy to memcpy, modified renaming
 *                      of yylook->getpar_yylook, etc. for Solaris
 */
#include <sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef HAVE_CTYPE_H
#include <ctype.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef __APPLE__
#include <sys/types.h>
#include <sys/uio.h>
#include <fcntl.h>
#endif
#include "../include/fastpar.h"
#include "../include/sep_file_types.h"
#include "../include/streamlist.h"
#define input DFDFFINPUT
#include <sep_main_external.h>
#undef input
extern char *alloc(size_t);
static int massage();


#define MAX_INPUT_DEPTH 10
#if !defined(FLEX) && !defined(FLEX_SCANNER)
#undef input
#define input() ((int) *(input_stack[input_depth]++))
#undef unput
/* The redundant  =(c) insures side effects of expressions occur */
#define unput(c) (*(--(input_stack[input_depth]))=(c))
#else 
YY_BUFFER_STATE input_state[MAX_INPUT_DEPTH];
#endif

#define yywrap getpar_pop_input
#define yylex getpar_lexscan
#define yylook getpar_yylook

static int input_depth = -1;
static char *input_stack[MAX_INPUT_DEPTH];
static char *dealloc_stack[MAX_INPUT_DEPTH];

static struct {
	char *tag;  int tlen;
	char *val;  int vlen;
       } yy;


#define WORDSIZE 8
static int  SMALLBLOCK = /*4096*/ -1;
static char
*suballoc (int size)
{
	static char *myblock = (char *) NULL; static int bytesleft = 0;
	char *ptr;
	int allocsize;

	if (size % WORDSIZE == 0) {
		allocsize = size;
	} else {
		allocsize = ( (size / WORDSIZE) + 1) * WORDSIZE;
		}
	if(allocsize > SMALLBLOCK) return(alloc(allocsize));
	else
	  { 
	    if(bytesleft < allocsize)
	       {
	        myblock = alloc (SMALLBLOCK);
		bytesleft = SMALLBLOCK - allocsize;
	       }
	   else
	       {
		bytesleft -= allocsize;
	       }
	   ptr = myblock;
	   myblock += allocsize;
	   return(ptr);
	  }
}

static int  prime[10] = {31,29,23,19,17,13,11,7,5,3};
int getpar_hash(
register char *array,
register int len)
{
  register int hash;
  register int i;
  if(len >10) len=10;
  hash=0;
  for(i=0; i<len; i++)
    hash += array[i]*prime[i];
  return(hash);
}

/* workhorse to decode par files; shell already parses command line */
void getpar_scan( register hash_item **queue, register int qlen)
{
 extern int yylex(void);

 while(yylex()) {
	getpar_hash_store(queue,qlen,yy.tag,yy.val,yy.tlen,yy.vlen);
	if(yy.tlen == 3 && 0 == memcmp(yy.tag,"par",3))
		/* this is the trouble spot */
		getpar_stack_par(yy.val);
	}
}

/* read parfile into core and put buffer on scan input stack */
void getpar_stack_par(char *val)
{
 register char *buffer;
 register int fd, len;
 extern int file();
/* extern se fsize();*/
 extern char *alloc(size_t);
 ssize_t rc;

    fd = file(val,0);
    len = (int)fsize(fd);
    buffer=alloc(len+3);
    buffer[0]='\n';
    rc = read(fd,buffer+1,len);
    if(rc != len) perror("getpar_stack_par() short read");
    buffer[len+1]='\n';
    buffer[len+2]='\0';

    getpar_push_input(buffer,1);
		/* this is the trouble spot */
    close(fd);
}

 /* return 1 if match; 0 otherwise */
#define getpar_hash_compare(next1,tag1,tlen1)  \
 ((next1)->tlen == (tlen1) && 0 == memcmp((next1)->tag,tag1,tlen1))

void getpar_hash_store(
hash_item **q,
int qlen,
register char *tag, register char *val,
register int tlen,
int vlen)
{
 register hash_item *hold, *next;
 static int storetime = 0;

 hold=(hash_item *) (q+getpar_hash(tag,tlen)%qlen);
 next=hold->next;

 while(next != ((hash_item *) NULL)) {
	if(getpar_hash_compare(next,tag,tlen) ) {
		next->val = val; next->vlen = vlen;
		next->timestamp = storetime++; return;
		}
	hold = next; next = next->next;
	}

 hold->next = next = (hash_item *) suballoc((int) (sizeof(hash_item)));
 next->next = (hash_item *) NULL;
 next->tlen = tlen;
 next->tag = tag;
 next->vlen = vlen;
 next->val = val;
 next->timestamp = storetime++;
}

hash_item *getpar_hash_lookup(
register hash_item **q,
register int qlen,
register char *tag,
register int tlen)
{
 register hash_item *next;

 next = *(q + getpar_hash(tag,tlen)%qlen);

 while(next != ((hash_item *) NULL) ) {
	if(getpar_hash_compare(next,tag,tlen)) break;
	next = next->next;
	}
 return(next);
}

%}
%S FOUNDTAG
%%
<FOUNDTAG>{SQ}		{
			 yy.vlen = yyleng-2; yy.val=suballoc(yy.vlen+1);
			 yy.vlen = massage(yytext+1,yy.val,yy.vlen,yytext[0]);
			 yy.val[yy.vlen]='\0'; BEGIN 0; return(FOUNDTAG);
			 }
<FOUNDTAG>{DQ}		{
			 yy.vlen = yyleng-2; yy.val=suballoc(yy.vlen+1);
			 yy.vlen = massage(yytext+1,yy.val,yy.vlen,yytext[0]);
			 yy.val[yy.vlen]='\0'; BEGIN 0; return(FOUNDTAG);
			 }
<FOUNDTAG>[^'"]{NONWHITE}*	{
				 yy.vlen=yyleng; yy.val=suballoc(yy.vlen+1);
		 		 memcpy(yy.val,yytext,yy.vlen+1); BEGIN 0;
				 return(FOUNDTAG);
				 }
^{TAG}/{NONWHITE}	{
			 yy.tlen=yyleng-1; yy.tag=suballoc(yy.tlen+1);
			 memcpy(yy.tag,yytext,yy.tlen);
			 yy.tag[yy.tlen]='\0'; BEGIN FOUNDTAG;
			 }
([ \t]{TAG})/{NONWHITE}	{
			 yy.tlen=yyleng-2; yy.tag=suballoc(yy.tlen+1);
			 memcpy(yy.tag,yytext+1,yy.tlen);
			 yy.tag[yy.tlen]='\0'; BEGIN FOUNDTAG;
			 }
^{TAG}/{WHITE}       	{
                         yy.tlen=yyleng-1; yy.tag=suballoc(yy.tlen+1);
                         memcpy(yy.tag,yytext,yy.tlen);
                         yy.tag[yy.tlen]='\0'; yy.val=NULL; return(FOUNDTAG);
                         }
([ \t]{TAG})/{WHITE}	{ yy.tlen=yyleng-2; yy.tag=suballoc(yy.tlen+1);
			  memcpy(yy.tag,yytext+1,yy.tlen);
			  yy.tag[yy.tlen]='\0'; yy.val=NULL; return(FOUNDTAG);
			}
^\#.*	/* skip comment lines */;
.	|
\n	;
%%
	void getpar_push_input(
	register char *buffer,
	register int dealloc)
	{
	  if(input_depth++ == MAX_INPUT_DEPTH)
		seperr("too many nested par files\n");
	  input_stack[input_depth] = buffer;
#if defined(FLEX) || defined(FLEX_SCANNER)
  input_state[input_depth] =yy_scan_string(input_stack[input_depth]);
#endif
	  if(dealloc) dealloc_stack[input_depth] = buffer;
	  else dealloc_stack[input_depth] = (char *) NULL;
	}

	int
	yywrap()
	{
#if defined(FLEX) || defined(FLEX_SCANNER)
 yy_delete_buffer(input_state[input_depth] );
  if(((char *) NULL) != dealloc_stack[input_depth]) {
	free(dealloc_stack[input_depth]);
	dealloc_stack[input_depth] = (char *) NULL;
	}
 input_depth--;
 if(input_depth >=0)
	 yy_switch_to_buffer(input_state[input_depth] );
#else
	  if(((char *) NULL) != dealloc_stack[input_depth]) {
		free(dealloc_stack[input_depth]);
		dealloc_stack[input_depth] = (char *) NULL;
		}
	  input_stack[input_depth--] = (char *) NULL;
#endif
	  if(input_depth < 0){ return(1);}
	  return(0);
	}

	static int
	massage(string,out,len,quote)
	register char *string, *out;
	register int len, quote;
	{
 	register int i,j;
	
	for(i=0,j=0; i<len-1; j++) {
		out[j]=string[i++];
		if(out[j]==quote) /* compress doubled quotes */
			if(string[i]==quote) i++;
		}
	if(i<len) out[j++] = string[i];
	return(j);
	}
