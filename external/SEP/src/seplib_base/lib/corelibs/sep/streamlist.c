
/*
 * LOGIC:
 *         Maintain a simple list of streaminf structures.
 *         Predefined tags are "in", 
 *	   "out", "head".  List is reordered on the
 *         fly to put most recently called tags at head of
 *         list.  New tags are added as necessary with a
 *         call to sepstr_open  used to initialize
 *         the entry.  No facility at present to close the
 *         files to make space.
 *
 * Revised: D.Nichols   1-8-93    SEP major revision to use new sepstream
 *                                logic. "head" is a pseudo-tag that actually
 *                                opens the "in" file.
 * Revised: Stew Levin  5/6/95    Don't use non-posix strdup()
 */
#include <sepConfig.h>
#include<stdlib.h>
#if defined(CRAY) && defined(__STDC__)
#undef __STDC__
#include <stdio.h>
#define __STDC__ 1
#else
#include <stdio.h>
#endif

#ifdef RS6000
#undef __STR__
#endif

#include <string.h>

#ifndef _NFILE
#if defined(vax) || defined(sun)
#include <sys/param.h>
#define _NFILE NOFILE
#endif
#if defined(CONVEX) 
#include <limits.h>
#define _NFILE OPEN_MAX
#endif
#endif

#include <sepcube.h>
#include "streamlist.h"


/* initial values for these pointers is zero (no streams initialized) */
static streaminf *streamlist=0;
static streaminf *streamtail=0;

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
streaminf * sepstr_head(void)
_XFUNCPROTOEND 
#else
streaminf *sepstr_head()
#endif

{ return streamlist; }

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void sepstr_del( streaminf *curr )
_XFUNCPROTOEND 
#else
void sepstr_del( curr )
     streaminf *curr;
#endif
{
    streaminf *tmp1,*tmp2;
    tmp1=curr->prev;
    tmp2=curr->next;

    if( curr == streamlist ) streamlist = tmp2;
    if( curr == streamtail ) streamtail = tmp1;

    if( tmp1 != 0 ) tmp1->next = tmp2;    
    if( tmp2 != 0 ) tmp2->prev = tmp1;

   
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void sepstr_addstart( streaminf *curr )
_XFUNCPROTOEND 
#else
void sepstr_addstart( curr )
     streaminf *curr;
#endif
{
    if( curr == streamlist )  return;

    curr->prev = 0;
    if( streamlist != 0 ) streamlist->prev=curr;
    curr->next=streamlist;
    streamlist=curr;
    if( streamtail == 0 ) streamtail=streamlist; /* the first entry */
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void sepstr_addend( streaminf *curr )
_XFUNCPROTOEND 
#else
void sepstr_addend( curr )
     streaminf *curr;
#endif
{
    if( curr == streamtail )  return;

    curr->prev=streamtail;
    if( streamtail != 0 ) streamtail->next = curr;
    curr-> next = 0;
    streamtail=curr;
    if( streamlist == 0 ) streamlist = curr;
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
streaminf * sepstr_new( char *tag, enum tagtype type )
_XFUNCPROTOEND 
#else
streaminf * sepstr_new( tag, type )
     char *tag; enum tagtype type;
#endif
{
  streaminf *newinf;

  /* temporary buffer for strings from par line, max length=4096 */
  char tmppar[4096];

  /* if this is just a query, generate an error */
  if( type == TAG_INQUIRE ){ seperr("Unknown tag \"%s\" \n",tag ); }

  newinf = (streaminf*) malloc(sizeof(streaminf));

  newinf->tagname = (char*)malloc(1+(int)strlen(tag));
  strcpy(newinf->tagname,tag);


  if( type == TAG_OUT ) newinf->entrytype = STREAMOUT;
  if( type == TAG_IN ) newinf->entrytype = STREAMIN;
  if( type == TAG_INOUT ) newinf->entrytype = STREAMINOUT;
  if( type == TAG_SOCKET ) newinf->entrytype = STREAMSOCKOUT;
  if( type == TAG_SCR ) newinf->entrytype = STREAMSCR;

  /* headername is either tagname or par associated with the tag */
  /* the special cases "in" and "out" have parnames "stdin" and "stdout" */
  /* headername is either a filename, a pipe command or a machine:port pair */
  
  if(strcmp(tag,"in") == 0){
      strcpy(tmppar,"stdin");
      getch("stdin","s",tmppar);
  }else if(strcmp(tag,"out") == 0){
      strcpy(tmppar,"stdout");
      getch("stdout","s",tmppar); 
  }else{
      strcpy(tmppar,tag); 
      if( type == TAG_IN || type == TAG_INOUT ){
         fetch(tag,"s",tmppar); 
      }else{
         getch(tag,"s",tmppar); 
      }
  }

  newinf->headername=(char*)malloc(1+(int)strlen(tmppar));
  strcpy(newinf->headername,tmppar);
/*  newinf->headfile=(FILE*)0;*/
/*  newinf->headfile=SEPPOINTNULL;*/
  newinf->headfile=0;

  /* initialize everything to sensible null values */
  newinf->valid = 1;
  newinf->headerbuf=0; newinf->hdrlen=0;
  newinf->hetch_queue=0; newinf->hqlen=0;
  newinf->tetch_queue=0; newinf->tqlen=0;
  newinf->backseekable = 0;
  newinf->header_title=0;
  newinf->header_outname=0;
/*  newinf->headerformatfile=alloc(1);*/
  newinf->headerformatfile=0;
  newinf->gridformatfile=0;
  newinf->key_bytes=0;
  newinf->ready_out=0;
/*  newinf->dataname=alloc(1);*/
  newinf->dataname=0;
  newinf->format_num=-1;
  newinf->iotype=FILE_IO;
  newinf->ioinf=0;
  newinf->streamfile=0;
/*  newinf->streamfile=SEPPOINTNULL;*/
  newinf->streamfd=-1;
  newinf->open_func=0;
  newinf->close_func=0;
  newinf->read_func=0;
  newinf->write_func=0;
  newinf->seek_func=0;
  newinf->size_func=0;
  newinf->isapipe=0;
  newinf->sockfd=0;


	newinf->n_key=-1;
  return newinf;

}
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void debug_tag( char *tag )
_XFUNCPROTOEND 
#else
void debug_tag( tag )
char *tag;
#endif
{
streaminf *info;
info = tag_info( tag, TAG_INQUIRE );
if(info==0) seperr("debug_tag request on an unintialized tag=%s \n", tag);
  print_streaminf(info);

return;
}

#if NeedFunctionPrototypes
extern char* get_format_name( int num);
_XFUNCPROTOBEGIN 
void print_streaminf( streaminf *info )
_XFUNCPROTOEND 
#else
extern char* get_format_name();
void print_streaminf( info )
streaminf* info;
#endif
{

fprintf(stderr," tag: %s \n",info->tagname);

switch( info->entrytype ){
    case( STREAMIN ):
       fprintf(stderr," entrytype: STREAMIN \n");
    break;
    case( STREAMOUT ):
       fprintf(stderr," entrytype: STREAMOUT \n");
    break;
    case( STREAMINOUT ):
       fprintf(stderr," entrytype: STREAMINOUT \n");
    break;
    case( STREAMSOCKOUT ):
       fprintf(stderr," entrytype: STREAMSOCKOUT \n");
    break;
    case( STREAMSCR ):
       fprintf(stderr," entrytype: STREAMSCR \n");
    break;
}

fprintf( stderr," headername: %s \n",info->headername);
fprintf( stderr," headerlen: %d \n",info->hdrlen);
fprintf( stderr," headerbuf: %s \n",info->headerbuf);

fprintf( stderr," headerformatfile: %s \n",info->headerformatfile);
fprintf( stderr," gridformatfile: %s \n",info->gridformatfile);

fprintf( stderr," backseekable: %d \n",info->backseekable);

fprintf( stderr," header_title: %d \n",info->header_title);
fprintf( stderr," ready_out: %d \n",info->ready_out);
fprintf( stderr," isapipe: %d \n",info->isapipe);

fprintf( stderr," dataname: %s \n",info->dataname);
switch( info->iotype ){
    case( FD_IO ):
       fprintf(stderr," iotype: FD_IO \n");
    break;
    case( MULTI_FD_IO ):
       fprintf(stderr," iotype: MULTI_FD_IO \n");
    break;
    case( FILE_IO ):
       fprintf(stderr," iotype: FILE_IO \n");
    break;
    case( CM_SDA_IO ):
       fprintf(stderr," iotype: CM_SDA_IO \n");
    break;
    case( CM_REG_IO ):
       fprintf(stderr," iotype: CM_REG_IO \n");
}


fprintf( stderr," format: %s \n", get_format_name(info->format_num) );
}
