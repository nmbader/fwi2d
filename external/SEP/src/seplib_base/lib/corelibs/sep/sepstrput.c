/*
 *      Append header parameter to  output header.  
 *      used by putch, auxpar etc.
 *      Opposite to  'sepstrpar()'
 *
 * Author: dave 7/21/92  derived from putch and putch2.
 * Author: dave 7/27/94  Made sepstrputlast more efficient.
 * Revised: stew 2/25/95 Deleted unused inline[256] variable
 */
#include <sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include    <stdio.h>
#include    <sepcube.h>
#include    "streamlist.h"
#include    <assert.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int sepstrput( streaminf *info, char *name, char *type, void* val )
_XFUNCPROTOEND 
#else
int sepstrput( info, name, type, val )
streaminf *info;
char *name;
char *type;
void* val;
#endif
{
  char outline[16374];
  /* char *filename; Never Used */
  MIXED var;

  assert( info != 0 );
  assert( info->entrytype == STREAMOUT || info->entrytype == STREAMINOUT || info->entrytype == STREAMSCR  || info->entrytype == STREAMSOCKOUT );
  assert( val != 0 );
/*  assert( name != 0 );*/
  assert( type != 0 );

  if( info->headfile == (FILE*)0 ) {
	seperr("sepstrput(): Attempt to putch to invalid or closed header for tag %s", info->tagname);
  }
  /* put a title in the header if it hasn't been done yet */
  if( !info->header_title ) write_title( info );

  switch(type[0]) {
	case 'g': var.g = (double *) val; break;
	case 'm': var.m = (long long *) val ; break;
	case 's': var.s = (char *) val ; break;
	case 'f': case 'r': var.f = (float *) val ; break;
	default: var.i = (int *) val ; break;
  }

  putch_format(outline,name,type,var);
  fputs(outline,info->headfile);
  fflush (info->headfile);
  if(ferror(info->headfile)) {
	    perror("sepstrput()");
	    seperr("sepstrput() I/O error on output header for tag %s\n",
			info->tagname);
   }

   if( ((info->entrytype == STREAMINOUT)||(info->entrytype == STREAMSCR)) && info->hetch_queue != 0 ){
      /* add it to the hash table because someone may want it later */
	getpar_push_input( outline, 0 );
	getpar_scan( info->hetch_queue,info->hqlen );
   }

   return 0;

}

/* a version of sepstrput that will overwrite the last line to update
 * a parameter
 */

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sepstrputlast( streaminf *info, char *name, char *type, void* val )
_XFUNCPROTOEND
#else
int sepstrputlast( info, name, type, val )
streaminf *info;
char *name;
char *type;
void* val;
#endif
{
  /* char inline[256]; */
  char outline[256];
  char lastname[256];
  char lastline[256];
  long endpos,searchlen;
  int found,nret;
  MIXED var;

  assert( info->entrytype == STREAMOUT || info->entrytype == STREAMINOUT || info->entrytype == STREAMSCR);
  assert( info->headfile != (FILE*)0 );
  assert( !info->isapipe );
  assert( info != 0 );
  assert( val != 0 );

  /* put a title in the header if it hasn't been done yet */
  if( !info->header_title ) write_title( info );

  switch(type[0]) {
	case 'g': var.g = (double *) val; break;
	case 'm': var.m = (long long *) val ; break;
	case 's': var.s = (char *) val ; break;
	case 'f': case 'r': var.f = (float *) val ; break;
	default: var.i = (int *) val ; break;
  }

  /* the "putlast" string */
  sprintf( lastname, "putlast: %s",name);
  /* the formatted output line */
  putch_format(outline,lastname,type,var);

  /* search backwards from the end of the header file for the name
   * if we find a newline before the match we will just append the
   * information. If we find a match we will overwrite the last line
   */

  found=0;


  /* seek backwards the length of the new output line plus a little */
  searchlen = (long)MIN( strlen(outline)+10, 256 );
  if( fseek( info->headfile, (long)-searchlen, 2 ) == 0 ){

      /* now read to the end */

      nret = (int)fread( lastline, 1, (int)searchlen, info->headfile) ;
      if( nret != searchlen ){
          perror("sepstrputlast()");
          seperr("sepstrputlast() I/O error on output header for tag %s\n",
		info->tagname);
      }
      /* search for the start of the last line */
      for( endpos = 1; 
	   endpos <searchlen && lastline[searchlen-endpos-1] != '\n' ; 
	   endpos++ );

      /* search for the putlast string */
      if( endpos<searchlen && 
          strstr( &lastline[searchlen-endpos], lastname ) != 0 ){
          /* found it so overwrite */
          fseek(info->headfile, (long)-endpos, 2 );
          fputs(outline,info->headfile);
     	  fflush (info->headfile);
	  found=1;
       }
  }else{
	clearerr(info->headfile); /* forget the bad seek */
  }

  if( found == 0 ){
     /* didn't find it so just append the info */
     fseek(info->headfile, 0, 2 );
     fputs(outline,info->headfile);
     fflush (info->headfile);
  }

  if(ferror(info->headfile)) {
    perror("sepstrputlast()");
    seperr("sepstrputlast() I/O error on output header for tag %s\n",
		info->tagname);
  }


   if( ((info->entrytype == STREAMINOUT)||(info->entrytype == STREAMSCR)) && info->hetch_queue != 0 ){
      /* add it to the hash table because someone may want it later */
	getpar_push_input( outline, 0 );
	getpar_scan( info->hetch_queue,info->hqlen );
   }
	return 0;
}
