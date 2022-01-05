#include <sepConfig.h>
#include "cfortran.h"
#include <unistd.h>
#include <sepcube.h>
#include "sep_fortran_internal.h"


/* cfortran.h handles everything except putch and auxputch which have 
 * variable argument type, we have hand crafted routines for them
 *
 * We should really handle getch etc. ourseleves as well but it never
 * makes any difference  
 */
 /*
  Modified: Dave Nichols 29/12/94, added uglified version of auxpar
            to fix string handling. 
  Modified: Stew Levin  30/12/94, added strpar Fortran-callable
    	    interface for saw to use for CHARACTER aux fetching
 Modified: Stew Levin  05/08/96, added alternate Fortran entry
             point setformat for set_format() C function.  This
             works around problem with g77 compiler generating an
             external with two trailing underscores.
  Modified: Paul hargrove 11/06/96, Added alternate FORTRAN entry points
             for all functions that contain an undescore for compatibilty
             with g77's behavior of placing a double underscore at the
             end of an external name that contains underscores.
	     FORTRAN code changed to use setformat() in place of set_format()
             can be changed back.
  Modified: Paul Hargrove 08Sep96, remove INITPAR().
  Modified: Robert G. Clapp 20July97 prototypes and erase some of the old 
                                     stufff (or at least not default)
  Modified: Bob Added auxsockout
  Modified: Bob Added sepwarn
 */




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int erexit(char*);
int ifsbrk(int);
void fexit(int);
static int sfdauxin(char*) ;
static int sfdauxout(char*) ;
static int sfdauxinout( char*); 
static int sfdauxsockout(char *); 
static int sfdauxscr(char *); 
static int sfdauxtmp(char *); 
static int sfdinput(void);
static int sfdoutput(void);
static int sfdhead(void);
static int sfdatapath(char *); 
static int sfauxputlin(char *,char *); 
#ifdef SEP_OLD
static int intalloc(int len) ;
static int sfcharpar(char*,char*,char*);
#endif /* END OF SEP_OLD STUFF */
_XFUNCPROTOEND
#else
int erexit();
int ifsbrk();
void fexit();
static int sfdinput();
static int sfdoutput();
static int sfdhead();
static int sfdauxin() ;
static int sfdauxout() ;
static int sfdauxinout( char*); 
static int sfdauxscr(char *); 
static int sfdauxtmp(char *); 
static int sfdauxsockout(char *); 
static int sfdatapath(char *); 
static int sfauxputlin(char *)
#ifdef SEP_OLD
static int intalloc(int len) ;
static int sfcharpar(char*,char*,char*);
#endif /*END OF SEP_OLD STUFF*/
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void fexit(int i) 
_XFUNCPROTOEND
#else
void fexit(i) int i;
#endif
{ exit(i); }

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int ifsbrk(int len) 
_XFUNCPROTOEND
#else
int ifsbrk(len) 
int len; 
#endif
{ return (int) fsbrk(len); }


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int erexit(char *a) 
_XFUNCPROTOEND
#else
int erexit(a) 
char* a;
#endif
{ 
seperr("%s\n",a); exit(-1); return(-1); }


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
/* wrapper functions that return (int) instead of (FILE*)  */

static int sfdinput(void)
_XFUNCPROTOEND
#else
static int sfdinput()
#endif
{
FILE *fil;
if((fil=input())==0) return(-1);
return(fileno(fil));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sfdoutput(void)
_XFUNCPROTOEND
#else
static int sfdoutput()
#endif
{
FILE *fil;
if((fil=output())==0) return(-1);
return(fileno(fil));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sfauxputlin(char *tag, char *line)
_XFUNCPROTOEND
#else
static int sfauxputlin(tag,line)
char *line,*tag;
#endif
{

return(auxputhead(tag,"%s\n",line));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sfdhead(void)
_XFUNCPROTOEND
#else
static int sfdhead()
#endif
{
FILE *fil;
if((fil=sep_head())==0) return(-1);
return(fileno(fil));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sfdauxin(char *fname) 
_XFUNCPROTOEND
#else
static int sfdauxin(fname) 
char*fname;
#endif
{
FILE *fil;
if((fil=auxin(fname))==0) return(-1);
return(fileno(fil));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sfdauxout(char *fname) 
_XFUNCPROTOEND
#else
static int sfdauxout(fname) 
char*fname;
#endif
{
FILE *fil;
if((fil=auxout(fname))==0) return(-1);
return(fileno(fil));
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sfdauxinout( char *fname) 
_XFUNCPROTOEND
#else
static int sfdauxinout(fname) 
char*fname;
#endif
{
FILE *fil;
if((fil=auxinout(fname))==0) return(-1);
return(fileno(fil));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sfdauxsockout(char *fname) 
_XFUNCPROTOEND
#else
static int sfdauxsockout(fname) 
char*fname;
#endif 
{
FILE *fil;
if((fil=auxsockout(fname))==0) return(-1);
return(fileno(fil));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sfdauxtmp(char *fname) 
_XFUNCPROTOEND
#else
static int sfdauxtmp(fname) 
char*fname;
#endif 
{
FILE *fil;
if((fil=auxtmp(fname))==0) return(-1);
return(fileno(fil));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sfdauxscr(char *fname) 
_XFUNCPROTOEND
#else
static int sfdauxscr(fname) 
char*fname;
#endif 
{
FILE *fil;
if((fil=auxscr(fname))==0) return(-1);
return(fileno(fil));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sfdatapath(char *path) 
_XFUNCPROTOEND
#else
static int sfdatapath(path) 
char*path; 
#endif
{ return (int)strlen(datapath(path)); }


#ifdef SEP_OLD
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int intalloc(int len) 
_XFUNCPROTOEND
#else
static int intalloc(len) 
int len; 
#endif
{ return (int)alloc(len); }
/* sfdclosed is defined in fdclosed.c */


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sfcharpar(char *lab,char *var,char *tag)
_XFUNCPROTOEND
#else
static int sfcharpar(lab,var,tag)
char *lab, *var, *tag;
#endif
{ return auxpar(lab,"s",var,tag); }

/* wrapper for message and exit */


FCALLSCFUN1(INT,intalloc,ALLOC,alloc,INT)
FCALLSCFUN3(INT,sfcharpar,AUXSTR,auxstr,STRING,STRING,STRING)
FCALLSCFUN2(INT,sep_send_msg,SEP_SEND_MSG,sep_send_msg,STRING,STRING)
FCALLSCFUN1(INT,sfdclosed,ISCLOSED,isclosed,INT)
#endif  /*END OF SEP_OLD STUFF DEFINED IN THIS ROUTINE*/
FCALLSCFUN1(INT,sfdauxout,AUXOUT,auxout,STRING)

FCALLSCFUN0(INT,sep_end_prog,SEP_END_PROG,sep_end_prog)
FCALLSCFUN0(INT,sep_end_prog,SEP_END_PROG_,sep_end_prog_)

FCALLSCFUN1(INT,sep_progress_prog,SEP_PROGRESS_PROG,sep_progress_prog,STRING)
FCALLSCFUN4(INT,sep_prog_stat,SEP_PROG_STAT,sep_prog_stat,STRING,INT,INT,INT)
FCALLSCFUN0(INT,sep_begin_prog,SEP_BEGIN_PROG,sep_begin_prog)
FCALLSCFUN1(INT,sep_progress_prog,SEP_PROGRESS_PROG_,sep_progress_prog_,STRING)
FCALLSCFUN4(INT,sep_prog_stat,SEP_PROG_STAT_,sep_prog_stat_,STRING,INT,INT,INT)
FCALLSCFUN0(INT,sep_begin_prog,SEP_BEGIN_PROG_,sep_begin_prog_)
FCALLSCFUN0(INT,sfdhead,HEAD,head)
FCALLSCFUN1(INT,sfdauxin,AUXIN,auxin,STRING)
FCALLSCFUN1(INT,sfdauxinout,AUXINOUT,auxinout,STRING)
FCALLSCFUN1(INT,sfdauxsockout,AUXSOCKOUT,auxsocketout,STRING)
FCALLSCFUN1(INT,sfdauxtmp,AUXTMP,auxtmp,STRING)
FCALLSCFUN1(INT,sfdauxscr,AUXSCR,auxscr,STRING)
FCALLSCFUN0(INT,sfdinput,INPUT,input)
FCALLSCFUN0(INT,sfdoutput,OUTPUT,output)
FCALLSCFUN1(INT,sfdatapath,DATAPATH,datapath,PSTRING)

FCALLSCSUB1(sep_add_doc_line,SEP_ADD_DOC_LINE,sep_add_doc_line,STRING)
FCALLSCSUB1(sep_add_doc_line,SEP_ADD_DOC_LINE_,sep_add_doc_line_,STRING)
FCALLSCSUB1(fexit,EXIT,exit,INT)
FCALLSCFUN1(INT,erexit,EREXIT,erexit,STRING)
FCALLSCFUN1(INT,tag_exists,TAG_EXISTS,tag_exists,STRING)
FCALLSCFUN2(INT,fortalloc,FORTALLOC,fortalloc,PVOID,INT)
FCALLSCFUN0(INT,noieeeinterupt,NOIEEEINTERUPT,noieeeinterupt)
FCALLSCFUN1(INT,fortfree,FORTFREE,fortfree,PVOID)
FCALLSCFUN1(INT,ifsbrk,FSBRK,fsbrk,INT)

/*PART SEP_MAIN_EXTERNAL */
FCALLSCFUN1(INT,erexit,SEPERR,seperr,STRING)
FCALLSCFUN2(INT,sepwarn,SEPWARN,sepwarn,INT,STRING)
FCALLSCFUN1(INT,auxclose,AUXCLOSE,auxclose,STRING)
FCALLSCFUN1(INT,doc,DOC,doc,STRING)
FCALLSCFUN1(INT,aux_unlink,AUX_UNLINK,aux_unlink,STRING);
FCALLSCFUN1(INT,aux_unlink,AUX_UNLINK_,aux_unlink_,STRING);
FCALLSCFUN0(INT,hclose,HCLOSE,hclose)
FCALLSCFUN5(INT,slice,SLICE,slice,PSTRING,INT,INT,INT,PVOID)
FCALLSCFUN1(INT,make_unpipe,MAKE_UNPIPE,make_unpipe,STRING)
FCALLSCFUN4(INT,snap,SNAP,snap,STRING,INT,INT,PFLOAT)
FCALLSCSUB2(mkrandom_string,MKRANDOM_STRING,mkrandom_string,STRING,PSTRING)
FCALLSCFUN3(INT,sreed,SREED,sreed,STRING,PVOID,INT)
FCALLSCFUN4(INT,sreed2,SREED2,sreed2,STRING,PVOID,INT,STRING)
FCALLSCFUN8(INT,sreed_window,SREED_WINDOW,sreed_window,STRING,PINT,PINT,PINT,PINT,PINT,INT,PVOID)
FCALLSCFUN8(INT,sreed_window,SREED_WINDOW_,sreed_window_,STRING,PINT,PINT,PINT,PINT,PINT,INT,PVOID)
FCALLSCFUN3(INT,srite,SRITE,srite,STRING,PVOID,INT)
FCALLSCSUB4(grab_history,GRAB_HISTORY,grab_history,STRING,PSTRING,INT,PINT)
FCALLSCSUB4(grab_history,GRAB_HISTORY_,grab_history_,STRING,PSTRING,INT,PINT)
FCALLSCFUN4(INT,srite2,SRITE2,srite2,STRING,PVOID,INT,STRING)
FCALLSCFUN8(INT,srite_window,SRITE_WINDOW,srite_window,STRING,PINT,PINT,PINT,PINT,PINT,INT,PVOID)
FCALLSCFUN8(INT,srite_window,SRITE_WINDOW_,srite_window_,STRING,PINT,PINT,PINT,PINT,PINT,INT,PVOID)
FCALLSCFUN4(INT,sseek_block,SSEEK_BLOCK,sseek_block,STRING,INT,INT,INT)
FCALLSCFUN4(INT,sseek_block,SSEEK_BLOCK_,sseek_block_,STRING,INT,INT,INT)
FCALLSCFUN3(INT,sseek,SSEEK,sseek,STRING,INT,INT)
FCALLSCFUN1(INT,sseekable,SSEEKABLE,sseekable,STRING)

/*PART SEP_PARS_EXTERNAL */
FCALLSCFUN3(INT,fetch,FETCH,fetch,STRING,STRING,PVOID)
FCALLSCFUN3(INT,getch,GETCH,getch,STRING,STRING,PVOID)
FCALLSCFUN1(INT,getch_add_string,GETCH_ADD_STRING,getch_add_string,STRING)
FCALLSCFUN1(INT,getch_add_string,GETCH_ADD_STRING_,getch_add_string_,STRING)
FCALLSCFUN3(INT,hetch,HETCH,hetch,STRING,STRING,PVOID)
FCALLSCFUN3(INT,putch,PUTHED,puthed,STRING,STRING,PVOID)
FCALLSCFUN3(INT,tetch,TETCH,tetch,STRING,STRING,PVOID)
FCALLSCFUN1(INT,putlin,PUTLIN,putlin,STRING)
FCALLSCFUN1(INT,sep_prog,SEP_PROG_,sep_prog_,PSTRING)
FCALLSCFUN1(INT,sep_prog,SEP_PROG,sep_prog,PSTRING)
FCALLSCFUN2(INT,separg,SEPARG,separg,INT,PSTRING)
FCALLSCFUN2(INT,sfauxputlin,AUXPUTLIN,auxputlin,STRING,STRING)
FCALLSCFUN2(INT,copy_history,COPY_HISTORY,copy_history,STRING,STRING);
FCALLSCFUN2(INT,copy_history,COPY_HISTORY_,copy_history_,STRING,STRING)
FCALLSCFUN1(INT,seploc,SEPLOC,seploc,PVOID)
FCALLSCFUN4(INT,file_position,FILE_POSITION,file_position,STRING,INT,PINT,PINT)


/*MORE OLD STUFF I NEED TO GO THROUGH */
#ifdef SEP_OLD
FCALLSCFUN1(INT,auxhclose,AUXHCLOSE,auxhclose,STRING)
FCALLSCFUN3(INT,fetpar,FETPAR,fetpar,STRING,STRING,PVOID)
FCALLSCFUN4(INT,getch2,GETCH2,getch2,STRING,STRING,PVOID,STRING)
FCALLSCFUN3(INT,getpar,GETPAR,getpar,STRING,STRING,PVOID)
FCALLSCFUN3(INT,getparin,GETPARIN,getparin,STRING,STRING,PVOID)
FCALLSCFUN1(INT,hcount,HCOUNT,hcount,INT)
FCALLSCFUN3(INT,lseek,LSEEK,lseek,INT,INT,INT)
FCALLSCFUN0(INT,modpar,MODPAR,modpar)
FCALLSCFUN1(INT,msleep,MSLEEP,msleep,FLOAT)
FCALLSCFUN0(INT,pipein,PIPEIN,pipein)
FCALLSCFUN3(INT,putch2,PUTCH2,putch2,STRING,STRING,PVOID)
FCALLSCFUN3(INT,reed,REED,reed,INT,PVOID,INT)
FCALLSCFUN3(INT,rite,RITE,rite,INT,PVOID,INT)
FCALLSCFUN4(INT,sreed2,SREED2,sreed2,STRING,PVOID,INT,STRING)
FCALLSCFUN4(INT,srite2,SRITE2,srite2,STRING,PVOID,INT,STRING)
FCALLSCFUN1(INT,putlin,PUTLIN,putlin,STRING);
#endif


/*	
 *     Wrapper for Fortran calling putch
 *
 *     Because putch doesn't have fixed argument types we can't rely on 
 *     cfortran.h to null teminate it when the argument is a string.
 *     We have to do it ourseleves.
 *
 */

#ifdef AbsoftProFortran
#define  auxpar_ AUXPAR
#define  auxputch_ AUXPUTCH
#define  putch_ PUTCH
#endif


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
#if !defined(CRAY)
/* most Fortrans pass string length as extra arguments */
int putch_( char *nptr, char *tptr, void* ptr, int nlen, int tlen, int plen )
#else
/* cray does something else, it passes a packed character descriptor */
int PUTCH( _fcd nptr, _fcd tptr, void* ptr )
#endif
_XFUNCPROTOEND
#else
/* assume sun type fortran, and char* is a suitable type for any pointer */
int putch_( nptr, tptr, ptr, nlen, tlen, plen )
   char *nptr, *tptr, *ptr;
   int nlen, tlen, plen;
#endif
{
   char *name, *type, *cvar; 
   int retcode;

#if !defined(CRAY)
   name = malloc( nlen+1); strncpy( name, nptr, nlen ); name[nlen]='\0';
   type = malloc( tlen+1); strncpy( type, tptr, tlen ); type[tlen]='\0';
#else
   int nlen, tlen, plen;
   nlen = (((((long)nptr))>>35)&0x7fffff);
   name = malloc( nlen+1 ); 
   strncpy( name,((char*)(((long)(nptr))&0xfc000000ffffffff)),nlen);
   name[nlen]='\0';
   tlen = (((((long)tptr))>>35)&0x7fffff);
   type = malloc( tlen+1 ); 
   strncpy( type,((char*)(((long)(tptr))&0xfc000000ffffffff)),tlen); 
   type[tlen]='\0';
#endif
   kill_trailing( name,' ');
   kill_trailing( type,' ');

   if( type[0] == 's' ) {  /* a string argument */ 
#if !defined(CRAY)
      cvar = malloc( plen+1); strncpy( cvar, ptr, plen ); cvar[plen]='\0';
#else
      plen = (((((long)ptr))>>35)&0x7fffff);
      cvar = malloc( plen+1 ); 
      strncpy( cvar,((char*)(((long)(ptr))&0xfc000000ffffffff)),plen); 
      cvar[plen]='\0';
#endif
      kill_trailing( cvar,' ');
      retcode =  putch( name, type, cvar );
      free( cvar );

   }else{ /* non string arguments, just pass on the pointer */
       retcode =  putch( name, type, ptr );
   }

   free( name ); free(type);
   return retcode;
}

/*	
 *     Wrapper for Fortran calling auxputch
 *
 *     Because auxputch doesn't have fixed argument types we can't rely on 
 *     cfortran.h to null teminate it when the argument is a string.
 *     We have to do it ourseleves.
 *
 */


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
#if !defined(CRAY)
/* most Fortrans pass string length as extra arguments */
int auxputch_( char *nptr, char *tptr, void* ptr, char* tagptr,
		int nlen, int tlen, int len3, int len4  )
#else
/* cray does something else, it passes a packed character descriptor */
int AUXPUTCH( _fcd nptr, _fcd tptr, void* ptr, _fcd tagptr )
#endif
_XFUNCPROTOBEGIN
#else /* not ansi */
/* assume sun type fortran, and char* is a suitable type for any pointer */
auxputch_( nptr, tptr, ptr, tagptr, nlen, tlen, len3, len4 )
   char *nptr, *tptr, *ptr, *tagptr ;
   int nlen, tlen, len3, len4;
#endif
{
   char *name, *type, *tag, *cvar; 
   int retcode;

#if !defined(CRAY)
   name = malloc( nlen+1); strncpy( name, nptr, nlen ); name[nlen]='\0';
   type = malloc( tlen+1); strncpy( type, tptr, tlen ); type[tlen]='\0';
#else
/* cray code here */
   int nlen, tlen, len3, len4;
   nlen = (((((long)nptr))>>35)&0x7fffff);
   name = malloc( nlen+1 );
   strncpy( name,((char*)(((long)(nptr))&0xfc000000ffffffff)),nlen);
   name[nlen]='\0';
   tlen = (((((long)tptr))>>35)&0x7fffff);
   type = malloc( tlen+1 );
   strncpy( type,((char*)(((long)(tptr))&0xfc000000ffffffff)),tlen);
   type[tlen]='\0';
#endif
   kill_trailing( name,' ');
   kill_trailing( type,' ');

   if( type[0] == 's' ) {  /* a string argument */ 
#if !defined(CRAY)
      cvar = malloc( len3+1); strncpy( cvar, ptr, len3 ); cvar[len3]='\0';
      tag = malloc( len4+1); strncpy( tag, tagptr, len4 ); tag[len4]='\0';
#else
      /* cray code here */
      len3 = (((((long)ptr))>>35)&0x7fffff);
      cvar = malloc( len3+1 );
      strncpy( cvar,((char*)(((long)(ptr))&0xfc000000ffffffff)),len3);
      cvar[len3]='\0';
      len4 = (((((long)tagptr))>>35)&0x7fffff);
      tag = malloc( len4+1 );
      strncpy( tag,((char*)(((long)(tagptr))&0xfc000000ffffffff)),len4);
      tag[len4]='\0';
#endif

      kill_trailing( cvar,' ');
      kill_trailing( tag,' ');
      retcode =  auxputch( name, type, cvar, tag );
      free( cvar );

   }else{ /* non string arguments, just pass on the pointer */
#if !defined(CRAY)
      tag = malloc(len3+1); strncpy( tag, tagptr, len3 ); tag[len3]='\0';
#else
      /* cray code here */
      len3 = (((((long)tagptr))>>35)&0x7fffff);
      tag = malloc( len3+1 );
      strncpy( tag,((char*)(((long)(tagptr))&0xfc000000ffffffff)),len3);
      tag[len3]='\0';
#endif
       kill_trailing( tag,' ');
       retcode =  auxputch( name, type, ptr, tag );
    }

    free( name ); free(type); free(tag);
    return retcode;
 }

 /*
  *     Wrapper for Fortran calling auxpar
  *
  *     Because auxpar doesn't have fixed argument types we can't rely on
  *     cfortran.h to null teminate it when the argument is a string.
  *     We have to do it ourseleves.
  *
  *     Here we ignore the string length (since we will fill the string in
  *     ourselves) and just grab the pointer.
  *
  */



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
#if !defined(CRAY)
/* most Fortrans pass string length as extra arguments */
int auxpar_( char *nptr, char *tptr, void* ptr, char* tagptr,
              int nlen, int tlen, int len3, int len4  )
#else
/* cray does something else, it passes a packed character descriptor */
int AUXPAR( _fcd nptr, _fcd tptr, void* ptr, _fcd tagptr )
#endif
_XFUNCPROTOBEGIN
#else /* not ansi */
/* assume sun type fortran, and char* is a suitable type for any pointer */
auxpar_( nptr, tptr, ptr, tagptr, nlen, tlen, len3, len4 )
   char *nptr, *tptr, *ptr, *tagptr ;
   int nlen, tlen, len3, len4;
#endif
{
   char *name, *type, *tag, *cvar;
   int retcode;

#if !defined(CRAY)
   name = malloc( nlen+1); strncpy( name, nptr, nlen ); name[nlen]='\0';
   type = malloc( tlen+1); strncpy( type, tptr, tlen ); type[tlen]='\0';
#else
/* cray code here */
   int nlen, tlen, len3, len4;
   nlen = (((((long)nptr))>>35)&0x7fffff);
   name = malloc( nlen+1 );
   strncpy( name,((char*)(((long)(nptr))&0xfc000000ffffffff)),nlen);
   name[nlen]='\0';
   tlen = (((((long)tptr))>>35)&0x7fffff);
   type = malloc( tlen+1 );
   strncpy( type,((char*)(((long)(tptr))&0xfc000000ffffffff)),tlen);
   type[tlen]='\0';
#endif
   kill_trailing( name,' ');
   kill_trailing( type,' ');

   if( type[0] == 's' ) {  /* a string argument */
     /* ignore the third length argument and just grab the pointer */
#if !defined(CRAY)
      cvar = ptr;
      tag = malloc( len4+1); strncpy( tag, tagptr, len4 ); tag[len4]='\0';
#else
      /* cray code here */
      cvar = ((char*)(((long)(ptr))&0xfc000000ffffffff));
      len4 = (((((long)tagptr))>>35)&0x7fffff);
      tag = malloc( len4+1 );
      strncpy( tag,((char*)(((long)(tagptr))&0xfc000000ffffffff)),len4);
      tag[len4]='\0';
#endif
      kill_trailing( tag,' ');
      retcode =  auxpar( name, type, cvar, tag );

   }else{ /* non string arguments, just pass on the pointer */
#if !defined(CRAY)
      tag = malloc(len3+1); strncpy( tag, tagptr, len3 ); tag[len3]='\0';
#else
      /* cray code here */
      len3 = (((((long)tagptr))>>35)&0x7fffff);
      tag = malloc( len3+1 );
      strncpy( tag,((char*)(((long)(tagptr))&0xfc000000ffffffff)),len3);
      tag[len3]='\0';
#endif
       kill_trailing( tag,' ');
       retcode =  auxpar( name, type, ptr, tag );
   }

   free( name ); free(type); free(tag);
   return retcode;
}
/*  $Id: fortstubs.c,v 1.3 2004/07/08 18:15:32 bob Exp $ */
