/* 

  sepstr_in_head( info )     
  
  info is a sepstream info structure with the tag field , entrytype,
  and the headername filled in.
  
  This routine prepare the stream for input. 
  
  There are four possibilites for the headername specification,
  1) If it contains a "|" (pipe symbol) it is a command to use with popen.
  2) If it contains a colon it is a hostname:portnumber pair for opensock.
     (Or a unix::pathname for Unix domain opensock.)
  3) It is "stdin".
  4) If none of the above are true it is a filename to be opened. 

  On return the following entries in the structure should be filled in:
	headfile, headerbuf, hdrlen, hetch_queue, tetch_queue
	valid, dataname, format_num and all the IO function pointers

  NB this does not actually open the data file, that is deferred until
     it is needed.

*/

/*
 modified Dave Nichols 9/1/94, deal with noheader properly 
 modified Stew Levin   5/6/95  Don't use non-posix strdup()
 modified Stew Levin   6/2/99  Replace ifdef to to GNU standard
 */

#include <sepConfig.h>
#ifdef RS6000
#undef __STR__
#endif

#if defined(CRAY) && defined(__STDC__)
#undef __STDC__
#include <stdio.h>
#define __STDC__ 1
#else
#include <stdio.h>
#endif
#include <string.h>

#if defined(HAVE_ERRNO_H) || defined(__APPLE__)
#include <errno.h>
#else

#ifndef STDC_HEADERS
extern int errno;
#endif
#endif


#include <sys/types.h>
#include <sys/stat.h>

#include <assert.h>
#include <unistd.h>

#if defined(__APPLE__) || defined(LINUX)
#define USE_SOCKETS
#endif
#include "streamlist.h"
#include "sep_main_internal.h"
#include <sepcube.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static void open_infile(streaminf*);
static void open_inpipe(streaminf*);
static void open_insok(streaminf*);
static void openstdin(streaminf*);
_XFUNCPROTOEND
#else
static void open_infile();
static void open_inpipe();
static void open_insok();
static void openstdin();
#endif


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void sepstr_in_head( streaminf* info )
_XFUNCPROTOEND
#else
void sepstr_in_head( info )
     streaminf *info;
#endif
{
    /* char titleline[256]; Never Used */

    assert( info->entrytype == STREAMIN ); 
   
    if( !noheader() || strcmp(info->tagname, "in")!=0 ){
       /* input header */
       if( strchr( info->headername, '|' ) != 0 )              /* pipe  */
         open_inpipe( info ); 
       else if ( strchr( info->headername, ':' ) != 0 )  {     /* socket */
         open_insok( info );
				}
       else if( strcmp( info->headername, "stdin" ) == 0 ){     /* stdin */
         openstdin( info );
       }else                                                    /* file */
         open_infile( info );
    
			
       if( info->headfile == 0 ){
   	 info->valid = 0;
	 return;
       }

    /* make sure the input stream is fully buffered */
    setvbuf(info->headfile,0,_IOFBF,BUFSIZ); 

    }

    /* read the header into the buffer for this input stream */
    readhdr( info );
    
    /* open the input stream */
    open_instream( info );
}


/* read the header from an input header file 
 * and store in the buffer (info->headerbuf), we don't initialize
 * the hashed table of parameters until the first parameter is
 * requested. (see sepstrpar.c)
 */

/* max length of header allowed down a (non file) stream */
#define STDINHEAD 3000
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void readhdr(streaminf* info)
_XFUNCPROTOEND
#else
void readhdr( info )
     streaminf *info;
#endif
{
    /* int count; Never Used */
    int flen, next;
    /* char *hptr; Never Used */
    int headfd, len;
    char *hdrbuf;
    int it1, it2, it3;

		len=STDINHEAD;
		hdrbuf=(char*)malloc(sizeof(char)*len);
    


    assert( info->entrytype == STREAMIN || info->entrytype == STREAMINOUT );

    it1 = noheader();
    it2 = strcmp(info->tagname, "in");
    if( it1  && (it2 == 0) ){ 
	info->headerbuf = (char*)malloc(18);
	strcpy( info->headerbuf, "No input header \n");
	info->hdrlen=16;
				free(hdrbuf); return;
    }

    headfd=fileno( (info->headfile) );

    it3 = isatty(headfd);
    if( it3 ){
	info->headerbuf = (char*)malloc(18);
	strcpy( info->headerbuf, "Input from a tty\n");
	info->hdrlen=16;
				free(hdrbuf); return;
    }

    
    

    /* the header may be coming down a stream so we can't know it's size */
    /* use the buffer "hdrbuf" this has a maximum size of len */
	
    /* read until we get an end of file or the EOT mark */
    /* we are using bufferd I/O so this shouldn't be too expensive */
    flen = 0;
    while( (next = getc(info->headfile)) != EOT && next != EOF ){
        /* read until we get an EOT or EOF */

				 if(flen==len-2){
					len=len*2;
					hdrbuf=(char*)realloc(hdrbuf,len*sizeof(char));
				}
        hdrbuf[flen]=next;
        flen++; 
    }
    hdrbuf[flen]='\0';
        	
    /* now allocate storage for this much header in the streaminf
     * structure and copy the header into it.
     */
    info->headerbuf = (char*)malloc( flen + 1 );
    memcpy( info->headerbuf, hdrbuf, flen+1 );
    info->hdrlen=flen;
	
    /* see if the last character was the EOT mark, if so we must
     * perform interprocess syncronization 
     */

    if(flen > 0 && next == EOT){
        syncin( info );
    }




	free(hdrbuf); return;
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static void open_inpipe(streaminf* info)
_XFUNCPROTOEND
#else
static void open_inpipe( info )
     streaminf *info;
#endif
{
    int len;
    len = (int)strlen( info->headername );
    if( info->headername[len-1] != '|'){
	seperr("input pipe command '%s' doesn't end with '|' for tag %s \n",
	       info->headername,info->tagname);
    }else{
	info->headername[len-1] = '\0';
	info->headfile = popen( info->headername, "r" );
	if( info->headfile == 0 ) {
	    perror( "sepstrhead:openpipe() ");
	    seperr( "error in opening input pipe \"%s\" for tag \"%s\" \n",
                   info->headername,info->tagname);
	}
    }
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static void open_insok(streaminf* info)
_XFUNCPROTOEND
#else
static void open_insok( info )
     streaminf *info;
#endif
{
char *hostname, *portstr;
char *tmp;
int infd,sock,unixdom;
    
   tmp = (char*)malloc( 1 + (int)strlen( info->headername)  );
   strcpy(tmp,info->headername  );

   /* info->headername is a colon separated host:portnumber pair */ 
   /* If there is no hostname we just open up a port and wait for
       someone else to connect */
   /* if there are two colons we should open in the unix domain */

    unixdom=0;

    if( tmp[0] == ':' ){
        tmp++; hostname = 0;
	if( *tmp == ':' ){ unixdom=1; tmp++;} /* two colons => unix domain */
        portstr = strtok( tmp, ":" );
    }else{
        hostname = strtok( tmp, ":" );
        portstr = strtok( 0, ":" );
    }

    if( hostname == 0 ){
        sock = opensock1( portstr, unixdom );
        assert( sock >0 );
        infd = socklisten( sock, 30 );
    }else{
	int i;
	/* try to connect 30 times, at 1 second intervals */
	for(i =0; i<30;i++) {
          infd = opensock2( hostname, portstr );
	  if( infd != -1 ) break;
	  sleep(1);
	}
	
    }

    if( infd == -1 ){
        seperr("open_insok(): connect to socket \"%s\" failed for tag \"%s\"\n",
                      info->headername,info->tagname);
    }

    info->headfile = fdopen( infd, "r" );
    if( info->headfile == 0 ) {
        perror( "sepstrhead:open_insok() ");
        seperr( "error in opening input socket \"%s\" for tag \"%s\" \n",
                      info->headername,info->tagname);
    }

}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static void open_infile(streaminf* info)
_XFUNCPROTOEND
#else
static void open_infile( info )
     streaminf *info;
#endif
{
    struct stat statbuf;
    
    /* first stat the header file name */
    if(-1 == stat(info->headername,&statbuf)) { 
	/* header does not exist */
	if(errno != ENOENT) {  /* invalid name */
	    perror("sepstrhead, openfile()");
	    seperr("Bad header name '%s' for tag %s \n",
		   info->headername,info->tagname);
	}
	/* the name is OK but it doesn't exist */
	info->headfile =0;
	
    }else{ /* header exists */
	/* open for read */
	info->headfile = fopen(info->headername,"r");
    }
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static void openstdin(streaminf* info)
_XFUNCPROTOEND
#else
static void openstdin( info )
     streaminf *info;
#endif
{
	
    if( isatty(fileno(stdin)) ) {
	info->headfile = stdin;
    }else{
	info->headfile = stdin;
    }
}
