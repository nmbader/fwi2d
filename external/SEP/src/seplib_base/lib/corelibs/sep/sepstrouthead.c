#include <sepConfig.h>
 
#ifdef RS6000
#undef __STR__
#endif

#include <stdio.h>
#include <string.h>

#if defined(HAVE_ERRNO_H) || defined(__APPLE__)
#include <errno.h>
#else

#ifndef STDC_HEADERS
extern int errno;
#endif
#endif

#if defined(__APPLE__) || defined(LINUX)
#define USE_SOCKETS
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>   
#include <unistd.h>




#include "strformats.h"
#include "streamlist.h"
#include "sep_main_internal.h"
#include <sepcube.h>

#if !defined (HAVE_STDLIB_H)
extern char *getenv ();
#endif /* HAVE_STDLIB_H  */

/* 
  sepstr_out_head( info )     
  
  info is a sepstream info structure with the tag field , entrytype,
  and the headername filled in.
  
  This routine will open the header for  output and 
  fill the headfile and the headbuf, members of the
  structure.
  
  It will also fill in a default value to info->format_num, this can be
  reset by "set_output_format( tag, format )".
  
  There are four possibilites for the headername specification,
  1) If it contains a "|" (pipe symbol) it is a command to use with popen.
  2) If it contains a colon it is a hostname:portnumber pair for opensock.
  3) It is "stdout".
  4) If none of the above are true it is a filename to be opened. 
  
  An output stream should not be used  until sepstr_hclose() is called. 
  This happens in the following cases.
  
  1) someone calls hclose() or auxclose()
  2) someone calls srite()
  
  No more info can be put to the output header after sepstr_hclose is called
  if the header and data go the same routes.
  
  */
  /* 
  Author: Dave??
  Modified:
     Bob - 7-97 :  Added new style prototypes
     Bob - 12-97:  Added error checking for STREAMSCR
     Bob - 6-99:   Change to GNU ifdef
  */

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static void open_outpipe( streaminf*);
static void  open_outsok( streaminf*);
static void  openstdout( streaminf*);
static void  open_outfile(streaminf*);
_XFUNCPROTOEND 
#else
static void  open_outpipe();
static void  open_outsok();
static void  openstdout();
static void  open_outfile();
#endif

static char headername[4096];
static char parambuf[4096];


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void sepstr_out_head( streaminf *info )
_XFUNCPROTOEND 
#else
void sepstr_out_head( info )
     streaminf *info;
#endif
{
    streaminf *infin;
    assert( info->entrytype == STREAMOUT || info->entrytype == STREAMINOUT || info->entrytype == STREAMSOCKOUT || info->entrytype == STREAMSCR );
    

    if( strchr( info->headername, '|' ) != 0 )         {    /* pipe  */
      open_outpipe( info ); 
}
    else if ( strchr( info->headername, ':' ) != 0 )        /* socket */
{
      open_outsok( info );
}
    else if( strcmp( info->headername, "stdout" ) == 0 )    /* stdout */
{
      openstdout( info ); 
}
    else                                                    /* file */
{
      open_outfile( info );
}
    
    if( info->headfile == (FILE*)0 ){
	info->valid = 0;
	return;  /* failure, but leave them to find out later 
		    if they actually try to do something on this tag */
    }
    
    /* at this point we set the default output format 
       it is "xdr_float" for everything except stdout
       for which it can be set by "data_format=..." on the command line 
       or it is assumed to be the same as on stdin */
    info->format_num = FMT_XDR_FLOAT;
       if( !strcmp(info->tagname,"out") ) {
        if( fetch("data_format","s,",parambuf) )
 	  info->format_num = get_format_num(parambuf);
    }
    
    /* now figure out what the output stream will be
       If the user doesn't remembers to do a hclose before writing to it
       the header will be closed by the first write to this tag.
     */  
    open_outstream( info );  

    /* for tag="out" copy the input header from "in" to the output */
    /* note sepstr_copyh also puts the output datapath is in the header */
    if( !strcmp( info->tagname, "out" )){
	infin = tag_info( "in", TAG_IN );
	sepstr_copyh( infin, info );
    }
    
}


/* write a title to the header file */
/* this happens before the first putch to this file */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void write_title( streaminf *info )
_XFUNCPROTOEND 
#else
void write_title( info )
     streaminf *info;
#endif
{
    char outline[4096];
    char  temp1_ch[4096],temp2_ch[4096];
    int doit;

   sprintf(temp1_ch,"%s.write_title",info->tagname);
   if(0==getch(temp1_ch,"s",&doit)) doit=1;
   if(doit==1){
    
    fputs(maketitle(outline),info->headfile);
    fflush(info->headfile);
    if(ferror(info->headfile)) {
	perror("write_title");
	seperr("write_title: unable to write to output header %s for tag %s",
	       info->headfile, info->tagname );
    }
    }
    info->header_title = 1;
}

/* write the output name in the header */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void write_outname( streaminf* info )
_XFUNCPROTOEND 
#else
void write_outname( info )
     streaminf *info;
#endif
{
    assert( info != 0 );
    assert( info->entrytype == STREAMOUT || info->entrytype == STREAMINOUT || info->entrytype == STREAMSOCKOUT || info->entrytype == STREAMSCR  );
/*    assert( strcmp(info->dataname , 0) );*/

    if( info->header_outname ) return;  /* already done it */

    if( !info->header_title ) write_title( info );

    if( !strcmp( info->dataname, "follow_hdr" )){
        sepstrput(info,"sepstr_ready_out() : sets next in", "s", "stdin" );
    }else{
        char *envname, *envval;
        sepstrput( info, "sets next: in", "s", info->dataname );
        /* write the current value of the environment variable */
        if( (envname = envhead(info->dataname)) != (char*) NULL ) {
            envval = getenv( envname );
            fprintf(info->headfile,
                    "\t\t current environment %s=\"%s\"\n",envname,envval);
        }
    }
    info->header_outname=1;
}



/*   Make sure the output is ready. */
/*   Close the header if the data is following it */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void sepstr_ready_out( streaminf* info )
_XFUNCPROTOEND 
#else
void sepstr_ready_out( info )
     streaminf *info;
#endif
{      


    assert( !info->ready_out );
    assert( info->entrytype == STREAMOUT || info->entrytype == STREAMINOUT || info->entrytype == STREAMSOCKOUT || info->entrytype == STREAMSCR );

    if( !info->header_title ) write_title( info );
    if( !info->header_outname ) write_outname( info );
    
    /* now is the time to write 
     * the output format to the header
     * the user has had their last chance to change it.
     * sreed will use whatever this format is now set to.
     */
    
    sepstrput( info, "data_format", "s", get_format_name(info->format_num) );
    
    if(  strcmp(info->dataname,"follow_hdr") == 0 ){
	syncout( info );
    }

    
    /* we are now ready for writing data */
    info->ready_out = 1;
}      

/*   close header file. */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void sepstr_hclose( streaminf* info )
_XFUNCPROTOEND 
#else
void sepstr_hclose( info )
     streaminf *info;
#endif
{      

    
    if( info->headfile != (FILE*)0 ){

	
	/* make sure it is ready for output */
        if( !info->ready_out ) sepstr_ready_out(info);
	
	if( strcmp(info->dataname,"follow_hdr")!=0 ){/* not following header */
	    assert( info->headfile != (FILE*)0 );
	    fputc( '\n', info->headfile );
	    fclose( info->headfile );
	    info->headfile = (FILE*)0;
	}
	
/*	info->headfile = 0;*/
	
    }else{
	fprintf(stderr,
		"WARNING: Multiple closes of output header with tag \"%s\"\n or close called after first write\n",
		info->tagname);
    }
    
}

/* copy the stored header from "infin" to 
 * the output header stream of "infout"
 */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void sepstr_copyh( streaminf* infin, streaminf* infout )
_XFUNCPROTOEND 
#else
void sepstr_copyh( infin, infout )
     streaminf *infin, *infout;
#endif
{    
    if( infout->header_title ){
	seperr(
	   "cannot copy header to output tag \"%s\" after it has been used \n",
	       infout->tagname);
    }
    
    if( infout->headfile == (FILE*)0 ||  
       isclosed(infout->headfile) ||
       infin->headerbuf == 0 ){
	seperr(" error copying header to output history\n");
    }else{
	fwrite( infin->headerbuf, 1,infin->hdrlen, infout->headfile);
	if ( ferror(infout->headfile) ){
	    perror("sepstr_copyh()");
	    seperr("I/O error writing output history\n");
	} 
        fflush(infout->headfile);
    }

    /*write the output filename in the header in case the
     * program crashes, we don't want the old in= in the header */
    write_outname( infout );

}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static void open_outpipe( streaminf *info )
_XFUNCPROTOEND 
#else
static void open_outpipe( info )
     streaminf *info;
#endif
{
    /* strip the leading pipe symbol */
    if( info->headername[0] != '|' ){
	seperr("output pipe command '%s' doesn't start with | for tag %s\n",
	       info->headername,info->tagname);
    }else{
	info->headfile = popen( &(info->headername)[1], "w" );
	if( info->headfile == (FILE*)0 ) {
	    perror( "sepstrhead:openpipe() ");
	    seperr( "error in opening output pipe %s for tag %s \n",
		   &(info->headername)[1],info->tagname);
	}
	info->isapipe = 1;
    }
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static void  open_outsok( streaminf *info )
_XFUNCPROTOEND 
#else
static void  open_outsok( info )
     streaminf *info;
#endif
{
    char *hostname, *portstr;
    char *tmp;
    int sock,outfd,unixdom;
    
    /* info->headername is a colon separated host:portnumber pair */
    
    tmp = (char *) alloc(1+(int)strlen(info->headername));
    strcpy(tmp,info->headername);
    unixdom = 0;
    
    /* If there is no hostname we just open up a port and wait for 
       someone else to connect */
    
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
        outfd = socklisten( sock, 30 );
    }else{
        int i;
        /* try to connect 30 times, at 1 second intervals */
        for(i =0; i<30;i++) {
          outfd = opensock2( hostname, portstr );
          if( outfd != -1 ) break;
          sleep(1);
        }

    }
    
    if( outfd == -1 ){
        seperr("open_outsok(): open socket \"%s\" failed for tag \"%s\"\n",
	       info->headername,info->tagname);
    }
    
    info->headfile = fdopen( outfd, "w" );
    if( info->headfile == (FILE*)0 ) {
        perror( "sepstrhead:open_outsok() ");
        seperr( "error in opening output socket \"%s\" for tag \"%s\" \n",
	       info->headername,info->tagname);
    }
    
    info->isapipe = 1;
    
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static void  open_outfile( streaminf *info )
_XFUNCPROTOEND 
#else
static void  open_outfile( info )
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
	
	/* create a file with the given name */
	info->headfile = fopen(info->headername,"w+");
	if( info->headfile==(FILE*)0) 
	  info->headfile = fopen(info->headername,"w");
	
	
    }else{ /* header exists */
	
	if(info->entrytype == STREAMOUT ){
	    info->headfile = fopen(info->headername,"w");
	}
    }
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static void openstdout( streaminf *info )
_XFUNCPROTOEND 
#else
static void openstdout( info )
     streaminf *info;
#endif
{
    
    /* The header file attached to stdout (associated with the "out" tag)
     * 
     * seplib allows the par head=... to specify the output headername 
     * stdout and stderr are valid values.
     * the default is equivalent to head=stdout
     */
    info->headfile = stdout; /* default value */
    

    if(getch("head","s",headername) ) { /* head is explicitly specified */
	if ( strcmp(headername,"stdout") == 0 ){
	    info->headername = (char *) alloc((int)strlen(headername)+1); 
	    strcpy(info->headername,headername);
	    info->headfile = stdout; 
            info->isapipe = isapipe( fileno(stdout) );
	}else if(strcmp(headername,"stderr") == 0 )  {
	    info->headername = (char *) alloc((int)strlen(headername)+1); 
	    strcpy(info->headername,headername);
	    info->headfile = stderr;
            info->isapipe = isapipe( fileno(stderr) );
	}else if(((int) strlen(headername)) > 0 ){
	    /* we got something that was a pipe, file or socket */
	    info->headername = (char *) alloc((int)strlen(headername)+1); 
	    strcpy(info->headername,headername);
            if( strchr( info->headername, '|' ) != 0 )            /* pipe  */
	      open_outpipe( info ); 
            else if ( strchr( info->headername, ':' ) != 0 )      /* socket */
	      open_outsok( info );
	    else{
		open_outfile( info );
	    }
	}
    }else{
        info->isapipe = isapipe( fileno(stdout) );
    }
    
}
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void sep_headername( char *tag, char *headername )
_XFUNCPROTOEND 
#else
void sep_headername( char *tag, char *heeadername )
#endif
{
     streaminf *info;
int ierr;

info = tag_info( tag, TAG_INQUIRE );


if( info == 0 )  seperr("requesting headername from tag(%s) where it is undefined \n");

ierr=findnm(fileno(info->headfile),headername,256);
if(ierr<0) seperr("trouble in findnm \n");
else if(ierr==0) strcpy(headername,info->headername);
}

