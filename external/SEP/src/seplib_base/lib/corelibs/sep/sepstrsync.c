/*
 * syncin( streaminf *info)
 * syncout( streaminf *info)
 *
 * Routines to perform synchronization of piped headers and data
 * in seplib. If the header and data are on the same stream we need to mark
 * the end of the header and the start of the data. The end of the data is
 * marked by a ^D. (EOT)
 *
 *
 * There are three synchronization methods:
 *
 * (the old method)
 *    The sending program writes its 9digit PID and a ^D and goes to sleep.
 *    The receiving program sends a SIGALRM to the process to let it know
 *    it has finished reading the header and it can start writing data.
 *    Advantages: fast (raw reads).
 *    Disadvantages: only works if jobs are on the same machine, much of our
 *    I/O is now fully buffered so mixing raw and buffered is dangerous.
 *
 *
 * ( the pipe method )
 *    The sending program writes "PIPE" followed by its hostname followed 
 *    by a port number followed by its 9 digit PID and a ^D and goes into
 *    accept().
 *    The receiving program connects to the port and sends a "GOTIT" message
 *    the sending program sends back an "ACK".
 *    Advantages: works across network, each job is left with a socket that
 *    could be used for interprocess communication in the future.
 *    Disadvantages: More complicated, requires berkley type sockets.
 *
 * ( the simple method )
 *   The sending program write EOL EOL EOT.
 *   The receiving program reads until it gets an EOT, if the two previous
 *   characters are EOL it is done.
 *   Advantages: simple, can be used to write header+data directly to any
 *   unix program (which doesn't need to understand a more complicated
 *   protocol.) e.g. write to tape or to a single file on disk.
 *   Disadvantages: Requires character at a time I/O, this realistically means
 *   buffered I/O (using FILE*) any raw I/O after this could be fatal.
 *
 *
 * Author: Dave Nichols (SEP)
 * Revised: Stew Levin (MOBIL) April 9, 1997
 *         Changed portnum to char string to match opensock changes
 *         supporting Unix domain sockets.
 * Revised: Robert Clapp 7/19/97 Prototypes
 */
#include <sepConfig.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#if defined(__APPLE__) || defined(LINUX)
#define USE_SOCKETS
#endif
#include <sepcube.h>
#include "sep_main_internal.h"
#include "streamlist.h"

#include <signal.h>
#include <sys/types.h>

#define EOT 004
#define EOL 014

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void syncin( streaminf *info )
_XFUNCPROTOEND
#else
void syncin( info )
streaminf *info;
#endif
{ 
    char* digits;
    int i,flen,prevpid;
    char remhost[255],msgbuf[10];
    char portstr[256];
    ssize_t rc1, rc2;
    flen= info->hdrlen;

    /* First check for previous two characters ctrl-L */
    if( info->headerbuf[flen-1] == EOL && info->headerbuf[flen-2] == EOL ){
	/* we are done, trim off the ctrl-l and return */
	(info->hdrlen) -=2;  info->headerbuf[info->hdrlen] = '\0';
	return;
    }

    if(flen > 9){
	digits = info->headerbuf+flen-1; i=9;
	while(i){
	    if(!isdigit(*digits)) break;
	    digits--; i--;
	}
	if(!i && (*digits) == ' '){
	    if(1 == sscanf(digits," %9d",&prevpid)) {
		/* found pid of sender */
		flen = digits - info->headerbuf + 1;

		
		/* now look backwards for PIPE synch */
		/* arbitrary max of 300 characters */
		i=MIN(flen,300);

		while(i != 0 ) {
		    if( *digits=='P' && !strncmp(digits,"PIPE",4))break;
		    digits--; i--;
		}

		if( i == 0 ){
		    /* didn't find the PIPE tag must be using old protocol */
		    kill(prevpid,SIGALRM);

		}else{
		    /* found the tag, read hostname and port number */
		    if( sscanf(digits,"PIPE %s %s ", remhost,portstr) == 2){
			/* connect to the sender */
			info->sockfd =  opensock2( remhost,portstr);
			if( info->sockfd != -1 ){
			    rc1 = write(info->sockfd,"GOTIT",6);
			    rc2 = read(info->sockfd,msgbuf,4);
			    if(strcmp(msgbuf,"ACK") || rc1 != 6 || rc2 != 4) {
			      seperr("syncin() sync err\n" );
                            }
			}
		    }
		    
		    flen = digits - info->headerbuf ;
		}
	    }
	}
    }
    info->hdrlen=flen;
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void syncout( streaminf *info )
_XFUNCPROTOEND
#else
void syncout( info )
streaminf *info;
#endif
{
    int sock;
   
    char mesgbuf[10];
    char hostnm[255];
    char portstr[256];
    int isoutstream,timeout;
    ssize_t rc1, rc2;

    /* If data and header are going down the same route but it is
     * an ordinary file we will just separate them with a ctrl-L ctrl-L ctrl-D
     *
     * We can also force this with isoutstream=y on the command line.
     */
  
    isoutstream=0; getch("isoutstream","1",&isoutstream); 


    if( !info->isapipe || isoutstream ){
	fprintf(info->headfile,"\n%c%c%c",EOL,EOL,EOT);
        fflush(info->headfile );
		
    }else{
     

    portstr[0] = '\0';
    sock=opensock1(portstr, 0);
    if( gethostname(hostnm,255)!=0 )
      seperr("syncout(): getting hostname\n");
    
    fprintf(info->headfile,"PIPE %s %s %09d%c",hostnm,portstr,getpid(),EOT);
    fflush(info->headfile );
    
    /* wait to be called back, or time out after "timeout" secs, default 5 */
    timeout=5; getch("timeout","i",&timeout);


   
    if( (info->sockfd=socklisten(sock, timeout )) != -1 ){
	rc1 = read( info->sockfd, mesgbuf, 6 );
        if(rc1 != 6) perror("syncout()");
	if( strcmp( mesgbuf, "GOTIT" ) )
	  seperr("syncout(): pipe synch failed!\n");
	rc2 = write( info->sockfd, "ACK", 4);
        if(rc2 != 4) {
          perror("syncout()");
	  seperr("syncout(): pipe synch write failed!\n");
        }
    }
    }
    
}
