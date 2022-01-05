
/*
 * Copyright 1988 the Board of Trustees of the Leland Stanford Junior
 * University. Official permission to use this software is included in
 * the documentation. It authorizes you to use this file for any
 * non-commercial purpose, provided that this copyright notice is not
 * removed and that any modifications made to this file are commented
 * and dated in the style of my example below.
 */

/*
 *
 *  source file:   ./remotepenlib/opensock.c
 *
 * Dave Nichols   (SEP), June 09 1988
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 *
 *  Edit History
 *
 *  Original, Dave Nichols (SEP) , June 09, 1988 
 *
 *  Revised, Dave Nichols (SEP) , Sep 12 , 1988 
 *  totally rewritten to do some basic user authentication. Opensock2 now
 *  passes the local hostname , local username and remote username down the
 *  socket so that we can use ruserok to check it.
 *  New entry socklisten added to do the authentication for you.
 *
 *  Revised, Dave Nichols (SEP), Jul. 16, 1990
 *  Added option to timeout the accept, uses alarm to interrrupt it.
 *
 *  Revised, Stew Levin , Feb 25, 1995
 *  Solaris strong-typing for socket call arguments, memcpy instead of
 *  bcopy.
 *
 *  Revised, Stew Levin (MOBIL),  Apr 9, 1997
 *  Use struct sockaddr_un for unix domain sockets.  bind() system
 *  call checks length on Solaris 5.5!  Also corrected unix domain
 *  test in opensock1(). Change portnum from int to char array to
 *  make Unix domain option work.
 *
 *  Revised Bob Clapp  6-2-99 Change to GNU ifdef standard
 */

/* 
 * These routines manage setting up a byte stream 2-way socket between two 
 * processes. The socket is opened in the internet/unix domain as as byte 
 * stream . If the socket is opened in the internet domain the
 * two processes can run on different machines.
 *
 * "opensock1" creates a socket it is passed the port number to open 
 * the returned value is the socket id to be used in the call to accept.
 * If the portnum pointed to is 0 on return the value will be changed to 
 * the actual port number opened. If "unixdom" is not zero the socket is opened
 * in the unix domain otherwise it is opened in the internet domain.
 *
 * "socklisten" is called by the process that called opensock1. It waits 
 * for the second process to bind to the other end and returns a file
 * descriptor for the use of the first process. Returns -1 if it cannot 
 * connect or it times out.
 *
 * "opensock2" called by the second process attatches the second process 
 * to the other end. It is passed the hostname and port number and returns
 * a file descriptor for the use of the second process.
 *
 *   e.g   in the first process
 *     int sock,fd1,secs,unixdom;
 *     char portstr[256];
 *
 *    sock = opensock1(  portstr , unixdom)
 *    ...
 *    now wait around for someone to call 
 *    if the value of secs is >0 it will timeout after that length of time
 *    and return -1.
 *
 *    fd1 = socklisten( sock, secs )
 *
 *
 *   e.g   in the second process
 *	int fd2;
 *      char *portstr;
 * 	char *hostname;
 *
 *    fd2 =  opensock2( hostname, portstr )
 *
 * The socket is opened in the internet or unix domain as as byte stream .
 * hostname should be "unix" for a socket opened in the unix domain.
 *
 */


#include <sepConfig.h>

#include <sys/types.h>

#if defined(HAVE_ERRNO_H) || defined(__APPLE__)
#include <errno.h>
#else

#ifndef STDC_HEADERS
extern int errno;
#endif
#endif

#include <sys/time.h>
#include <sys/ioctl.h>
#include <sys/param.h>
#include <sys/socket.h>

#include <netinet/in.h>
#include <sys/un.h>

#include <arpa/inet.h>

#include <stdio.h>
#include <pwd.h>
#include <netdb.h>
#include <string.h>
#include <unistd.h>
#if defined (HAVE_STDLIB_H) || defined(__APPLE__)
#include<stdlib.h>
#else
extern int atoi();
#endif /* HAVE_STDLIB  */

#include "sep_main_internal.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void strget(int ,char*,int,char*);
void strput( int,char*);
_XFUNCPROTOEND
#else
void strget();
void strput();
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int opensock1( char *portstr, int unixdom )
_XFUNCPROTOEND
#else
int opensock1( portstr, unixdom )
     char *portstr;
     int unixdom;
#endif
{
    int sock;
    socklen_t length;
    /* int len; Not Used */
    /* int buflen; Not Used */
    struct sockaddr_in server_in ;
    struct sockaddr_un server_un ;
    struct sockaddr *server;
    int server_namelen;
    struct linger ling;
    
    /* create socket */
    if( unixdom ){
        sock = socket( AF_UNIX, SOCK_STREAM, 0 );
    }else{
        sock = socket( AF_INET, SOCK_STREAM, 0 );
    }

    if( sock < 0 ) {
	perror("opening stream socket ");
	exit(1);
    }

    ling.l_onoff = 0;
    ling.l_linger = 1;
    setsockopt( sock, SOL_SOCKET, SO_LINGER, (char *) &ling, sizeof(ling) );

    /* Name socket bound to the port */
    if( unixdom ){
        server_un.sun_family = AF_UNIX;
	strcpy(server_un.sun_path,portstr);
	server = (struct sockaddr *) &server_un;
	server_namelen = (int)sizeof(server_un.sun_family)+(int)strlen(server_un.sun_path);
    }else{
        server_in.sin_family = AF_INET;
	server_in.sin_addr.s_addr = INADDR_ANY;
	server_in.sin_port = htons(atoi(portstr)) ;
	server = (struct sockaddr *) &server_in;
	server_namelen = sizeof(server_in);
    }
    if( bind( sock,  server, server_namelen) )  {
	perror("binding stream socket ");
	exit(1);
    }
    
    if (!unixdom) {
	/* find out assigned port number if we passed a zero */
	if( atoi(portstr) == 0 ){
	    length = server_namelen;
	    if( getsockname( sock, server, &length)) {
		perror("getting socket name");
		exit(1);
	    }
	    
	    sprintf(portstr,"%d", ntohs(server_in.sin_port));
	}
    }
    
    /* start accepting connections */
    listen( sock, 5 );
    return(sock);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int opensock2( char* remhost, char* portstr )
_XFUNCPROTOEND
#else
int opensock2( remhost, portstr )
     char *portstr;
     char *remhost;
#endif
{
    int sock;
    /* int buflen; Not Used*/
    struct sockaddr_in server_in;
    struct sockaddr_un server_un;
    struct sockaddr *server;
    int server_namelen;
    struct hostent *hp;
#ifdef AUTH
    char locuser[16],remuser[16],lochost[64];
#endif
    int unixdom;
    struct linger ling;

    unixdom = !strcmp( remhost, "unix" );
    
    /* connect socket using name passed  */
   
    if( unixdom ){
	/* open in the unix domain */ 
        sock = socket( AF_UNIX , SOCK_STREAM, 0 );
    }else{
        hp = gethostbyname(remhost);
        if( hp == 0 ) {
	    fprintf( stderr, "%s: unknown host\n " , remhost );
	    exit(2);
        }
    
        /* create socket */
        sock = socket( AF_INET , SOCK_STREAM, 0 );
    }

    if( sock < 0 ) { perror("opening stream socket ");
	exit(1);
    }
    
    ling.l_onoff = 0;
    ling.l_linger = 1;
    setsockopt( sock, SOL_SOCKET, SO_LINGER, (char *) &ling, sizeof(ling) );
    
    /* connect socket using port number passed  */
    
    if( unixdom ){
        server_un.sun_family = AF_UNIX;
	strcpy(server_un.sun_path,portstr);
	server = (struct sockaddr *) &server_un;
	server_namelen = (int)sizeof(server_un.sun_family)+(int)strlen(portstr);
    }else{
        memcpy( (char *) &server_in.sin_addr, hp->h_addr, hp->h_length);
        server_in.sin_family = AF_INET;
	server_in.sin_port = htons( atoi(portstr) );
	server = (struct sockaddr *) &server_in;
	server_namelen = sizeof(server_in);
    }
    
    if ( connect(sock, server, server_namelen) < 0) {
/*	perror("connecting stream socket");*/
	fprintf(stderr," server not responding on %s \n",remhost);
	return(-1);
    }
#ifdef AUTH
    /* send the local hostname the user name on this machine
       and the user name on the remote machine */
    
    gethostname( lochost, 64);
    strcpy( locuser,  getlogin() );
    strcpy( remuser,  locuser );
    
    strput( sock, lochost );
    strput( sock, locuser );
    strput( sock, remuser );
#endif
    
    return(sock);
}

#include <signal.h>
static int ringring;

#define RETSIGTYPE void

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static RETSIGTYPE catchem( int signo ) { ringring = 0; }
_XFUNCPROTOEND
#else
static RETSIGTYPE catchem(signo) int signo; { ringring = 0; }
#endif/*END OF PROTO TYPING*/
     
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int socklisten(  int sock, int secs )
_XFUNCPROTOEND
#else
int socklisten(  sock, secs )
     int sock;
     int secs;
#endif   
{
    
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
    RETSIGTYPE (*oldsig)(int);
_XFUNCPROTOEND
#else
    RETSIGTYPE (*oldsig)();
#endif

 int omask;

#ifdef AUTH
    char locuser[16], remuser[16],remhost[64];
    struct passwd *pwd;
#endif

#if defined(SOLARIS) || defined(CYGWIN) || defined(LINUX)
    sigset_t newset, oldset;
#endif
       
    int s;
    
#if defined(SOLARIS) || defined(CYGWIN) || defined(LINUX)
    sigemptyset(&newset);
#endif
    /* wait for someone to call us */
    
    ringring = 1;
    
    if( secs > 0 ){  
	oldsig = signal(SIGALRM,catchem); /* clobber old timer trap */
#if defined(SOLARIS) || defined(CYGWIN) || defined(LINUX)
	sigprocmask(SIG_SETMASK, &newset, &oldset);
#else
	omask=sigsetmask(0);
#endif
	alarm(secs);
    }
    
    /* If there is no timer this blocks until someone calls */
     s = accept( sock , 0, 0 );
    
    if( secs > 0 ){  
#if defined(SOLARIS) || defined(CYGWIN) || defined(LINUX)
	sigprocmask(SIG_SETMASK, &oldset, (sigset_t *) NULL);
#else
	sigsetmask(omask);
#endif
	alarm(0); /* shut timer down */
	signal(SIGALRM,oldsig);
    }
    
    if( s == -1 && ringring==1 ) {
	perror("accept");
	return(-1);
    }
    
    if( ringring == 0 ){
	return(-1); /* timed out */
    }
    
     
#ifdef AUTH
    /* check the authorization */
    strget(s,remhost, sizeof(remhost), "remhost");
    strget(s,remuser, sizeof(remuser), "remuser");
    strget(s,locuser, sizeof(locuser), "locuser");
    
    
    /* do they exist in the passwd file */
    setpwent();
    pwd = getpwnam(locuser);
    if (pwd == NULL) {
	fprintf(stderr,"Login incorrect.\n");
	return(-1);
    }
    endpwent();
    
    /* can we chdir to their directory  */
    if (chdir(pwd->pw_dir) < 0) {
	(void) chdir("/");
	fprintf(stderr,"No remote directory.\n");
	return(-1);
    }
    
    /* is the hosts.equiv , .rhosts setup OK */
    if (pwd->pw_passwd != 0 && *pwd->pw_passwd != '\0' &&
	ruserok(remhost, pwd->pw_uid == 0, remuser, locuser) < 0) {
	fprintf(stderr,"Permission denied.\n");
	return(-1);
    }
#endif
    
    return(s);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void strget(int f ,char *buf, int cnt, char *err)
_XFUNCPROTOEND
#else
void strget(f ,buf, cnt, err)
     char *buf;
     int cnt;
     char *err;
     int f;
#endif
{
    char c;
    
    do {
	if (read(f, &c, 1) != 1)
	  exit(1);
	*buf++ = c;
	if (--cnt == 0) {
	    fprintf(stderr,"%s too long\n", err);
	    exit(1);
	}
    } while (c != 0);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void strput( int fdes , char *buf)
_XFUNCPROTOEND
#else
void strput( fdes , buf)
     int fdes;
     char *buf;
#endif
{
  ssize_t rc;
    do { 
	rc = write(  fdes, buf++, 1 );
        if(rc < 1) {
           perror("strput: problem writing string");
           break;
        }
    } while ( *(buf-1) != '\0' );
    return;
}
