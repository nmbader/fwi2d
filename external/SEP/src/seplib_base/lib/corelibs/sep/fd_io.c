/* modules to implement basic I/O functions using raw
   I/O to file descriptors 

   fd_open(), fd_close(), fd_seek(), fd_read(), fd_write(), fd_size()

*/

/* Author Dave Nichols 9-94 
   Modified Dave NIchols 9-1-94  Don't create files with execute permission.
   Modified Dave Nichols 9/8/94 Changes some asserts into calls to seperr
   Modified Bob  Clapp 7/3/96 Changed size to sep_file_size_t 
   Modified Bob  Clapp 7/3/96 Changed seek to sep_file_size_t (and offset
                              to sep_off_t) 
   Modified Bob  Clapp 7/4/96 Changed curpos to sep_file_size_t  
   Modified Stew Levin 4/10/97 Use ssize_t and size_t for read/write
   Modified Robert  Clapp 7/20/97 Changed STDC to XFUNCPROTO  
   Modified Robert  Clapp 12/17/97 Added STREAMSOCKOUT options
   Modified Robert  Clapp 10/6/98 replaced tell call because now non-standard
   Modified Robert  Clapp 6/1/99 Introduce GNU defs, include possixstat.h
 */
#include <sepConfig.h>
#if defined (HAVE_FCNTL_H) || defined(__APPLE__)
#include <fcntl.h>
#endif

#include <unistd.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "streamlist.h"
#include "sep_main_internal.h"
#include <sepcube.h>
#ifdef _POSIX_SOURCE
#include <unistd.h>
#else /* not posix source */
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

#include "../include/posixstat.h"

struct fd_info {
    int fd;
    sep_file_size_t curpos;
};

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void fd_open( streaminf*, void**);
void fd_close( streaminf*, void*);
ssize_t fd_read(streaminf*,void*,void*,size_t);
ssize_t fd_write(streaminf*,void*,void*,size_t);
sep_file_size_t fd_seek( streaminf*,void*, sep_off_t, int);
sep_file_size_t fd_size( streaminf*,void*);
_XFUNCPROTOEND
#else
void fd_open();
void fd_close();
ssize_t fd_read();
ssize_t fd_write();
sep_file_size_t fd_seek();
sep_file_size_t fd_size();
#endif


/* FROM GNU BASH */

#if !defined (S_IRWXU)
#  if !defined (S_IREAD)
#    define S_IREAD 00400
#    define S_IWRITE  00200
#    define S_IEXEC 00100
#  endif /* S_IREAD */
#endif




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void fd_open( streaminf *info , void** ioinfo )
_XFUNCPROTOEND
#else
void fd_open( info , ioinfo )
     streaminf *info;
     void **ioinfo;
#endif
{
    struct fd_info *inf;
    char *exp_name;

    /* make an info structure for FILE* I/O and set it in the streaminf */
    inf = (struct fd_info*) malloc( sizeof(struct fd_info) );
    *ioinfo=inf;
    inf->curpos=0;
    inf->fd=-1;

    
    /* open the dataset */
    switch( info->entrytype ){
	
	case( STREAMIN ):
	  if( !strcmp(info->dataname,"follow_hdr") ||  
	     !strcmp(info->dataname,"stdin") ){
	      seperr( "%s is invalid for file descriptor I/O for tag \"%s\"\n",
		     info->dataname,info->tagname);
	  }else{
	      /* open up a file */
	      exp_name = expand_info(info->dataname,info);
				if(0!=strcmp(exp_name,"-1")){
	      inf->fd= open(exp_name,O_RDONLY);
	      if(inf->fd == -1 ){
		  perror("fd_open");
		  fprintf(stderr,
			"unable to obtain input file \"%s\" for tag \"%s\"\n",
			 exp_name,info->tagname);
		  info->valid = 0; return;
	      }
	      if( fdordinary( inf->fd ) ) info->backseekable=1;
			}
	      free(exp_name);
	  }
	break;
	
	case( STREAMOUT ):
	  if( !strcmp( info->dataname,"follow_hdr") || 
	     !strcmp( info->dataname,"stdout") ){
	      seperr( "%s is invalid for file descriptor I/O for tag \"%s\"\n",
		     info->dataname,info->tagname);
	  } else {
	      exp_name = expand_info(info->dataname,0);
   inf->fd= open(exp_name,O_WRONLY|O_CREAT|O_TRUNC,S_IRUSR|S_IRGRP|S_IROTH|S_IWUSR|S_IWGRP|S_IWOTH);
	      if(inf->fd  == -1 ){
		  perror("fd_open");
		  fprintf(stderr, 
			" unable to open output datafile %s for tag %s \n",
			 info->dataname,info->tagname);
		  info->valid = 0; return;
	      }
	      if( fdordinary( inf->fd ) ) info->backseekable=1;
	      free(exp_name);
	  }
	  break;
	  
	case( STREAMINOUT ):
	  exp_name = expand_info(info->dataname,info);
	  if( info->headerbuf != 0 ){
	    inf->fd = open(exp_name,O_RDWR );
          }
	  if(inf->fd == -1 ){
	      inf->fd=open(exp_name,O_RDWR|O_CREAT|O_TRUNC,S_IRUSR|S_IRGRP|S_IROTH|S_IWUSR|S_IWGRP|S_IWOTH);
	  }
	
	  if( inf->fd  == -1 ){
	      perror("fd_open");
	      fprintf(stderr,
		"unable to open input/output datafile %s for tag %s \n",
		     info->dataname,info->tagname);
	      info->valid = 0; return;
	  }
	  if( fdordinary( inf->fd ) ) info->backseekable=1;
	  free(exp_name);
	  break;

	case( STREAMSCR ):
	  exp_name = expand_info(info->dataname,info);
	      inf->fd=open(exp_name,O_RDWR|O_CREAT|O_TRUNC,S_IRUSR|S_IRGRP|S_IROTH|S_IWUSR|S_IWGRP|S_IWOTH);
	
	  if( inf->fd  == -1 ){
	      perror("fd_open");
	      fprintf(stderr,
		"unable to open input/output datafile %s for tag %s \n",
		     info->dataname,info->tagname);
	      info->valid = 0; return;
	  }
	  if( fdordinary( inf->fd ) ) info->backseekable=1;
	  free(exp_name);
	  break;



	case( STREAMSOCKOUT ):
		/*I THINK THAT YOU CAN'T DO THE SOCKET TRICK WHEN USING FD IO */
	  if( !strcmp(info->dataname,"follow_hdr") ||  
	     !strcmp(info->dataname,"stdin") ){
	      seperr( "%s is invalid for file descriptor I/O for tag \"%s\"\n",
		     info->dataname,info->tagname);
    }
	}
    
    info->streamfd = inf->fd;

    switch( info->entrytype ){
	case( STREAMIN ):
    		info->streamfile = fdopen(inf->fd,"r");
	break;
	case( STREAMOUT ):
    		info->streamfile = fdopen(inf->fd,"w");
        break;
	case( STREAMINOUT ):
    		info->streamfile = fdopen(inf->fd,"r+");
        break;
    }
    
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void fd_close( streaminf *info, void* ioinfo )
_XFUNCPROTOEND
#else
     void fd_close( info , ioinfo )
     streaminf *info;
     void* ioinfo;
#endif
{
    struct fd_info *inf;

    if(ioinfo==NULL) return;
    
    inf=(struct fd_info *)ioinfo;
    
    if( inf->fd != -1  ) close( inf->fd );
    
    free(inf);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
ssize_t fd_read( streaminf *info, void* ioinfo, void *buffer, size_t nbytes )
_XFUNCPROTOEND
#else
ssize_t fd_read( info, ioinfo,  buffer, nbytes )
     streaminf *info;
     void* ioinfo;
     void* buffer;
     size_t nbytes;
#endif
{
    ssize_t total, nread;
    struct fd_info *inf;
    size_t sep_bufsiz ;
    char *cptr;
    int infd;
//fprintf(stderr,"CHECK READ %d  \n",(int)nbytes);
    
    cptr=(char*)buffer;
    
    if(!info->isapipe){
	sep_bufsiz = SEP_BUFSIZ;
    }else{
	sep_bufsiz = 1024;
    }
    
    inf=(struct fd_info *)ioinfo;
    if( inf->fd == -1 ) seperr("Invalid dataset \"%s\" for tag \"%s\"\n",
				info->dataname,info->tagname);
    
    infd = inf->fd;
    

    total=0;
 
    do{
	nread=read(infd,cptr, MIN(nbytes-total,sep_bufsiz));
	if( nread == 0 ) break;
	
	total += nread;
	cptr += nread;
	
	if (nread<0) {
	    perror("fd_read()");
	    seperr ("fd_read() : I/O error, tag \"%s\"\n",info->tagname);            }
    }while( total <nbytes );
    
    inf->curpos += (sep_file_size_t)total;
    
    return total;
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
ssize_t fd_write( streaminf *info, void* ioinfo, void *buffer, size_t nbytes )
_XFUNCPROTOEND
#else
ssize_t fd_write( info, ioinfo, buffer, nbytes )
     streaminf *info;
     void* ioinfo;
     void* buffer;
     size_t nbytes;
#endif
{
    ssize_t total, nwrite;
    struct fd_info *inf;
    size_t sep_bufsiz ;
    char *cptr;
    int outfd;


    
//fprintf(stderr,"CHECK WRITE %d  \n",(int)nbytes);
    cptr=(char*)buffer;
    
    if(!info->isapipe){
	sep_bufsiz = SEP_BUFSIZ;
    }else{
	sep_bufsiz = 1024;
    }
    
    inf=(struct fd_info *)ioinfo;
    if( inf->fd == -1 ) seperr("Invalid dataset \"%s\" for tag \"%s\"\n",
				info->dataname,info->tagname);
    
    outfd = inf->fd;
    
    total=0;
    for (total=0; (total < nbytes) &&
	 (0 != (nwrite=write(outfd,cptr+total,MIN(nbytes-total,sep_bufsiz))));
	 total += nwrite) {
	
	if (nwrite<0) {
	    perror("info_rite()");
	    seperr ("info_rite() : I/O error, tag \"%s\"\n",info->tagname);
	}
    }
    
    inf->curpos += (sep_file_size_t)total;
    
    return total;
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
     sep_file_size_t fd_seek( streaminf *info, void* ioinfo, sep_off_t offset, int whence )
_XFUNCPROTOEND
#else
     sep_file_size_t fd_seek( info, ioinfo, offset, whence )
     streaminf *info;
     void* ioinfo;
     int whence;
     sep_off_t offset;
#endif
{
    struct fd_info *inf;
/*    long seekpos;*/
    off_t seekpos;
    int ierr;
    
//fprintf(stderr,"CHECK SEEK %d %d \n",(int)offset,(int)whence);
    inf=(struct fd_info *)ioinfo;
    
    ierr=0;
    
    switch( whence ){
      case SEEK_SET:
	seekpos = (off_t) offset ;
	if( lseek(inf->fd,seekpos,whence)< 0){
	    perror("fd_seek");
	    ierr=1;
	}else{
	    inf->curpos =(sep_file_size_t) offset;
	}
        break;
      case SEEK_CUR:
	seekpos = (off_t) offset ;
	if( lseek(inf->fd,seekpos,whence) < 0 ){
	    perror("fd_seek");
	    ierr=1;
	}else{
	    inf->curpos +=(sep_file_size_t) offset;
	}
	break;
      case SEEK_END:
	seekpos = (off_t) offset ;
	if( lseek(inf->fd,seekpos,whence) < 0 ){
	    perror("fd_seek");
	    ierr=1;
	}else{
			/*this really ugly and annoying, but I can't think of another
       way around tell now being non-standard (linux/gcc*/
	    inf->curpos =(sep_file_size_t) (fsize(inf->fd))  - 
                  (sep_file_size_t) (seekpos);
	}
	break;
      default:
	fprintf( stderr,"unknown seek type requested\n");
	break;
    }
    
    if( ierr == 0 ){
	return(inf->curpos);
    }else{
	return(-1);
    }
    
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
sep_file_size_t fd_size( streaminf *info, void* ioinfo )
_XFUNCPROTOEND
#else
sep_file_size_t fd_size( info, ioinfo )
streaminf *info;
void* ioinfo ;
#endif
{
struct fd_info *inf;
sep_file_size_t size;

inf=(struct fd_info *)ioinfo;
assert( inf != 0 );

if( inf->fd == -1 ) seperr("Invalid dataset \"%s\" for tag \"%s\"\n",
				info->dataname,info->tagname);

	size = (sep_file_size_t)fsize( inf->fd ) ;
return size;
}



/* Set up the function pointers in an info structure to point
 * at the fd I/O functions
 */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void init_fd_io( streaminf* info )
_XFUNCPROTOEND
#else
     void init_fd_io( info )
     streaminf * info;
#endif
{
    info->open_func = fd_open;
    info->close_func = fd_close;
    info->read_func = fd_read;
    info->write_func = fd_write;
    info->seek_func = fd_seek;
    info->size_func = fd_size;
}

