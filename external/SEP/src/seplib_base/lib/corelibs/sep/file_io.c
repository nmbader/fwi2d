/* modules to implement basic I/O functions using buffered
   I/O to FILE* 

 file_open(), file_close(), file_seek(), file_read(), file_write(), file_size()

*/

/* Author Dave Nichols 9-94
   Modified Dave Nichols 9/8/94 Changes some asserts into calls to seperr
   Modified by David Lumley 9-29-94, commented out some seek error prints
   Modified by David Lumley 10-01-94, modified backwards seek logic, 
	 "if( offset >" --> "if( offset >=", converted errors to sepperr().
   Modified Biondo 7/3/96 Changed size to sep_file_size_t 
   Modified Biondo 7/3/96 Changed seek to sep_file_size_t (and offset
                              to sep_off_t) 
   Modified Biondo 7/4/96 Changed curpos to sep_file_size_t  
   Modified Stewart A. Levin 4/10/97 Use ssize_t and size_t
   Modified Robert G. Clapp 7/20/97 Changed STDC to ANSI_DECL
   Modified Robert G. Clapp 12/18/97 Added STREAMSOCKOUT
  
 */
#include <sepConfig.h>
#include <fcntl.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <sepcube.h>
#include "sep_main_internal.h"
#include "streamlist.h"

#ifdef _POSIX_SOURCE
#include <unistd.h>
#else /* not posix source */
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif


struct file_info {
      FILE* file;
      sep_file_size_t curpos;
      off_t seekoffset;
};


/*INTERNAL PROTOTYPES */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void file_open( streaminf*, void**  );
void file_close( streaminf*, void * );
ssize_t file_read( streaminf *, void* , void *, size_t  );
ssize_t file_write( streaminf *, void* , void *, size_t  );
static sep_off_t discard(streaminf*, FILE*, sep_off_t);
sep_file_size_t file_seek( streaminf *, void* , sep_off_t , int );
static sep_off_t discard( streaminf* , FILE* disc, sep_off_t );
sep_file_size_t file_size( streaminf *, void*  );
_XFUNCPROTOEND
#else
void file_open();
void file_close();
ssize_t file_read();
ssize_t file_write();
static sep_off_t discard();
sep_file_size_t file_seek();
static sep_off_t discard();
sep_file_size_t file_size();
#endif


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void file_open( streaminf *info, void** ioinfo )
_XFUNCPROTOEND
#else
void file_open( info, ioinfo )
streaminf *info;
void** ioinfo;
#endif
{
struct file_info *inf;
char *exp_name;

/* make an info structure for FILE* I/O and set it in the streaminf */
inf = (struct file_info*) malloc( sizeof(struct file_info) );
*ioinfo = inf;
inf->curpos=0;


/* open the dataset */
switch( info->entrytype ){

    case( STREAMIN ):
	if( !strcmp(info->dataname,"follow_hdr") ||  
	    !strcmp(info->dataname,"stdout") ||
	    !strcmp(info->dataname,"stdin") ){

	    if( !strcmp(info->dataname,"stdin") ){
	      inf->file = stdin;
	      inf->file = info->headfile;
	    }else{
	      inf->file = info->headfile;
	    }

	    info->headfile = 0; /* we can't read the input header any more */
	    info->isapipe = isapipe( fileno(inf->file) );
	    /* if it is a regular file we can seek with a seek offset */
	    if( !fdordinary( fileno(inf->file) ) ){
	       info->backseekable=  0;
	       inf->seekoffset = 0;
	    }else{
	       info->backseekable= 1;
	       inf->seekoffset = ftell(inf->file);
	    }
	}else{
	    /* open up a file */
	    exp_name = expand_info(info->dataname,info);
			if(0!=strcmp(exp_name,"-1")){
	    inf->file= fopen(exp_name,"r+");
            if(inf->file == 0 ){
              inf->file= fopen(exp_name,"r");
            }
            if(inf->file == 0 ){
                perror("file_open");
                seperr( "unable to obtain input file \"%s\" for tag \"%s\"\n",
                info->dataname,info->tagname);
            }
            if( fdordinary( fileno(inf->file) ) ) info->backseekable=1;
					}
	    free(exp_name);
	}
      break;

    case( STREAMOUT ):
	if( !strcmp( info->dataname,"follow_hdr") ){
	    inf->file = info->headfile;
	    if( !fdordinary( fileno(inf->file) ) ){
	       info->backseekable=  0; inf->seekoffset = 0;
	    }else{
	       info->backseekable = 1; inf->seekoffset = ftell(inf->file);
	    }
	}else if( !strcmp( info->dataname,"stdout") ){
	    inf->file = stdout;
	    info->isapipe = isapipe( fileno(stdout) );
	    if( !fdordinary( fileno(inf->file) ) ){
	       info->backseekable=  0; inf->seekoffset = 0;
	    }else{
	       info->backseekable = 1; inf->seekoffset = ftell(inf->file);
	    }
	}else{
	   exp_name = expand_info(info->dataname,0);
	   inf->file= fopen(exp_name,"w+");
            if(inf->file == 0 ){
              inf->file= fopen(exp_name,"w");
	    }

            if(inf->file  == 0 ){
                perror("file_open");
                seperr( " unable to open output datafile %s for tag %s \n",
                     info->dataname,info->tagname);
            }
            if( fdordinary( fileno(inf->file) ) ) info->backseekable=1;
	    free(exp_name);
        }
      break;

    case( STREAMINOUT ):
	exp_name = expand_info(info->dataname,info);
	if( info->headerbuf == 0 ){
          inf->file= fopen(exp_name,"w+"); /* new header file */
        }else{
				inf->file= fopen(exp_name,"r+");
        }
        if(inf->file == 0 ){
            inf->file=fopen(exp_name,"w+");
	}

        if( inf->file  == 0 ){
            perror("file_open");
            seperr("unable to open input/output datafile %s for tag %s \n",
                 info->dataname,info->tagname);
        }
        if( fdordinary( fileno(inf->file) ) ) info->backseekable=1;
	free(exp_name);
      break;

		
		 case( STREAMSCR ):
  exp_name = expand_info(info->dataname,info);
          inf->file= fopen(exp_name,"w+"); /* new header file */

        if( inf->file  == 0 ){
            perror("file_open");
            seperr("unable to open scratch datafile %s for tag %s \n",
                 info->dataname,info->tagname);
        }
        if( fdordinary( fileno(inf->file) ) ) info->backseekable=1;
  free(exp_name);
      break;


    case( STREAMSOCKOUT ):
	 if( !strcmp( info->dataname,"follow_hdr") ){
      inf->file = info->headfile;
      if( !fdordinary( fileno(inf->file) ) ){
         info->backseekable=  0; inf->seekoffset = 0;
      }else{
         info->backseekable = 1; inf->seekoffset = ftell(inf->file);
      }
  }else{
			exp_name = expand_info(info->dataname,info);
     exp_name = expand_info(info->dataname,0);
     inf->file= fopen(exp_name,"w+");
            if(inf->file == 0 ){
              inf->file= fopen(exp_name,"w");
      }
  if(inf->file  == 0 ){
                perror("file_open");
                seperr( " unable to open output datafile %s for tag %s \n",
                     info->dataname,info->tagname);
            }
            if( fdordinary( fileno(inf->file) ) ) info->backseekable=1;
      free(exp_name);
        }
      break;
}

info->streamfile = inf->file;
info->streamfd = fileno(inf->file);

}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void file_close( streaminf *info, void *ioinfo )
_XFUNCPROTOEND
#else
void file_close( info, ioinfo )
streaminf *info;
void* ioinfo;
#endif
{
struct file_info *inf;

inf=(struct file_info *)ioinfo;
assert( inf != 0 );

if( inf->file != 0  ) fclose( inf->file );

free(inf);

}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
ssize_t file_read( streaminf *info, void* ioinfo, void *buffer, size_t nbytes )
_XFUNCPROTOEND
#else
ssize_t file_read( info, ioinfo, buffer, nbytes )
streaminf *info;
void* ioinfo;
void* buffer;
size_t nbytes;
#endif
{
ssize_t total, nread;
struct file_info *inf;
size_t sep_bufsiz ;
char *cptr;


cptr=(char*)buffer;

if(!info->isapipe){
   sep_bufsiz = 65536;
}else{
   sep_bufsiz = 1024;
}

inf=(struct file_info *)ioinfo;
if( inf->file == 0 ) seperr("Invalid dataset \"%s\" for tag \"%s\"\n",
                                info->dataname,info->tagname);


total=0;



do{
    nread=(long)fread(cptr,1,MIN(nbytes-total,sep_bufsiz),inf->file);


    total += nread;
    cptr += nread;
    if ( ferror(inf->file) ) {
        perror("file_read()");
        seperr ("file_read() : I/O error, tag \"%s\"\n",info->tagname);
    }

    if( nread == 0 ) break;

}while( total < nbytes );

inf->curpos += (sep_file_size_t)total;


return total;
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
ssize_t file_write( streaminf *info, void* ioinfo, void *buffer, size_t nbytes )
_XFUNCPROTOEND
#else
ssize_t file_write( info, ioinfo, buffer, nbytes )
streaminf *info;
void* ioinfo;
void* buffer;
size_t nbytes;
#endif
{
ssize_t total, nwrite;
struct file_info *inf;
size_t sep_bufsiz ;
char *cptr;

cptr=(char*)buffer;

if(!info->isapipe){
   sep_bufsiz = 65536;
}else{
   sep_bufsiz = 1024;
}

inf=(struct file_info *)ioinfo;
if( inf->file == 0 ) seperr("Invalid dataset \"%s\" for tag \"%s\"\n",
                                info->dataname,info->tagname);

total=0;

do{
    nwrite=(long)fwrite(cptr,1,MIN(nbytes-total,sep_bufsiz),inf->file);

    total += nwrite;
    cptr += nwrite;

    if ( ferror(inf->file) ) {
        perror("file_write()");
        seperr ("file_write() : I/O error, tag \"%s\"\n",info->tagname);
    }

    if ( feof(inf->file) ) {
        fprintf(stderr,"unexpected eof for tag %s\n",info->tagname);
    }
    if( nwrite == 0 ) break;

}while( total < nbytes );

fflush( inf->file );

inf->curpos += (sep_file_size_t)total;

return total;
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static sep_off_t discard(streaminf*, FILE*, sep_off_t); 
_XFUNCPROTOEND
#else
static sep_off_t discard();
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
sep_file_size_t file_seek( streaminf *info, void* ioinfo, sep_off_t offset, int whence )
_XFUNCPROTOEND
#else
sep_file_size_t file_seek( info, ioinfo, offset, whence )
streaminf *info;
void* ioinfo;
int whence;
sep_off_t offset;
#endif
{
struct file_info *inf;
off_t seekpos;
int ierr;

inf=(struct file_info *)ioinfo;
assert( inf != 0 );
if( inf->file == 0 ) seperr("Invalid dataset \"%s\" for tag \"%s\"\n",
                                info->dataname,info->tagname);

ierr=0;

if( !info->backseekable ){

    switch( whence ){

      case SEEK_SET:
        if( (sep_file_size_t)offset > inf->curpos ){
          inf->curpos += (sep_file_size_t)discard( info, inf->file, (sep_off_t)((sep_file_size_t)offset-inf->curpos));
        }else if( (sep_file_size_t)offset < inf->curpos ){
            fprintf(stderr,"backwards absolute seek on stream device \n");  
            ierr=1;
        }
        break;
      case SEEK_CUR:
        if( offset > 0 ){
            inf->curpos += (sep_file_size_t)discard( info, inf->file, offset ); 
        }else if( offset < 0 ){
            fprintf( stderr,"backwards relative seek on stream device \n");
            ierr=1;
        }
        break;
      case SEEK_END:
        fprintf( stderr,"seek relative to end on stream device \n");
        ierr =1;
        break;
      default:
        fprintf( stderr,"unknown seek type requested\n");
        ierr =1;
        break;
    }

}else{
    /* a regular backseekable file !  */
    /* we have to take care of a possible seekoffset */
    switch( whence ){
      case SEEK_SET:
        seekpos = ((off_t)offset + inf->seekoffset);
	if( !fseek(inf->file,seekpos,whence) ){
	   inf->curpos = (sep_file_size_t)offset;
        }else{
	   perror("file_seek");
	   ierr=1;
        }
        break;
      case SEEK_CUR:
        seekpos = (off_t)offset ;
	if( !fseek(inf->file,seekpos,whence) ){
	   inf->curpos += (sep_file_size_t)offset;
        }else{
	   perror("file_seek");
	   ierr=1;
        }
        break;
      case SEEK_END:
        seekpos = (off_t)offset ;
	if( !fseek(inf->file,seekpos,whence) ){
	   inf->curpos = (sep_file_size_t)ftell( inf->file );
        }else{
	   perror("file_seek");
	   ierr=1;
        }
        break;
      default:
        fprintf( stderr,"unknown seek type requested\n");
        break;
    }
}

if( ierr == 0 ){
    return(inf->curpos);
}else{
    return(-1);
}

}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
sep_off_t discard( streaminf* info, FILE* disc, sep_off_t nbytes)
_XFUNCPROTOEND
#else
sep_off_t discard( info, disc, nbytes )
streaminf *info;
FILE* disc;
sep_off_t nbytes;
#endif
{
static char     skip[BUFSIZ];
sep_off_t             total, len, nread ;

    total = 0;

    do {
       len = MIN( BUFSIZ, nbytes - total );
       nread = (long)fread( skip, 1, len, disc );
       if( nread < 0 ) {
        seperr("sseek(): I/O error in seek for tag \"%s\"\n",info->tagname);
       }
       if( nread == 0 ) break;
       total += nread;
    }while(total<nbytes);

/*    return ( nbytes - total );*/
    return (  total );
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
sep_file_size_t file_size( streaminf *info, void* ioinfo )
_XFUNCPROTOEND
#else
sep_file_size_t file_size( info, ioinfo )
streaminf *info;
void* ioinfo ;
#endif
{
struct file_info *inf;
sep_file_size_t size;

inf=(struct file_info *)ioinfo;
assert( inf != 0 );

if( info->isapipe ) return 0;

if( inf->file == 0 ) seperr("Invalid dataset \"%s\" for tag \"%s\"\n",
                                info->dataname,info->tagname);
size = (sep_file_size_t)fsize( fileno( inf->file )) -(off_t)inf->seekoffset ;

return size;

}


/* Set up the function pointers in an info structure to point
 * at the FILE* I/O functions
 */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void init_file_io( streaminf* info )
_XFUNCPROTOEND
#else
void init_file_io( info )
streaminf* info;
#endif
{
info->open_func = file_open;
info->close_func = file_close;
info->read_func = file_read;
info->write_func = file_write;
info->seek_func = file_seek;
info->size_func = file_size;
}
