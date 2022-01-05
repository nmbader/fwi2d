/* modules to implement basic I/O functions using raw
   I/O to multiple file descriptors

   multifd_open(), multifd_close(), multifd_seek(), multifd_read(),
   multifd_write(), multifd_size()

*/

/* Author Dave Nichols 9-94
 * Revised: stew 2-25-95 Solaris hacks
 * Revised: Ray Abma 30 May 95
 *  Added Dave's changes that stew missed -
 *      Modified Dave Nichols 21/12/94
 *        Make the seek function return -1 if you try to seek off the end.
 *        Correctly update buffer pointer in read and write functions.
 * Revised: Stewart A. Levin (MOBIL) 4-10-97 
 *      Changed to use and support off_t, sep_off_t and sep_file_size_t 
 *      correctly.
 * Revised: Bob Added Prototypes 7-97
 * Revised: Bob:Added error conditions on STREAMSOCKOUT
 * Revised:  10/7/98 Bob: Got rid of tell because not standard (linux problem)
 * Revised:  6/2/98 Bob:  Added GNU ifedfdefs
 */

#include <sepConfig.h>
#if defined(HAVE_FCNTL_H) || defined(__APPLE__)
#include <fcntl.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <string.h>


#if defined(HAVE_ERRNO_H) || defined(__APPLE__)
#include <errno.h>
#endif

#ifndef STDC_HEADERS
extern int errno;
#endif

#include <unistd.h>
#include "posixstat.h"


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

struct fd_info {
    int fd;
    sep_off_t curpos;
    sep_off_t size;
};

struct multifd_info {
	int nfile;
	int active;
	struct fd_info* infos[MAX_SEP_FILE];
  char* name[MAX_SEP_FILE];
	int used[MAX_SEP_FILE];
	int last_used;
};

/*PROTOTYPES FOR FUNCTIONS USED ONLY INSIDE THIS FILE */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void multifd_open( streaminf * , void**  );
void multifd_close( streaminf *, void*  );
static ssize_t onefd_read( int , void *, size_t  );
ssize_t multifd_read( streaminf *, void* , void *, size_t  );
static ssize_t onefd_write( int , void *, size_t  );
ssize_t multifd_write( streaminf *, void* , void *, size_t  );
sep_file_size_t multifd_seek( streaminf *, void* , sep_off_t , int);
sep_file_size_t multifd_size( streaminf *, void*  );
void update_out(streaminf*, void*);
_XFUNCPROTOEND
#else
void multifd_open();
void multifd_close();
static ssize_t onefd_read();
ssize_t multifd_read();
static ssize_t onefd_write();
ssize_t multifd_write();
sep_file_size_t multifd_seek();
sep_file_size_t multifd_size();
void update_out();
#endif


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void multifd_open( streaminf *info , void** ioinfo )
_XFUNCPROTOEND
#else
void multifd_open( info , ioinfo )
     streaminf *info;
     void **ioinfo;
#endif
{
    struct multifd_info *inf;
    struct fd_info *one_inf;
    char *all_names,*one_name,*exp_name;
    off_t sizes[MAX_SEP_FILE];
		int filesize[MAX_SEP_FILE];
    int i, nsize, last_used=0;

    /* make an info structure for FILE* I/O and set it in the streaminf */
    inf = (struct multifd_info*) alloc( sizeof(struct multifd_info) );
    for( i=0; i< MAX_SEP_FILE; i++ ){
	inf->infos[i]=0; inf->name[i]=0; inf->used[i]=0;
    }
    *ioinfo=inf;

    /* first file is the active one by default */
    inf->active = 0;
		inf->last_used=-1;

    /* make a fd_io structure for each filename specified */
    all_names = strcpy((char*)malloc(strlen(info->dataname)+1), info->dataname);
    inf->nfile = 0;
    one_name = strtok(all_names,";");
    do{
	inf->name[inf->nfile] = strcpy((char*)malloc((int)strlen(one_name)+1), one_name);
	one_inf = (struct fd_info*) alloc(sizeof(struct fd_info));
        one_inf->curpos=0;
        one_inf->fd=-1;
	(inf->infos)[inf->nfile] = one_inf;
        inf->nfile++;
    } while( (one_name = strtok(0,";") ) != 0 );

    /* get default size of sections in MB and convert to bytes */
		filesize[0]=SEP_FILE_SIZE;
    nsize = getch("filesize","d",filesize);
    for ( i = 0 ; i <nsize ; i++ ){
        sizes[i] = filesize[i] * (double) 1024. *  (double)1024. ;
    }
    for ( i = nsize ; i <MAX_SEP_FILE ; i++ ){
        sizes[i] = (double) filesize[MAX(0,nsize-1)] *(double) 1024. * (double)1024. ;
    }

    
    /* open each dataset in turn */
    for( i=0; i< inf->nfile; i++ ){

	one_inf = inf->infos[i];

        switch( info->entrytype ){
	
	case( STREAMIN ):
	  if( !strcmp(inf->name[i],"follow_hdr") ||  
	     !strcmp(inf->name[i],"stdin") ){
	      seperr( "%s is invalid for file descriptor I/O for tag \"%s\"\n",
		     inf->name[i],info->tagname);
	  }else{
	      /* open up a file */
	      exp_name = expand_info(inf->name[i],info);
	      one_inf->fd= open(exp_name,O_RDONLY);
	      if(one_inf->fd == -1 ){
		  perror("multifd_open");
		  seperr("unable to obtain input file \"%s\" for tag \"%s\"\n",
			 exp_name,info->tagname);
	      }
	      free(exp_name);
	  }
          /* get its size */
	  one_inf->size = (sep_off_t)fsize(one_inf->fd);
	break;
	
	case( STREAMOUT ):
	  if( !strcmp( inf->name[i],"follow_hdr") || 
	     !strcmp( inf->name[i],"stdout") ){
	      seperr( "%s is invalid for file descriptor I/O for tag \"%s\"\n",
		     inf->name[i],info->tagname);
	  } else {
	      exp_name = expand_info(inf->name[i],0);
#ifdef S_IRWXU
	      one_inf->fd= open(exp_name,O_WRONLY|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO);
#else
	      one_inf->fd= open(exp_name,O_WRONLY|O_CREAT|O_TRUNC,0777);
#endif
	      if(one_inf->fd  == -1 ){
		  perror("fd_open");
		  seperr( " unable to open output datafile %s for tag %s \n",
			 inf->name[i],info->tagname);
	      }
	      free(exp_name);
	      one_inf->size = sizes[i];
	  }
	  break;
	  
	case( STREAMINOUT ):
	  exp_name = expand_info(inf->name[i],info);
	  if( info->headerbuf != 0 ){
	      one_inf->fd = open(exp_name,O_RDWR );
          }
	  if(one_inf->fd != -1 ){
	      /* use existing size */
	      one_inf->size = (sep_file_size_t)fsize( one_inf->fd );
	      if( one_inf->size != 0 ) last_used = i;
          }else{
	      /* create the file (with truncation) */
#ifdef S_IRWXU
	      one_inf->fd=open(exp_name,O_RDWR|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO); 
#else
	      one_inf->fd=open(exp_name,O_RDWR|O_CREAT|O_TRUNC,0777); 
#endif
	      /* use existing size */
	      one_inf->size = sizes[i];
	  }
	
	  if( one_inf->fd  == -1 ){
	      perror("fd_open");
	      seperr("unable to open input/output datafile %s for tag %s \n",
		     inf->name[i],info->tagname);
	  }
	  free(exp_name);
	  break;
	case( STREAMSCR ):
	  exp_name = expand_info(inf->name[i],info);
	  if( info->headerbuf != 0 ){
	      one_inf->fd = open(exp_name,O_RDWR );
          }
	  if(one_inf->fd != -1 ){
	      /* use existing size */
	      one_inf->size = (long)fsize( one_inf->fd );
	      if( one_inf->size != 0 ) last_used = i;
          }else{
	      /* create the file (with truncation) */
#ifdef S_IRWXU
	      one_inf->fd=open(exp_name,O_RDWR|O_CREAT|O_TRUNC,S_IRWXU|S_IRWXG|S_IRWXO); 
#else
	      one_inf->fd=open(exp_name,O_RDWR|O_CREAT|O_TRUNC,0777); 
#endif
	      /* use existing size */
	      one_inf->size = sizes[i];
	  }
	
	  if( one_inf->fd  == -1 ){
	      perror("fd_open");
	      seperr("unable to open scratch datafile %s for tag %s \n",
		     inf->name[i],info->tagname);
	  }
	  free(exp_name);
	  break;

		case(STREAMSOCKOUT):
			seperr("Stream socket can not be broken into multiple files \n");
    break;

       }

/*       if( !fdordinary( one_inf->fd ) ) {*/
/*	  seperr("Multi file datasets must be ordinary files\n untrue for file number: %d, called \"%s\" for tag \"%s\"\n",i,inf->name[i],info->tagname);*/
/*       }*/
    }

    /* fix up file sizes for STREAMINOUT.
     * the last file written and the unused files can be expanded
     * all the others must remain their current size */
    if( info->entrytype == STREAMINOUT ){
        for( i=last_used; i<inf->nfile ; i++ ){
	    inf->infos[i]->size = MAX( inf->infos[i]->size, sizes[i]);
        }
    }
   
    /* bogus entries for streamfd and streamfile */ 
    info->streamfd = (inf->infos)[0]->fd;
    switch( info->entrytype ){
	case( STREAMIN ):
    		info->streamfile = fdopen(inf->infos[0]->fd,"r");
	break;
	case( STREAMOUT ):
    		info->streamfile = fdopen(inf->infos[0]->fd,"w");
        break;
	case( STREAMINOUT ):
    		info->streamfile = fdopen(inf->infos[0]->fd,"r+");
        break;
	case( STREAMSOCKOUT ):
        seperr("STREAMSOCKOUT can not be of multi-file type \n");
        break;
	case( STREAMSCR ):
    		info->streamfile = fdopen(inf->infos[0]->fd,"r+");
        break;
    }

    /* it must be backward seekable since they are all regular files ! */
    info->backseekable = 1;
    
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void multifd_close( streaminf *info, void* ioinfo )
_XFUNCPROTOEND
#else
     void multifd_close( info , ioinfo )
     streaminf *info;
     void* ioinfo;
#endif
{
    struct multifd_info *inf;
    struct fd_info *one_inf;
    int i;
    
    if(ioinfo==NULL) return;

    inf=(struct multifd_info *)ioinfo;
   
    for( i=0; i<inf->nfile; i++) {
	one_inf = inf->infos[i];
        if( one_inf->fd != -1  ) close( one_inf->fd );
        free( inf->name[i] );
        free( one_inf );
    }
    
    free(inf);
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static ssize_t onefd_read( int infd, void *buffer, size_t nbytes )
_XFUNCPROTOEND
#else
static ssize_t onefd_read( infd, buffer,  nbytes )
int infd; 
void *buffer; 
size_t nbytes ;
#endif
{
    size_t total;
    ssize_t nread;
    char *cptr;
    size_t BUF = SEP_BUFSIZ;

    cptr=(char*)buffer;
    
  
    
    total=0;
    do{
	nread=read(infd,cptr, MIN(nbytes-total,BUF));
	if( nread == 0 ) break;
	
	if (nread<0) {
	    return -1;
        }

	total += nread;
	cptr += nread;

    }while( total <nbytes );
    
    return (long)total;
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
ssize_t multifd_read( streaminf *info, void* ioinfo, void *buffer, size_t nbytes )
_XFUNCPROTOEND
#else
ssize_t multifd_read( info, ioinfo,  buffer, nbytes )
     streaminf *info;
     void* ioinfo;
     void* buffer;
     size_t nbytes;
#endif
{
    struct multifd_info *inf;
    struct fd_info *one_inf;
    char* cptr;
    ssize_t total,nread;
    size_t numtry;
    /* int infd; This variable never used. */


    if( nbytes == 0 ) return 0;
    cptr = (char*)buffer;
    total = 0;

    inf=(struct multifd_info *)ioinfo;

    if( inf->active >= inf->nfile ){ 
	/* no more files left to read */
	return 0; 
    }


    do {

        one_inf = inf->infos[inf->active];

        /* either read all we need or up to the end of this file */
        numtry = MIN( nbytes-total,  one_inf->size - one_inf->curpos );

	nread = onefd_read( one_inf->fd, cptr, numtry );

        if( nread < 0 ){
	   seperr ("multifd_read() : I/O error, tag \"%s\" file number %d\n",
			info->tagname,inf->active);
        }
        one_inf->curpos += nread;

        total += nread;
        cptr  += nread;

        if( nread == 0){ 
            /* unexpected EOF (since we only read upto file size)*/
	    if( total != nbytes && inf->active == inf->nfile-1 ){
 		/* no more files and not at the end */
/*	        perror("multifd_read()");*/
                return(0); /* to be consistent with sreed behavior*/
/*	        seperr ("multifd_read() EOF for tag \"%s\" no more files \n",*/
/*			info->tagname);*/
            }
	} 

    
        if( one_inf->curpos == one_inf->size ){
            /* end of one file  go to the next */
	    inf->active++; 

            if( inf->active < inf->nfile ){
		/* next file */
                one_inf = inf->infos[inf->active];

           /* if not at the start seek to the start of the new file */
	        if( one_inf->curpos != 0 ){
		    if( lseek(one_inf->fd,0,SEEK_SET) < 0 ){
	                perror("seeking in multi_fd_read() ");
			seperr ("multifd_read() unable to seek to file start for tag \"%s\", file %d \n",info->tagname,inf->active);
	            }else{
	    	        one_inf->curpos = 0;
	            }
 	        }

            }else{ 
         /*special case when we have read to the end of the data but not over.
         We haven't committed an error and we need the option to go backwards
         so the active pointer should not go beyond the real files */
         if(total==nbytes)  inf->active--;
              
		break; /* no more files to read */
            }
       }

      /* keep going for more data */ 
    } while( total < nbytes );

    return total;

}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static ssize_t onefd_write( int outfd, void *buffer, size_t nbytes )
_XFUNCPROTOEND
#else
static ssize_t onefd_write(  outfd, buffer, nbytes )
int outfd; 
void *buffer; 
size_t nbytes ;
#endif
{

    ssize_t total;
    ssize_t nwrite;
    char *cptr;    
    size_t BUF = SEP_BUFSIZ;
    
    cptr=(char*)buffer;
 
    for (total=0; (total < nbytes) &&
	 (0 != (nwrite=write(outfd,cptr+total,MIN(nbytes-total,BUF))));
	 total += nwrite) {
	
	if (nwrite<0) {
	    return -1;
	}
    }
    
    return total;
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
ssize_t multifd_write( streaminf *info, void* ioinfo, void *buffer, size_t nbytes )
_XFUNCPROTOEND
#else
ssize_t multifd_write( info, ioinfo,  buffer, nbytes )
     streaminf *info;
     void* ioinfo;
     void* buffer;
     size_t nbytes;
#endif
{
    struct multifd_info *inf;
    struct fd_info *one_inf;
    char* cptr;
    ssize_t total,nwrite;
    size_t numtry;
    /* int infd; Never Used. */


    if( nbytes == 0 ) return 0;

    cptr = (char*)buffer;
    total = 0;

    inf=(struct multifd_info *)ioinfo;

    if( inf->active >= inf->nfile ){ 
	/* no more files left to write */
	return 0; 
    }

    one_inf = inf->infos[inf->active];
		inf->used[inf->active]=1;

    do {

        /* either write all we need or up to the end of this file */
        numtry = MIN( nbytes-total,  one_inf->size - one_inf->curpos );

	nwrite = onefd_write( one_inf->fd, cptr, numtry );

        if( nwrite < 0 ){

#ifndef EDQUOT
#ifdef CRAY
#define EDQUOT EQACT
#else
#define EDQUOT ENOSPC
#endif
#endif
#if (defined CRAY)
           if( errno == EQACT || errno == ENOSPC ){
#else
           if( errno == EDQUOT || errno == ENOSPC ){
#endif


	       /* if we ran out of space, treat it like an EOF 
                * after truncating the file to where it was before this write */
		if (-1 == ftruncate(  one_inf->fd, one_inf->curpos )) {
                    perror("multifd write()");
                    seperr("multifd write(): Unable to truncate tag \"%s\" file number %d.\n",info->tagname, inf->active);
                }
		nwrite=0;

	   }else{
	      perror("multifd_write()");
	      seperr("multifd_write() : I/O error, tag \"%s\" file number %d\n",
			info->tagname,inf->active);
	   }
        } 

        one_inf->curpos += nwrite;
        total += nwrite;
        cptr  += nwrite;

        if( one_inf->curpos == one_inf->size || nwrite == 0 ){
            /* end of space requested or EOF, go to the next file*/
	    inf->active++; 

            if( inf->active < inf->nfile ){
		/* next file */
                one_inf = inf->infos[inf->active];

                /* if not at the start seek to the start of the new file */
	        if( one_inf->curpos != 0 ){
		    if( lseek(one_inf->fd,0,SEEK_SET) < 0  ){
	                perror("seeking in multi_fd_write() "); 
			seperr ("multifd_write() unable to seek to file start for tag \"%s\", file %d \n",info->tagname,inf->active);
	            }else{
	    	        one_inf->curpos = 0;
	            }
 	        }

            }else{ 
		break; /* no more files to write */
            }
       }

      /* keep going for more data */ 
    } while( total < nbytes );
		update_out(info,ioinfo);

    return total;

}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
sep_file_size_t multifd_seek( streaminf *info, void* ioinfo, sep_off_t offset, int whence )
_XFUNCPROTOEND
#else
     sep_file_size_t multifd_seek( info, ioinfo, offset, whence )
     streaminf *info;
     void* ioinfo;
     int whence;
     sep_off_t offset;
#endif
{
    struct multifd_info *inf;
    struct fd_info *one_inf;
    sep_off_t coffset;
    sep_file_size_t position;
    off_t seekpos;
    off_t testpos;
    int seekok;
    int i;
    int ierr;

    
    inf=(struct multifd_info *)ioinfo;
    coffset = offset;
    
    ierr=0;
    
    switch( whence ){
      case SEEK_SET:
	/* figure out which file this point is in */
		
	for(inf->active=0; 
	     (inf->infos[inf->active])->size < coffset ; 
	     coffset -= (double)((inf->infos[inf->active])->size), inf->active++ ){
	    if( inf->active == inf->nfile ) {
                return ((sep_file_size_t) (-1));
	    }
        }

        one_inf = inf->infos[inf->active];
        seekpos = coffset;
	if( seekpos == one_inf->curpos ) break;

	if( lseek(one_inf->fd,seekpos,whence) < 0 ){
	    perror("multifd_seek()");
	    ierr=1;
	}else{
	    one_inf->curpos = seekpos;
	}
        break;

      case SEEK_CUR:

	/* seek relative to current point */


       if( inf->active == inf->nfile ) {
               return -1;
       }

        one_inf = inf->infos[inf->active];
	      testpos = one_inf->curpos;
        seekok = 0;
	
        do{/*switch files */
	   if( (double)coffset + (double)testpos < (double)0 ){
	      /* goto previous file */
	      inf->active--;
	      coffset = coffset + testpos;
              one_inf = inf->infos[inf->active];
	      testpos = one_inf->curpos;
	      coffset = coffset + (one_inf->size - testpos);
	      
	   }else if( (double)coffset + (double)testpos > (double)one_inf->size ){
	      /* goto next file */
	      inf->active++;
				if(inf->active==inf->nfile){
         return -1;
        }
	      coffset = coffset - ( one_inf->size - testpos);
        one_inf = inf->infos[inf->active];
	      testpos = one_inf->curpos;
	      coffset = coffset - testpos;

	   }else{
		seekok=1;
 	   }
        }while(!seekok);

        seekpos = coffset;

	if( seekpos == 0 ) break;

	if( lseek(one_inf->fd,seekpos,SEEK_CUR) <0 ){
	    perror("multifd_seek");
	    ierr=1;
	}else{
	    one_inf->curpos += seekpos;
	}
	break;

      case SEEK_END:
	assert( info->entrytype != STREAMOUT );

	/* figure out which file this point is in */
	for( inf->active=(inf->nfile-1); 
	     (inf->infos[inf->active])->size < coffset ; 
	     coffset -= (inf->infos[inf->active])->size, inf->active-- ){
	    if( inf->active <0 ){
		seperr("multifd_seek(): seeking past start of all available files for tag \"%s\" \n",info->tagname);
	    }
        }

        seekpos = coffset;
        one_inf = inf->infos[inf->active];

	if( lseek(one_inf->fd,seekpos,whence)<0 ){
	    perror("multifd_seek");
	    ierr=1;
	}else{
	    one_inf->curpos = (sep_file_size_t)fsize( one_inf->fd)-
               (sep_file_size_t)seekpos;
	}
	break;
      default:
	fprintf( stderr,"unknown seek type requested\n");
	break;
    }
    
    if( ierr == 0 ){ /* calculate current location */
	position = inf->infos[inf->active]->curpos;
	for (i = inf->active-1; i>=0 ; i--){
		 position = (double)position+ (double)(inf->infos[i]->size);
	}
	return((sep_file_size_t) (position));
    }else{
	return((sep_file_size_t) (-1));
    }
    
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
sep_file_size_t multifd_size( streaminf *info, void* ioinfo )
_XFUNCPROTOEND
#else
sep_file_size_t multifd_size( info, ioinfo )
streaminf *info;
void* ioinfo ;
#endif
{
struct multifd_info *inf;
struct fd_info *one_inf;
int i;
sep_file_size_t size;

inf=(struct multifd_info *)ioinfo;
assert( inf != 0 );

size=0;
   
for( i=0; i<inf->nfile; i++) {
    one_inf = inf->infos[i];
    size += fsize( one_inf->fd ) ;
}

return size;
}


/* Set up the function pointers in an info structure to point
 * at the fd I/O functions
 */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void init_multifd_io( streaminf* info ) 
_XFUNCPROTOEND
#else 
void init_multifd_io( info )
     streaminf * info;
#endif 
{
    info->open_func = multifd_open;
    info->close_func = multifd_close;
    info->read_func = multifd_read;
    info->write_func = multifd_write;
    info->seek_func = multifd_seek;
    info->size_func = multifd_size;
}
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void update_out( streaminf* info, void *ioinfo ) 
_XFUNCPROTOEND
#else
void update_out(info ,ioinfo) 
streaminf *info;
void *ioinfo;
#endif
{
struct multifd_info *inf;
int last,i;
char *new_data;
int new_size;

inf=(struct multifd_info *)ioinfo;last=inf->last_used;
for(i=inf->last_used+1; i < inf->nfile;i++){
	if(inf->used[i]==1) last=i;
}

if(last>inf->last_used){
	inf->last_used=last;
	new_size=0;
	for(i=0; i<= last; i++)new_size+=1+strlen(inf->name[i]);
	free(info->dataname);
	info->dataname=(char*) malloc(sizeof(char*)*new_size);
	strcpy(info->dataname,inf->name[0]);
	for(i=1;i<=last;i++){
		strcat(info->dataname,";");
		strcat(info->dataname,inf->name[i]);
	}
	info->dataname[new_size-1]='\0';
/*
	if(info->headfile != (FILE*)0 ) 
    auxputch("in","s",info->dataname,info->tagname);
*/
}

}


