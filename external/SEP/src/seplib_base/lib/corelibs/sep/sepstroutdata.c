/*
 *
 * LOGIC:
 *         module for preparing sepstream output data file.
 *         External entries defined here are:
 *
 *     open_outstream( info )        initialize output data stream.
 *     outname( info )  generate output file name (also used by sepstrinoutdata)
 *
 * entries in info filled in by this routine are:
 *          dataname, format_num and all the IO function pointers
 *
 *
 * Author:  D. Nichols  8-94   SEP
 * Revised: S.A. Levin  5-6-95 Mobil/SEP  Don't use non-posix strdup()
 * Recised: R.G. Clapp  7-19-97 Bob Prototype=good
 *
 */
#include <sepConfig.h>
#include <stdio.h>
#include <fcntl.h>
#include <limits.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <assert.h>
#include<time.h>

#include <sepcube.h>
#include "sep_main_internal.h"
#include "streamlist.h"

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




/* INTERNAL PROTOTYPES */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static char* expand_name( char*);
static char* makenames( char*,char*,int);
static void stdoutname( streaminf*);
void outname(streaminf*);
_XFUNCPROTOEND 
#else
static char* expand_name();
static char* makenames();
static void stdoutname();
void outname();
#endif


char parambuf4[4096];
char datapth[4096];
extern char** sepxargv;

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void open_outstream( streaminf *info )
_XFUNCPROTOEND
#else
void open_outstream( info )
     streaminf* info;
#endif
{
    
    assert( info->entrytype == STREAMOUT );
    assert( info->headfile  != (FILE*)0 );
    
    if( !strcmp( info->tagname, "out" )){
	stdoutname(info);
    }else{
	outname( info );
    }

    /* set up the pointers to the I/O routines */
    init_io( info );
}

/* expand a name to a full path, if required */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static char* expand_name( char* file )
_XFUNCPROTOEND 
#else
static char* expand_name( file )
char* file;
#endif
{
 static char expbuf[4096];
 char *rc;

   strcpy( expbuf, file );

   if(expbuf[0] != '.' || expbuf[1] != '/'){
     fullnm(expbuf,sizeof(expbuf));
   }

   rc = (char *) malloc((int)strlen(expbuf)+1);
   strcpy(rc, expbuf);
   return ( rc );
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
char* sep_tail(char* name)
_XFUNCPROTOEND 
#else
char * sep_tail(name)
     char *name;
#endif
{
    char *result;
    
    result=strrchr(name,'/');
    
    if(result == (char *) NULL) return(name);
    else return(result+1);
}
#if !defined (HAVE_STDLIB_H)
extern char *putenv ();
#endif /* HAVE_STDLIB_H  */
/* datapath may be a semi-colon separated list of paths */
/* we must construct a semi-colon seprated list of filenames from it */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static char* makenames( char* path, char* tail, int randit )
_XFUNCPROTOEND 
#else
static char* makenames( path, tail,randit )
char* path;
char* tail;
int randit
#endif
{
    char one_name[4096];
    char *all_paths;
    char *one_path;
    char *exp_name;
    char *all_names;
    char my_name[4982];
    int file_count=0,fs=0,ds=0,id=0;
		int filesize[MAX_SEP_FILE],dirsize[MAX_SEP_FILE],this_dir;
		int ncount,i,nfiles;

    char path_use[19999],tt[9999];
    int ierr=getch("env_path","s",tt);
    if(ierr==1) {
       sprintf(path_use,"%s{%s}","$",tt);
       sprintf(tt,"%s=%s",tt,path);
       putenv(tt);
       getch_add_string(tt);
    }
    else  strcpy(path_use,path);
    all_paths = strcpy((char*) malloc((int)strlen(path_use)+1), path_use );
    all_names =(char*)malloc(((int)strlen(all_paths)+(int)strlen(tail)+4)*MAX_SEP_FILE+1); /* the +4 is file number + ; + @ */


	/* 
		The general idea for multi-file writes is to grab the list of directories 
    to write to and the size of files to write. We default to SEP_FILE_SIZE per drive
    and file size of SEP_FILE_SIZE (which can be overwritten by "dirsize" and 
    "filesize".  We  then create a new datapath containing an expanded
    list of names. The multi-file option can be initialized by either a
    ';' in the datapath list or by a filesize < dirsize
  */

    if(0==getch("nfiles","d",&nfiles)) nfiles=1;
    if(nfiles>MAX_SEP_FILE)  seperr("nfiles to large, max %d \n",MAX_SEP_FILE);
		dirsize[0]=SEP_DIR_SIZE; filesize[0]=SEP_FILE_SIZE;
		ncount=getch("dirsize","d",&dirsize);
		for(i=ncount;i<nfiles;i++) dirsize[i]=dirsize[MAX(ncount-1,0)];
		ncount=getch("filesize","d",&filesize);
		for(i=ncount;i<nfiles;i++) filesize[i]=filesize[MAX(ncount-1,0)];
	  
	
    one_path = strtok( all_paths,";");
		id=-1; 
    do{
			this_dir=0; /* initialize how much of this directory we have covered*/
			id++;
			ds=dirsize[id];
			do{
     		if( file_count == 0 ){
	   			/* no number for the first file */
	   			sprintf(one_name,"%s%s@", one_path, tail );
        }else{
	   			/* the subsequent ones get a number for uniqueness */
	   			sprintf(one_name,"%s%s@%d", one_path, tail, file_count );
        }
					exp_name = expand_name( one_name );
                if(randit==1){
                     mkrandom_string(exp_name,my_name);
                 }
                else strcpy(my_name,exp_name);
        	if( file_count == 0 ){
	   			strcpy( all_names,my_name);
        }else{
	   			strcat(strcat( all_names,";"),my_name);
        }
        free(exp_name);
				this_dir+=filesize[file_count];
				file_count++;
			}while(this_dir< ds && file_count < nfiles);
/*    }while(( ( one_path = strtok(0,";")) != 0 ) && file_count < nfiles);*/
    }while(( ( one_path = strtok(0,";")) != 0 ) );

   free(all_paths);
    return all_names;
} 

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static void stdoutname( streaminf* info )
_XFUNCPROTOEND 
#else
static void stdoutname( info )
     streaminf *info;
#endif
{
    
    char *front,hold[300],prog_name[256],ttemp[256];
   /* char *tail; Never Used */
		time_t tempt;	
		int put_end,i,j;
		/* int num; Never Used */
    
    
    /* now identify the output data stream. If we are going down a pipe
     * then it is "follow_hdr".
     * otherwise it must be a file.
     */
    
    if( getch("out","s",parambuf4)){
        /* if explicitly specified then leave it alone, unless it is "stdout" */
	if( !strcmp(parambuf4,"stdout") ){
	  if( info->headfile == stdout ){
	    /* set to stdout and header going down stdout */
	    info->dataname = (char *) malloc((int)strlen("follow_hdr")+1);
	    strcpy(info->dataname,"follow_hdr");
	  }else{
	    /* explicitly going down stdout but header has gone elsewhere */
	    info->dataname = (char *) malloc((int)strlen("stdout")+1);
	    strcpy(info->dataname,"stdout");
          }
	}else{
	    /* its a file name, get a full path */
	    info->dataname = expand_name(parambuf4);
	}
    }else{
	
	/* Send data to stdout if header going down pipe */
	if ( info->isapipe || isapipe(fileno(info->headfile)) ){
	    info->dataname = (char *) malloc((int)strlen("follow_hdr")+1);
	    strcpy(info->dataname,"follow_hdr");
	    
	}else{   
  	    /* retrieve or construct a default path name */
	    /* make out= from datapath+command name */
	    
	    if ( !strcmp(datapath(datapth), "stdout")){
	      info->dataname = (char *) malloc((int)strlen("follow_hdr")+1);
	      strcpy(info->dataname,"follow_hdr");
	    }else{

                datapath(datapth);
		if( (getch("head","s",parambuf4) && isordinary(parambuf4))
		   || ( 0 < findnm(fileno(info->headfile), 
				   parambuf4, sizeof(parambuf4) )  )){
		    /* modify slightly */
		    front = strrchr(parambuf4,'/');
		    if(front == 0 ) front = parambuf4; else front++;
		    
		    if((*front) == 'H') front++;

                    info->dataname = makenames( datapth, front,0 );

		} else { 
		    /* couldn't find (use) header name */
				  
		    /* check if command was pathnamed */
				i=0;j=0;
			 	while(sepxargv[0][i] != '\0'){
					if(sepxargv[0][i] == '/')  j=0;
					else{
						prog_name[j]=sepxargv[0][i]; j++;
					}
					i++;
				}
				prog_name[j]='\0';
/*        mkrandom_string(prog_name,ttemp);*/
         datapath(datapth); 
        info->dataname = makenames( datapth, prog_name,1 );

		}

	    }
	}
    }
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int get_data_name(char *tag,char *string_out)
_XFUNCPROTOEND 
#endif
{
streaminf *info;
char *name;
char junk[PATH_MAX];
int i;

info = tag_info( tag, TAG_INQUIRE );
if(info==0) info = tag_info( tag, TAG_OUT);
if( strchr( info->headername, '|' ) != 0 ) return(1);
else if ( strchr( info->headername, ':' ) != 0 )  return(1);
else if( strcmp( info->headername, "stdout" ) == 0 )  {
  i=findnm(fileno(info->headfile),junk,sizeof(junk));
  if(i <=0) return(1);
}
else{
  name=expand_name(info->headername);
  strncpy(junk,sep_tail(name),sizeof(junk)-1); junk[PATH_MAX-1] = '\0';
  strcpy(string_out,junk);i=0;
  free(name);
 return(0);
}

return(i);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void mkrandom_string(char *string_in,char *string_out)
_XFUNCPROTOEND 
#endif
{
		time_t tempt;	
		int put_end,i,j;
    char temp3_file[1080];
    char temp2file[1080];
    char temp_file[1080];
    char *junk;

    /* couldn't find (use) header name */

        time(&tempt);
        srand(tempt);
/*        put_end=rand();*/
/*        sprintf(string_out,"%s%d",string_in,put_end);*/
/*        sprintf(string_out,"%s",string_in);*/
        sprintf(string_out,"%sXXXXXX",string_in);
        if(-1 == mkstemp(string_out)) {
           perror("mkrandom_string(): mkstemp()");
        }
        unlink(string_out);

            
/*         mkstemp(string_out);*/
/*        sprintf(string_out,"%s@",string_out);*/
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
void outname(streaminf* info)
_XFUNCPROTOEND 
#else
void outname( info )
     streaminf* info;
#endif
{
char temp_ch[4096],temp2_ch[4096];

   sprintf(temp_ch,"%s.out",info->tagname);
   if(1==getch(temp_ch,"s",temp2_ch)){
	     info->dataname = (char *) alloc((int)strlen(temp2_ch)+1);
	     strcpy(info->dataname,temp2_ch);
     return;
   }
   

    if( info->isapipe ){
	info->dataname = (char *) alloc((int)strlen("follow_hdr")+1);
	strcpy(info->dataname,"follow_hdr");
	return;
    }

    /* get the datapath specifier */
    datapath(datapth);
    info->dataname = makenames( datapth, sep_tail(info->headername),0);
}

