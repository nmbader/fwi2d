/*$
=head1 NAME

doc - program self documentation


=head1 SYNOPSIS

doc(SOURCE)

=head1 INPUT PARAMETERS 

=over 4

=item SOURCE - char* 

      Location of source file

=back

=head1 DESCRIPTION

SEPlibs documentation mechanism. When a program is called
with no arguments and no in or out the begining of the file
is searched for documentation.  For C, the first comment

=head1 COMMENTS

Standard practice is to give program usage description in comment
lines at the beginning of the source file.  Doc reads such comments
from the file name it is passed and prints them (on the stderr file)
when the program is invoked without inputs or arguments.  Thus the
prospective user need only type the program name to obtain a parameter
description.

If the filename is a relative path the documentation will be searched for
in the list of directories defined by the compile time 
variable DEFAULT_DOC_PATH followed by the list of directories in 
the environment variable "SEP_DOC_PATH". Both of these should be 
colon separated lists of directories.

Copying is done up to the end of the first comment block or
preprocessor control symbol (# for C, % for Ratfor).  Output is piped
through more so that long descriptions can be viewed.

In the interest of increased portability,  provision is made to input the
comments directly as a (long) character string. (This is fairly messy,
however, and compiler length limitations may interpose.)
  
The standard include file <sep.startup> generates the call
`doc(SOURCE)'. The standard gnu-makefile rules for building seplib
will automatically define SOURCE.

 If the programmer predefines SOURCE via a #define
they should only do so if it is not previously defined., e.g.

 #ifndef SOURCE 

 #define SOURCE "./segy/programs/Segy.c"

 #endif

=head1 SEE ALSO

 seplib, more

=head1 DIAGNOSTICS

If SOURCE has not been defined,  the cc compiler will yell.

=head1 BUGS

Unpredictable results if comments don't appear at the beginning of the
program.  Fortran version needs % to recognize end of comment stream.

=head1 KEYWORDS 

	self doc documentation manuals

=head1 LIBRARY

B<sep>

=cut
*/

/*
 *      C program self-documentation
 *
 *	prints first comment statements to a # symbol or end of comment
 *	useful for self-documenting programs
 *	ascii name of the source file must be specified
 *
 * Modified 1/26/83  S. Levin : Added search for initial c-com  in input string
 *                              signifying inline literal documentation
 * Modified 1/31/83  S. Levin : Incorporated fix and rewrite communicated by
 *				Jeff Thorson. Put both C and Ratfor code in
 *				same file.
 * Modified 2/7/83   S. Levin : Incorporated suggestion by Dave Hale to pipe
 *				documentation through  more  command.
 * Modified 7/14/83  S. Levin : added err and perror diagnostic calls.
 *				Changed more flag from -f to -c .
 * Modified 8/8/83   S. Levin : Split fortran routine into separate file
 *				to permit independent linkage.
 * Modified 9/16/85  S. Levin : Convex changes to handle inline doc. Apparently
 *				Convex fprintf does fail on long >128 strings!
 * Modified 10/9/91  D. Nichols : make doc.c handle c and ratfor comments
 *				so that doc_.c is just a stub.
 * Modified 7/15/97  R. Clapp : Changed dup2 to typing to STD_C 
 *                              While ago added the ability to read
 *                              Fortran90 comments
 *                              
 * Modiefied 4/1/99   R.Clapp   Added support for POD (Perl) in F90
 * Modiefied 6/1/99   R.Clapp   CHangd to GNU standard for ifdef
 *  
 */

#include <sepConfig.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* FROM GNU Bash */
#if !defined (MAXPATHLEN) && defined (HAVE_LIMITS_H)
#include <limits.h>
#endif /* !MAXPATHLEN && HAVE_LIMITS_H */

#if !defined (MAXPATHLEN) && defined (HAVE_SYS_PARAM)
#include <sys/param.h>
#endif /* !MAXPATHLEN && HAVE_SYS_PARAM */

#if !defined (MAXPATHLEN) && defined (PATH_MAX)
#define MAXPATHLEN PATH_MAX
#endif /* !MAXPATHLEN && PATH_MAX */

/* Yecch!  Who cares about this gross concept in the first place? */
#if !defined (MAXPATHLEN)
#  define MAXPATHLEN 1024
#endif /* MAXPATHLEN */





#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <sep_main_external.h>
#include "../include/sep_main_internal.h"
#include "../include/sepstream.h"
static sep_self_doc *prog_list=0;


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
char* find_doc( char* );
char* get_doc_path( int );
sep_self_doc *get_last_doc_line (void);
void page_doc(void);
_XFUNCPROTOEND
#else
char* find_doc();
char* get_doc_path();
sep_self_doc *get_last_doc_line ();
void page_doc();
#endif

_XFUNCPROTOBEGIN
void sep_add_doc_line (const char* line)
_XFUNCPROTOEND
{
                                                                                    
sep_self_doc *newinf;
                                                                                    
  newinf=get_last_doc_line();
  newinf->line=(char*)malloc(sizeof(char)*((int)strlen(line)+1));
  strcpy(newinf->line,line);
  newinf->next=0;
}
_XFUNCPROTOBEGIN
sep_self_doc *get_last_doc_line (void)
_XFUNCPROTOEND
{
sep_self_doc *cur,*newinf,*last;
int i;
                                                                                    
                                                                                    
  cur=prog_list;
 i=0;
  while(cur!=0){
    last=cur;
    cur=cur->next;
  i++;
  }
                                                                                    
if(i==0){
  newinf = (sep_self_doc*) malloc(sizeof(sep_self_doc));
  prog_list=newinf;
}
else{
  newinf = (sep_self_doc*) malloc(sizeof(sep_self_doc));
  last->next=newinf;
}
  return(newinf);
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int doc (char* name)
_XFUNCPROTOEND
#else
int doc (name)
char *name;
#endif
{
	FILE *s,*morein;
	int c,c1 ;
	extern int sepxargc;

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
_XFUNCPROTOEND
#else
extern FILE *popen();
#endif 


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
_XFUNCPROTOEND
#else
extern int dup2();
#endif
char first[1024];



  if ( sepxargc == 1  &&  !redin() ){
  if(prog_list!=0)  page_doc();
else{
		if( strncmp(name,"/*",2) == 0 ||    /* C comment */
                    strncmp(name,"#",1) == 0 ||    /* ratfor comment */
                    strncmp(name,"!",1) == 0 ||    /* fortran90 comment */
		    strncmp(name,"c",1) == 0  ||
	  	    strncmp(name,"C",1) == 0  )      /* fortran comments */
		{
		/* like all shell scripts this takes a bit of explanation */
		/* the problem here is that the default stdout for the popen */
		/* command is the stdout of the calling program which is not */
		/* the screen. The 1>&2 construction under  sh (NOT csh!) forces */
		/* file  1 (stdout) to be a direct copy of file 2 (stderr). */
		/* In this way we make the more command output to stderr */
		/* as desired.						*/
		/*							*/
		/*	morein = popen("sh -c 'more 1>&2'","w");  */
		/* 							*/
		/* Better method from Jeff Thorson:			*/
			if(-1 == dup2(2,1)) {
				perror("doc()");
				seperr("doc(): dup2 bomb\n");
				}
			if(NULL == (morein = popen("more","w")))
			  {
			   perror("doc()");
			   seperr("doc() unable to pipe to more command\n");
			  }
			while (c = *(name++)) putc(c,morein);
			if(ferror(morein))
			   {
			    perror("doc()");
			    seperr("doc() I/O error piping to more command\n");
			   }
			pclose(morein);
			exit(1);
		}
		else 
		{
			char	*doc_name = NULL;

			doc_name = find_doc( name );

			if ( (s=fopen(doc_name,"r")) == NULL) {
				perror("doc()");
				seperr("doc() source not at  %s \n",name);
				}
                                        first[0] = '\0';
					if(((char *) NULL) ==  fgets(first,5,s)) {
                                           perror("doc read problem ");
					}
					if(0==strncmp(first,"!!$",3) ||
					 0==strncmp(first,"#$",2) ||
					 0==strncmp(first,"/*$",3) )
						strcpy(first,"ExtractPOD | /usr/bin/pod2text | more");
					else strcpy(first,"more");
					fclose (s);
					s=fopen(doc_name,"r");
					
				
			
			c1 = ' ';

		/**** morein = popen("sh -c 'more 1>&2'","w"); ****/

			dup2(2,1);
/*			if(NULL == (morein = popen("more","w")))*/
			if(NULL == (morein = popen(first,"w")))
			  {
			   perror("doc()");
			   seperr("doc() unable to pipe to more command\n");
			  }
			while ( (c=getc(s)) != EOF  &&   
			    !(c1 == 'C' && c == '%') &&  /* end of fortran */
			    !(c1 == 'c' && c == '%') &&  /* end of fortran */
			    !(c1 == '!' && c == '%') &&  /* end of fortran90 */
			    !(c1 == '#' && c == '%') &&  /* end of ratfor */
			    !(c1 == '*' && c == '/')  ) /* end of C */
				{
				if(ferror(s))
				  { 
				   perror("doc()");
   				   seperr("doc() I/O error reading source file %s\n",name);
				  }
				c1 = c;
				putc (c,morein);
				if(ferror(morein))
				  { 
				   perror("doc()");
				   seperr("doc() I/O error piping to more command\n");
				  }
				}
			pclose(morein);
			fclose (s);
			exit(1);
		}
  }

}
return 0;
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
char *find_doc( char *name )
_XFUNCPROTOEND
#else
char *find_doc( name )
char    *name;
#endif
{
        static char     file_name[ MAXPATHLEN ];
        struct stat     s;
        int             path_index;
        char            *path;
				int   myl,i,j;
				char            base[1024],name2[4096];

			
				
				strcpy(base,DEFAULT_DOC_PATH);

				

        /* name begins with / */
        if( name[0] == '/'){
						/*check to see if the base name is the same as DEFAULT_DOC_PATH*/
						myl=0;i=0;
						while(i<MIN(strlen(base),strlen(name)) && myl==0){
							if(base[i]!=name[i]) myl=1;
							i+=1;
						}
						if(myl==1) return name; /*not the path listed in DEFAULT_DOC_PATH*/
						while(name[i]=='/' && i < strlen(name)) i+=1;
						for(j=i;j< strlen(name);j++) name2[j-i]=name[j];
							name2[j-i]='\0';
				}
				else strcpy(name2,name);
        /* search path */
        path_index = 0;
        while( ( path = get_doc_path( path_index ) ) != NULL )
        {
                strcpy( file_name, path );
                strcat( file_name, "/" );
                strcat( file_name, name2);
                if( stat( file_name, &s ) != -1 )
                        return file_name;
                path_index++;
        }

        return name;
}

#ifndef DEFAULT_DOC_PATH
#define DEFAULT_DOC_PATH "/usr/local/SEP/src/"
#endif
#ifndef DEFAULT_SEP_PATH
#define DEFAULT_SEP_PATH "/usr/local/SEP/src/"
#endif
#if !defined (HAVE_STDLIB_H)
extern char *getenv ();
#endif /* HAVE_STDLIB_H  */

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
char* get_doc_path( int index )
_XFUNCPROTOEND
#else
char* get_doc_path( index )
int     index;
#endif
{
        static char path[ MAXPATHLEN ];
        char    *doc_path = NULL;
        char    *p;
        char    *q;


        /* the compiled in path is always first */
	strcpy( path, DEFAULT_DOC_PATH );
   strcat( path, ":"  ); strcat( path, DEFAULT_SEP_PATH );

        if( ( doc_path = getenv( "SEP_DOC_PATH" ) ) != NULL ){
              /* add the environment variable */
	      strcat( path, ":"  );
              strcat( path, doc_path );
	}

        q = p = path;

        while( index-- >= 0 )
        {
                q = p;
                p = strchr( p, ':' );
                if( p == NULL )
                        break;
                else
                        p++;
        }

        if( index >= 0 )
                return NULL;

        if( p != NULL ) p--;
        if( p != NULL && *p == ':' ) *p = '\0';

        return q;

}

                                                                                    
/*FROM SU's docpkge
                                                                                    
****************************************************************************
Authors: Jack Cohen, Center for Wave Phenomena, 1993, based on on earlier
                                                                                    

versions by:
SEP: Einar Kjartansson, Stew Levin CWP: Jack Cohen, Shuki Ronen
HRC: Lyle
****************************************************************************/
                                                                                    
                                                                                    
void page_doc(void)
{
        FILE *fp;
    char *pager;
    char cmd[32];
    sep_self_doc *line;
                                                                                    
    if ((pager=getenv("PAGER")) != (char *)NULL)
      sprintf(cmd,"%s 1>&2", pager);
    else
      sprintf(cmd,"more 1>&2");
                                                                                    
     line=prog_list;
                                                                                    
        fflush(stdout);
       /*  fp = popen("more -22 1>&2", "w"); */
       /*  fp = popen("more  1>&2", "w"); */
        fp = popen(cmd, "w");
       while(line!=0){
         (void)fprintf(fp, "%s\n", line->line);
         line=line->next;
       }
        pclose(fp);
        exit(EXIT_FAILURE);
}


