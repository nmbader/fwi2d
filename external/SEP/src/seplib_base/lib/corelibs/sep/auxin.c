/*$

=head1 NAME

auxin, auxout, auxinout, auxscr,auxsockout,copy_history - returns I/O pointer 
to auxiliary data file

=head1 SYNOPSIS
	From C:
	#include <seplib.h>

	FILE *auxin(tag)
	FILE *auxtmp(tag);
	FILE *auxout(tag)
	FILE *auxinout(tag)
	FILE *auxsockout(tag)
	FILE *auxscr(tag)
	void copy_history( intag, outtag )
	char *tag, *intag, *outtag;


	From Fortran:

	integer auxin(tag)
	integer auxout(tag)
	integer auxinout(tag)
	integer auxsockout(tag)
	integer auxscr(tag)
	integer auxtmp(tag)
	integer copy_history(intag,outtag)

=head1 DESCRIPTION

These routines open an auxillary seplib dataset for input, output, or
both. For historical reasons they return a stream pointer 
(or a file descriptor for Fortran users).  These return values should not
be used for I/O. Always use the tag name to read and write data. The
return values should only be used to check for errors.

For auxin() the location of this data is determined by a line "in=filename" 
in the auxiliary input history.  The history is in turn located according the
following priorities:
Look for 'tag=history' on the command line.
Look for 'tag=history' on standard input.
Look for the file `tag' in the current directory.

For auxout(), an output history is created and
initialized as necessary. The default name `tag' of this
output history may be overrrided by specifying `tag=history'
on the command line. If the output
history already exists it will be overwritten.
The location of the data file is constructed
automatically following rules similar to those used by
output(). If the file exists it will be truncated to zero length
before starting output.

For auxinout(), the history file will be searched for using the same
rules as auxin and appended to if it exists.
If the history already exists, the data file it points at will be
reused. If you wish to append to the end of that file you should
seek to the end before writing.  If the history file doesn't exist the 
history file and data file will be created as for auxout.
The call to auxinout() must be the first use of the tag in your program.
Any other call will implicity open the dataset as either an input or
output dataset.

copy_history() is used to copy the input history from the stream
defined by tagin to the stream defined by tagout. 


=head1 EXAMPLE
	A binary file "elevations" has been generated for a 
	seismic section. The history, say Helev, describing it 
	would contain

		in="elevations"

		ne=120 units=feet datum=250

	A program written to use these elevations would be invoked

		<Hin Prog elev=Helev >Hout

	and the source for Prog might contain code such as

		auxpar("ne","d",&ne,"elev");

		auxpar("esize","d",&esize,"elev");

		auxpar("units","s",units,"elev");

		auxpar("datum","f",&datum,"elev");

		sreed("elev",elevations,ne*esize)

		for(i=0; i<ne; ++i) elevations[i] -= datum; 


=head1 DIAGNOSTICS

	Calls to auxin()/auxout()/auxinout()/auxsockout()/auxscr() from C return 
	NULL if there is no auxiliary history or data whilst
	the fortran routine returns -1 in the same situation.
	Other errors will cause program termination with a
	suitable message.

=head1 SEE ALSO

	L<auxclose>, L<auxpar>, L<auxputch>, input, L<sreed>, L<srite>, 
	L<sseek>

=head1 BUGS

	Don't forget to declare auxin an integer in Fortran 
	and Ratfor programs.

=head1 KEYWORDS   

auxillary dataset input output

=head1 LIBRARY

B<sep>

=cut

*/

/* 
 Author Dave Nichols 9/94
	New routine using new seplib IO architecture
 */

#include <sepConfig.h>
#include <stdio.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "streamlist.h"
#include "sep_main_internal.h"
#include <sep_main_external.h>
#include <sep_pars_external.h>


/*< 
auxin

Usage
	FILE *auxin(tag)-C ; int auxin(tag)

Input Parameters
tag - char* name of tag

Description
Opens auxilary input

Category
Lib:Sep:History File

Compile Level
DISTR

>*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
FILE *auxin(const char* name) 
_XFUNCPROTOEND 
#else
FILE *auxin(name) 
char *name;
#endif
{
streaminf* info;

info = tag_info( name, TAG_IN );

if( info->entrytype == STREAMOUT ){
 seperr("auxin(\"%s\"): Already opened for output only\n",name);
}
if( info->entrytype == STREAMINOUT ){
 seperr("auxin(\"%s\"): Already opened for input/output\n",name);
}
if( info->entrytype == STREAMSOCKOUT ){
 seperr("auxin(\"%s\"): Already opened for socket file\n",name);
}
if( info->entrytype == STREAMSCR ){
 seperr("auxin(\"%s\"): Already opened for scratch file\n",name);
}

if( info->valid && info->ioinf == 0 ){
    (*info->open_func)( info, &(info->ioinf) );
    if( !info->valid ){
       return 0;
   }
}

   return(info->streamfile);
}
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int fauxin(const char* name) 
_XFUNCPROTOEND 
#else
int fauxin(name) 
char *name;
#endif
{
if(auxin(name)==NULL) return -1;
else return 1;

}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
FILE *auxout( const char *name) 
_XFUNCPROTOEND 
#else
FILE *auxout(name) 
char *name;
#endif
{
streaminf* info;
info = tag_info(name,TAG_OUT);
if( info->entrytype == STREAMIN ){
 seperr("auxout(\"%s\"): Already opened for input only\n",name);
}
if( info->entrytype == STREAMINOUT ){
 seperr("auxout(\"%s\"): Already opened for input/output\n",name);
}
if( info->entrytype == STREAMSCR ){
 seperr("auxout(\"%s\"): Already opened for scratch file\n",name);
}


if( info->entrytype == STREAMSOCKOUT ){
 seperr("auxout(\"%s\"): Already opened for socket out file\n",name);
}


if( info->valid && info->ioinf == 0 ){
    (*info->open_func)( info, &(info->ioinf) );
    if( !info->valid ) return 0;
}

return(info->streamfile);
}
/*< 
auxinout

Usage
	FILE *auxinout(tag)-C ; int auxinout(tag)

Input Parameters
tag - char* name of tag

Description
Opens file for auxilary input/output

Category
Lib:Sep:History File

Compile Level
DISTR

>*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
FILE *auxinout( const char *name) 
_XFUNCPROTOEND 
#else
FILE *auxinout(name) 
char *name;
#endif
{
streaminf* info;
info = tag_info(name,TAG_INOUT);
if( info->entrytype == STREAMOUT ){
 seperr("auxinout(\"%s\"): Already opened for output only\n",name);
}
if( info->entrytype == STREAMIN ){
 seperr("auxinout(\"%s\"): Already opened for input only\n",name);
}
if( info->entrytype == STREAMSCR ){
 seperr("auxinout(\"%s\"): Already opened for scratch file\n",name);
}
if( info->entrytype == STREAMSOCKOUT ){
 seperr("auxinout(\"%s\"): Already opened for temporary file file\n",name);
}

if( info->valid && info->ioinf == 0 ){
    (*info->open_func)( info, &(info->ioinf) );
    if( !info->valid ) return 0;
}

return(info->streamfile);
}
/*< 
auxsockout

Usage
	FILE *auxsockout(tag)-C ; int auxsockout(tag)

Input Parameters
tag - char* name of tag

Description
Opens auxilary socket output file.

COMMENTS
Enables multiple sockets to be communicated between two programs by delaying
syncing of input and output until AFTER the history file, containing
the socket number has been passed.

Category
Lib:Sep:History File

Compile Level
DISTR

>*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
FILE *auxsockout( const char *name) 
_XFUNCPROTOEND 
#else
FILE *auxsockout(name) 
char *name;
#endif
{
streaminf* info;
info = tag_info(name,TAG_SOCKET);
if( info->entrytype == STREAMOUT ){
 seperr("auxsocket(\"%s\"): Already opened for output only\n",name);
}
if( info->entrytype == STREAMIN ){
 seperr("auxiscr(\"%s\"): Already opened for input only\n",name);
}
if( info->entrytype == STREAMINOUT ){
 seperr("auxsocket(\"%s\"): Already opened for input/output\n",name);
}

if( info->entrytype == STREAMSCR ){
 seperr("auxsocket(\"%s\"): Already opened for scratch file \n",name);
}

return(info->streamfile);
}


/*< 
auxscr

Usage
	FILE *auxscr(tag)-C ; int auxscr(tag)

Input Parameters
tag - char* name of tag

Description
Opens auxilary scratch file.


Category
Lib:Sep:History File

Compile Level
DISTR

>*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
FILE *auxscr( const char *name) 
_XFUNCPROTOEND 
#else
FILE *auxscr(name) 
char *name;
#endif
{
streaminf* info;
info = tag_info(name,TAG_SCR);
if( info->entrytype == STREAMOUT ){
 seperr("auxscr(\"%s\"): Already opened for output only\n",name);
}
if( info->entrytype == STREAMIN ){
 seperr("auxiscr(\"%s\"): Already opened for input only\n",name);
}
if( info->entrytype == STREAMINOUT ){
 seperr("auxscr(\"%s\"): Already opened for input/output\n",name);
}
if( info->entrytype == STREAMSOCKOUT ){
 seperr("auxscr(\"%s\"): Already opened as auxilary output socket \n",name);
}


if( info->valid && info->ioinf == 0 ){
    (*info->open_func)( info, &(info->ioinf) );
    if( !info->valid ) return 0;
}

return(info->streamfile);
}





#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
FILE *auxtmp( const char *name) 
_XFUNCPROTOEND 
#else
FILE *auxtmp(name) 
char *name;
#endif
{
streaminf* info;
 char *all_names,*one_name,*exp_name;


info = tag_info(name,TAG_SCR);
if( info->entrytype == STREAMOUT ){
 seperr("auxscr(\"%s\"): Already opened for output only\n",name);
}
if( info->entrytype == STREAMIN ){
 seperr("auxiscr(\"%s\"): Already opened for input only\n",name);
}
if( info->entrytype == STREAMINOUT ){
 seperr("auxscr(\"%s\"): Already opened for input/output\n",name);
}
if( info->entrytype == STREAMSOCKOUT ){
 seperr("auxscr(\"%s\"): Already opened as auxilary output socket \n",name);
}


if( info->valid && info->ioinf == 0 ){
    (*info->open_func)( info, &(info->ioinf) );
    if( !info->valid ){
 			return 0;
		}
}

aux_unlink(name);

return(info->streamfile);

}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int aux_unlink( const char *name) 
_XFUNCPROTOEND 
#else
int aux_unlink(name) 
char *name;
#endif
{
streaminf* info;
char *all_names,*one_name,*exp_name;

info = tag_info(name,TAG_INQUIRE);
if(info==0) seperr("no structure for tag(%s)\n",name);

all_names = strcpy((char*)malloc(strlen(info->dataname)+1), info->dataname);
one_name = strtok(all_names,";");
 do{
	unlink(one_name);
 } while( (one_name = strtok(0,";") ) != 0 );

unlink(info->headername);
free(all_names);
return(0);
}









/*< 
copy_history

Usage
	int copy_history(from,to)

Input Parameters
from - char* name of tag to copy from
to - char* name of tag to copy ro

Description
Copies content of history file to new tag

Category
Lib:Sep:History File

Compile Level
DISTR

>*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int copy_history( const char *from, const char *to )
_XFUNCPROTOEND 
#else
int copy_history( from, to )
char *from, *to;
#endif
{
streaminf *frominf, *toinf;

frominf = tag_info( from, TAG_IN );
toinf = tag_info( to, TAG_OUT );

sepstr_copyh( frominf, toinf );
return(0);
}
void grab_history(const char *tag,char *buf, int nmax,int *nsize)
{
  streaminf *infin;
  int ncpy=nmax;
  
  infin = tag_info( tag, TAG_IN );
 
  if(infin->hdrlen < ncpy) ncpy=infin->hdrlen;
   
  strncpy(buf,infin->headerbuf,ncpy);
   
  *nsize=ncpy;
 
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int tag_exists( const char *tag)
_XFUNCPROTOEND 
#endif
{
streaminf *info;
int i;
char temp1[256],temp2[256];
FILE *f;

i=1;
info = tag_info( tag, TAG_INQUIRE );
if(NULL==info){
  strcpy(temp1,tag);
  if(0==getch(temp1,"s",temp2)) strcpy(temp2,temp1);
  f=fopen(temp2,"rb");
  if(f==NULL) i=0;
  else{
    fclose(f);
  }
}

return(i);
}
