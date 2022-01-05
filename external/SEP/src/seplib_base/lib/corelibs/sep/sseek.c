/*$

=head1 NAME

sseek - seek to a position in a SEPlib dataset


=head1 SYNOPSIS

	int sseek (tag,offset,whence)

=head1 INPUT PARAMETERS

=over 4

=item char* - tag 

      name of history file;

=item int - offset 

      offset from whence given in machine types size for file i/o;

=item int   - whence 

      from where 

=back

=head1 RETURN VALUES

=over 4

 -1 =  if fails

 x = current position if successful

=back

=head1 DESCRIPTION

Moves file pointers inside a given SEPLIB dataset


=head1 COMMENTS

The tag argument is either the string "in","out", or any tag appropriate for
use with auxin() or auxout().  
This means either an explicit filename or a command
line redirect parameter tag=filename.  

sseek() sets the seek pointer associated with the open seplib dataset
or  device referred to by the tag according to the
value supplied for whence.  whence must be one of  the  following 
constants defined in <unistd.h>:

               SEEK_SET
               SEEK_CUR
               SEEK_END

If whence is SEEK_SET, the seek pointer  is  set  to  offset
bytes.   If  whence  is SEEK_CUR, the seek pointer is set to
its current location plus offset.  If  whence  is  SEEK_END,
the seek pointer is set to the size of the file plus offset.


=head1 SEE ALSO
seplib, L<sreed>, L<srite>, L<auxclose>, L<sseek_block>

=head1 DIAGNOSTICS

If an error occurs the return value will be -1. A diagnostic
error should be printed.


=head1 KEYWORDS 
seek position

=head1 LIBRARY

B<sep>

=cut


*/

/*
 * Author Dave Nichols (SEP) 9/94
	 Modified Bob Clapp 7/3/96 Added sep_block to get around 2GB limit and
                             lack of 2GB+ ints on some machines
	 Modified Bob Clapp 7/3/96 Switched to using sep_file_size_t and
                              sep_off_t for the same reason
   
 */

#include <sepConfig.h>
#include <stdio.h>
#include <string.h>
#include<math.h>
#include <sys/types.h>

#include "streamlist.h"
#include <sep_main_external.h>
#include "sep_main_internal.h"
#include <assert.h>

/*static int discard();*/
#ifdef _POSIX_SOURCE

#include <unistd.h>
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int sseek(const char *tag, int offset, int whence)
_XFUNCPROTOEND 
#else
int sseek( tag, offset, whence )
char *tag;
int offset;
int whence;
#endif

#else /* not posix source */

#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int sseek(const char *tag, int offset, int whence)
_XFUNCPROTOEND
#else
int sseek( tag, offset, whence )
char *tag;
int offset;
int whence;
#endif

#endif

{

 streaminf *info;
 /* int ierr; Never used */
 /* int positon; Never Referenced, neither is position */
 sep_off_t offset_send;
 

assert( tag != 0 );

info = tag_info(tag,TAG_INQUIRE);  /* get the info about this tag */

assert( info != 0 );

/* check if this is the first I/O for this tag */
if( info->ioinf == 0 ){
    /* if it is then we open the dataset */
    (*info->open_func)( info , &(info->ioinf) );
    if( !info->valid ){seperr("sseek(): invalid tag %s\n",tag);}
}

	offset_send = (sep_off_t)offset;
	
return (int)((*info->seek_func)(info,info->ioinf, offset_send, whence ) );


}

/*$

=head1 NAME

sseek_block - seek to a position in a SEPlib dataset by blocks


=head1 SYNOPSIS

int sseek_block (tag,offset,blocksize, whence)

=head1 INPUT PARAMETERS

=over 4

=item char* - tag 

      name of history file;

=item int - offset 

      offset from whence given in machine types size for file i/o;

=item int   - block 

      size 

=item int   - whence 

      from where 

=back

=head1 RETURN VALUES

 -1 =  if fails

 x = current position if successful

=head1 DESCRIPTION

Moves file pointers inside a given SEPLIB dataset by blocks. Has the
advantage of handling file sizes bigger than 2 GB


=head1 COMMENTS

The tag argument is either the string "in","out", or any tag appropriate for
use with auxin() or auxout().  
This means either an explicit filename or a command
line redirect parameter tag=filename.  

sseek() sets the seek pointer associated with the open seplib dataset
or  device referred to by the tag according to the
value supplied for whence.  whence must be one of  the  following 
constants defined in <unistd.h>:

               SEEK_SET
               SEEK_CUR
               SEEK_END

If whence is SEEK_SET, the seek pointer  is  set  to  offset
bytes.   If  whence  is SEEK_CUR, the seek pointer is set to
its current location plus offset.  If  whence  is  SEEK_END,
the seek pointer is set to the size of the file plus offset.


=head1 SEE ALSO
seplib, L<sreed>, L<srite>, L<auxclose>, L<sseek>

=head1 DIAGNOSTICS

If an error occurs the return value will be -1. A diagnostic
error should be printed.


=head1 KEYWORDS 
seek position

=head1 LIBRARY

B<sep>

=cut


*/




#ifdef _POSIX_SOURCE
#include <unistd.h>
#else
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int sseek_block(const char *tag, int offset, int block_size,int whence)
_XFUNCPROTOEND 
#else
int sseek_block( tag, offset, block_size,whence )
char *tag;
int offset;
int block_size;
int whence;
#endif
{

return(int)sseek_block_d(tag,offset,block_size,whence);

}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
long long sseek_l_l(const char *tag, long long offset,int whence)
_XFUNCPROTOEND 
#endif
{

 streaminf *info;
 

assert( tag != 0 );
info = tag_info(tag,TAG_INQUIRE);  /* get the info about this tag */
assert( info != 0 );


/* check if this is the first I/O for this tag */
if( info->ioinf == 0 ){
    /* if it is then we open the dataset */
    (*info->open_func)( info , &(info->ioinf) );
    if( !info->valid ){seperr("sseek(): invalid tag %s\n",tag);}
}

return(((*info->seek_func)(info,info->ioinf, (sep_off_t)offset,whence)));
}






#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
double sseek_block_d(const char *tag, int offset, int block_size,int whence)
_XFUNCPROTOEND 
#else
double sseek_block_d( tag, offset, block_size,whence )
char *tag;
int offset;
int block_size;
int whence;
#endif

{

 streaminf *info;
  int ierr;
 int position,big_block;
 double position_d;
 int whence_send,blocks,remainder;
 sep_off_t offset_send;
 sep_file_size_t big_size,position_back=0,current_pos,temp,seek_request;
 /* sep_file_size_t position_add; Never Used */
 

/*fprintf(stderr,"SEEKING %d %d \n",whence,offset*block_size);*/
assert( tag != 0 );

info = tag_info(tag,TAG_INQUIRE);  /* get the info about this tag */

assert( info != 0 );


/* check if this is the first I/O for this tag */
if( info->ioinf == 0 ){
    /* if it is then we open the dataset */
    (*info->open_func)( info , &(info->ioinf) );
    if( !info->valid ){seperr("sseek(): invalid tag %s\n",tag);}
}


	temp=(sep_off_t)block_size;
	seek_request=(sep_file_size_t)temp*(sep_file_size_t)offset;
	current_pos = ((*info->seek_func)(info,info->ioinf, (sep_off_t)0,1));
  if(current_pos < 0) {
			return((double)sepwarn(-1,"SEEK ERROR: trouble obtaining current position\n"));
  }



	if(whence==SEEK_SET){ /*change to a relative seek */
		seek_request=seek_request-current_pos;
		whence_send=SEEK_CUR;
	}
	else whence_send=whence;


	while(fabs(seek_request) > (sep_file_size_t)MAX_INT_SIZE){
		if(seek_request >0.) offset_send=(sep_off_t)MAX_INT_SIZE;
		else offset_send=-(sep_off_t)MAX_INT_SIZE;
		position_back = ((*info->seek_func)(info,info->ioinf, offset_send, whence_send ) );
		if(position_back<0){ /*seek error */
			return((double)sepwarn(-1,"SEEK ERROR: tag=%s offset=%d whence=%d \n        current_pos=%g blocksize=%d \n ",
				tag,offset_send,whence_send,current_pos,block_size));
		 }
	  	seek_request = seek_request - (sep_file_size_t) offset_send;
		whence_send=SEEK_CUR;
	}

	offset_send= (sep_off_t) seek_request;


	position_back = ((*info->seek_func)(info,info->ioinf, offset_send,whence_send ));


	if(position_back<0){ /*seek error */
			return((double)sepwarn(-1,"SEEK ERROR(2): tag=%s offset=%d whence=%d  \n                  block=%d curent position=%g \n \n",
				tag,(int)offset_send,(int)whence_send,block_size,current_pos));
	}

	position_d = (double)((double)position_back / (sep_file_size_t)block_size);


	return(position_d);

}











/*<
file_position


USAGE 
ierr=file_position(tag,block_size,blocks,remainder)


RETURN VALUES
0   if successfual
1   if block size > MAX_INT_SIZE
-1  if seek fails

INPUT PARAMETERS
tag          -  char*   tag of the file
block_size   -  int     size of blocks

OUTPUT PARAMETERS
blocks       -  int*     number of blocks
remainder    -  int*     remainder  


DESCRIPTION
Get the position in a file

>*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int file_position(const char *tag, int block_size, int *blocks,int *remainder)
_XFUNCPROTOEND 
#else
int file_position( tag, block_size, blocks,remainder )
char *tag;
int block_size;
int *blocks;
int *remainder;
#endif
{

 streaminf *info;
 int position;
 int return_it=0;
 int whence_send,big_block;
 sep_off_t offset_send;
 sep_file_size_t big_size,position_back=0;
 /* sep_file_size_t position_add; Never Used */
 
assert( tag != 0 );

info = tag_info(tag,TAG_INQUIRE);  /* get the info about this tag */

assert( info != 0 );

/* check if this is the first I/O for this tag */
if( info->ioinf == 0 ){
    /* if it is then we open the dataset */
    (*info->open_func)( info , &(info->ioinf) );
    if( !info->valid ){seperr("sseek(): invalid tag %s\n",tag);}
}

position_back = ((*info->seek_func)(info,info->ioinf, (sep_off_t)0,1));
if(position_back <0) return((int)position_back);



big_size=(position_back/(sep_file_size_t)block_size);
if(big_size>MAX_INT_SIZE){ /* we have more than MAX_INT_SIZE blocks 
                              get the remainder accuarte anyway */
	big_block=MAX_INT_SIZE/block_size;
	return_it+=1;
	while(big_size > MAX_INT_SIZE){
		big_size-=(sep_off_t)big_block;
		position_back-=(sep_off_t)block_size*(sep_off_t)big_block;
	}
}
		

position=(int)big_size;

big_size=(position_back-(sep_off_t)position*(sep_file_size_t)block_size);
if(big_size>MAX_INT_SIZE){
	 return(sepwarn(-1,
  "(file_position), remainder is larger than MAX_INT_SIZE\n"));
}
*remainder=(int)big_size;
*blocks=position;
return(return_it);
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int sseekable(const char *tag)
_XFUNCPROTOEND 
#else
int sseekable( tag)
char *tag;
#endif
{
streaminf *info;

assert( tag != 0 );

info = tag_info(tag,TAG_INQUIRE);  /* get the info about this tag */

assert( info != 0 );

if(info->isapipe) return(1);
else return(0);

}
