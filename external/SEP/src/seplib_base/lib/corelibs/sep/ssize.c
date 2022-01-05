/*$

=head1 NAME

ssize - obtain the size of a seplib file

=head1 SYNOPSIS

C<int ssize(tag)>

=head1 INPUT PARAMETERS

=over 4

=item tag - char*  

      tag for dataset

=back

=head1 DESCRIPTION

The tag argument is either the string "in" or any tag appropriate for
use with auxin().  This means either an explicit filename or a command
line redirect parameter tag=filename.  

=head1 COMMENTS

The return value is the size in bytes of the seplib data file or -1 on error.

=head1 SEE ALSO

seplib, L<sreed>, L<srite>, L<ssize>

=head1 RETURN VALUES

 -1 = is returned if there is an error.

 x =  number of bytes in file

=head1 BUGS

=head1 KEYWORDS 

size

=head1 LIBRARY

B<sep>

=cut



>*/

/*
  Dave Nichols 12/22/95  First version separate from fsize 
	Modified: Bob Clapp 7/3/96 Added ssize_block to get around 2GB limit
	Modified: Bob Clapp 7/3/96 Changed interface to functions to use
                             sep_off_t and sep_file_size_t
*/

#include <sepConfig.h>
#include <assert.h>
#include "streamlist.h"
#include <sep_main_external.h>
#include "sep_main_internal.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int ssize (char * tag)
_XFUNCPROTOEND 
#else
int ssize (tag)
char* tag ;
#endif
{
streaminf* info;

info = tag_info( tag, TAG_IN );
assert( info != 0 );

if( !info->valid ){ seperr("ssize(): invalid input tag %s\n",tag);}

/* check if this is the first I/O for this tag */
if( info->ioinf == 0 ){
    /* if it is then we open the dataset */
    (*info->open_func)( info, &(info->ioinf) );
     if( !info->valid ){seperr("ssek(): invalid input tag %s\n",tag);}
}

if( info->size_func != 0 ){
   return (int)(*info->size_func)( info, info->ioinf );
}else{
   return (int)fsize( info->streamfd );
}

}

/*$
=head1 NAME

ssize_block - obtain the number of blocks in a seplib file

=head1 SYNOPSIS

C<int ssize_block(tag,block)>

=head1 INPUT PARAMETERS

=over 4

=item tag - char*  

      tag for dataset

=item block - int

      number of bytes in a block

=back

=head1 DESCRIPTION

The tag argument is either the string "in" or any tag appropriate for
use with auxin().  This means either an explicit filename or a command
line redirect parameter tag=filename.  

=head1 COMMENTS

The return value is the size in bytes of the seplib data file or -1 on error.

=head1 SEE ALSO

seplib, L<sreed>, L<srite>, L<ssize>

=head1 RETURN VALUES

 -1 = is returned if there is an error.

 x  = number of blocks

=head1 BUGS

=head1 KEYWORDS 

size

=head1 LIBRARY

B<sep>

=cut

*/



_XFUNCPROTOBEGIN 
double ssize_info (char * tag )
_XFUNCPROTOEND 
{
streaminf* info;

info = tag_info( tag, TAG_IN );
assert( info != 0 );

if( !info->valid ){ seperr("ssize(): invalid input tag %s\n",tag);}

/* check if this is the first I/O for this tag */
if( info->ioinf == 0 ){
    /* if it is then we open the dataset */
    (*info->open_func)( info, &(info->ioinf) );
     if( !info->valid ){seperr("ssek(): invalid input tag %s\n",tag);}
}

  return((double)((*info->size_func)(info,info->ioinf )));
}




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int ssize_block (char * tag , int block_size)
_XFUNCPROTOEND 
#else
int ssize_block (tag , block_size)
char* tag ;
int block_size;
#endif
{
streaminf* info;

info = tag_info( tag, TAG_IN );
assert( info != 0 );

if( !info->valid ){ seperr("ssize(): invalid input tag %s\n",tag);}

/* check if this is the first I/O for this tag */
if( info->ioinf == 0 ){
    /* if it is then we open the dataset */
    (*info->open_func)( info, &(info->ioinf) );
     if( !info->valid ){seperr("ssek(): invalid input tag %s\n",tag);}
}

if( info->size_func != 0 ){
   return (int)((*info->size_func)( info, info->ioinf )/  
		(sep_file_size_t)block_size);
}else{
   return (int)(fsize( info->streamfd ) /(sep_file_size_t) block_size);
}
}


