/*$

=head1 NAME

sep_copy_header_keys - copy keys from one SEP3d file to another

=head1 SYNOPSIS

#include <sep3d.h>

C<int sep_copy_header_keys(char *tag_in, char *tag_out)>

=head1 INPUT PARAMETER

=over 4

=item    char* - tag_in  

         tag of input History File

=item    char* - tag_out 

         tag of output History File

=back

=head1 RETURN VALUE

 0 if successful

 +1 if tag_history is an Sep77 History File

 -1 if fails for other reasons

=head1 DESCRIPTION

    copy keys from one SEP3d file to another  

=head1 KEYWORDS

   header keys

=head1 SEE ALSO

   sep3d

=head1 LIBRARY

B<sep3d>

=cut



AUTHOR
   Bob , December ... 1995

>*/
#include "sep3d.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_copy_header_keys(char *tag_in, char *tag_out)
_XFUNCPROTOEND
#else 
int sep_copy_header_keys(tag_in,tag_out)
char *tag_in,*tag_out;
#endif
{
    int i1,nkeys, ierr;
    char keyname[128],keytype[128],keyformat[128];


		

    /* Get tag for header format file */
	ierr=sep_get_number_keys(tag_in,&nkeys);
    	  if(ierr!=0)  return ierr;
	ierr=sep_put_number_keys(tag_out,&nkeys);
    	  if(ierr!=0)  return ierr;

    /* Check for errors in getting number of keys */
     for (i1=1;i1<=nkeys; i1++){
  	ierr=sep_get_key_name(tag_in,&i1, keyname)!=0;
	  if(ierr!=0) return ierr;
  	ierr=sep_get_key_type(tag_in,&i1, keytype);
	  if(ierr!=0) return ierr;
  	ierr=sep_get_key_fmt(tag_in,&i1, keyformat);
	  if(ierr!=0) return ierr;
  	ierr=sep_put_key(tag_out, keyname, keytype, keyformat,&i1);
	  if(ierr!=0) return ierr;
      }
    return 0;
}






