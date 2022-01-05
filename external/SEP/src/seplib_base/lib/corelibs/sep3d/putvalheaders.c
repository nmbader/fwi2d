/*$

=head1 NAME

   sep_put_val_headers  - write a header block

=head1 SYNOPSIS

   int sep_put_val_headers (tag_history,record_number,n_headers,header_values)

=head1 INPUT PARAMETER

=over 4

=item   char*   - tag_history    

        tag of History File

=item   int*    - record_number   

        Header Record Number computed using sep_get_subspace_pointers

=item   int*    - n_headers       

        number of headers to be retrieved

=item   float*  - header_values 

        Header values for n_headers

=back

=head1 RETURN VALUE

 -1= if fails for other reasons

 0= if successful

 +1= if tag_history is a Sep77 History File

 +2= if n_headers is too large

=head1 DESCRIPTION

   The Header Record Number are computed using the function
   sep_get_numbers.

=head1 COMMENTS

   return=sep_get_header_bytes(char *tag_history, int *n_bytes)
   i_byte=(record_number-1)*n_bytes
   sseek into the Header Format File at position i_byte
   srite2 n_headers from Header Format File into value


=head1 SEE ALSO

L<sep_get_val_headers>

=head1 LIBRARY

B<sep3d>

=cut


>*/

/*


KEYWORDS
   header keys

SEE ALSO
   sep3d

AUTHOR
   SEP , July ... 1995

*/
#include "sep3d.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_put_val_headers(char *tag_history, int *record_number,                                         int *n_headers, void *header_values)
_XFUNCPROTOEND
#else 
int sep_put_val_headers(tag_history, record_number, n_headers, header_values)
char *tag_history;
int  *record_number;
int  *n_headers;
void *header_values;
#endif
{
    char *tag_header[1];
    int n_bytes_value, n_head_bytes, seek_byte;
    int ierr;


    /* Get tag_header and check for returning errors */
    ierr=sep_get_header_format_tag(tag_history, tag_header);
    if(ierr!=0) {
        return ierr;
    }



    
    /* Read number of bytes per header record */
    ierr=sep_get_header_bytes(tag_history, &n_head_bytes);
    if(ierr!=0) {
        return ierr;
    }


    n_bytes_value=n_head_bytes*(*n_headers);
    seek_byte = n_head_bytes*((*record_number)-1);


/*
   if(-1==sseek_block(*tag_header,n_head_bytes,((*record_number)-1),0)){
     fprintf(stderr,"call %s %d %d \n",tag_history,n_head_bytes,*record_number);
     seperr("sep_put_val_headers: sseek error \n");
   }
*/
     


    /* rite Header Records */
    if( srite(*tag_header,header_values,n_bytes_value) != n_bytes_value) {
	seperr("sep_put_val_headers: srite error \n" );
    }

        
    free(*tag_header);
    return 0;
}

