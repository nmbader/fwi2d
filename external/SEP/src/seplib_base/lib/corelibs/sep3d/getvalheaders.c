/*$

=head1 NAME

   sep_get_val_headers  -get header values


=head1 SYNOPSIS

C<ierr=sep_get_val_headers (tag_history,record_number,n_headers,header_values)>

=head1 INPUT PARAMETER

=over 4

=item   char*-  tag_history       

        tag of History File

=item   int* -  record_number      

        Header Record Number computed using sep_get_grid_window 

=item   int* -  n_headers          

        number of headers to be retrieved

=back

=head1 OUTPUT PARAMETER

=over 4

=item   void *header_values    

        Header values for n_headers

=back


=head1 RETURN VALUE

 -1 = if fails for other reasons

 0 = if successful

 +1 = if tag_history is a Sep77 History File

 +2 = if n_headers is too large

=head1 DESCRIPTION

   The values of the headers are read from the Header Values File and stored in header_values array.

=head1 COMMENTS

   return=sep_get_header_bytes(char *tag_history, int *n_bytes)
   i_byte=(record_number-1)*n_bytes
   sseek into the Header Format File at position i_byte
   sreed n_headers from Header Format File into value


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
int sep_get_val_headers(const char *tag_history, int *record_number,                                         int *n_headers, void *header_values)
_XFUNCPROTOEND
#else 
int sep_get_val_headers(tag_history, record_number, n_headers, header_values)
const char *tag_history;
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
    ierr=sep_get_header_bytes((char *) tag_history, &n_head_bytes);
    if(ierr!=0) {
        return ierr;
    }

    n_bytes_value=n_head_bytes*(*n_headers);
    seek_byte = n_head_bytes*((*record_number)-1);




     /* seek to beginning of first Header Record */
    if( sseek_block(*tag_header,n_head_bytes,((*record_number)-1),0) == -1) {
	seperr("sep_get_val_headers: sseek error seek_byte=%d\n",seek_byte);
    }

     /* read Header Records */
    if( sreed(*tag_header,header_values,n_bytes_value) != n_bytes_value) {
	seperr("sep_get_val_headers: sreed error \n" );
    }

    free(*tag_header);
            
    return 0;
}






