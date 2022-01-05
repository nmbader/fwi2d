/*$

=head1  NAME

   sep_get_val_by_index  - get key value by index

=head1 SYNOPSIS

C<ierr=sep_get_val_by_index(tag_history,record_number,key_index,n_values,values)>

=head1 INPUT PARAMETER

=over 4

=item    char* - tag_history    

         tag of History File

=item    int* - record_number   

         Header Record Number computed using sep_get_subspace_pointers

=item    int* - key_index       

         Key Index of Header Key

=item    int* - n_values        

         Number of Values to be written in contiguous records.  

=back

=head1 OUTPUT PARAMETER

=over 4

=item    void*- values          

         Header Values of type Key Type

=back

=head1 RETURN VALUE

 -1 = if fails for other reasons

 0 = if successful

 +1 = if tag_history is a Sep77 History File

 +2 = if tag_history is a Sep3d History File but
          no matching Key Index is in Header Format File

=head1 DESCRIPTION

	Get header values by given key index

=head1 COMMENTS

   return=sep_get_header_bytes(char *tag_history, int *n_bytes)
   i_byte=(record_number-1)*n_bytes+(key_index-1)*n_bytes_value
   sseek into the Header Format File at position i_byte
   sreed2 from Header Format File into values


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
#include <stdlib.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_get_val_by_index(char *tag_history, int *record_number,                                         int *key_index, int *n_values, void *values)
_XFUNCPROTOEND
#else 
int sep_get_val_by_index(tag_history, record_number, key_index, n_values, values)
char *tag_history;
int  *record_number;
int  *key_index;
int *n_values;
void *values;
#endif
{
    char *tag_header[1];
    char key_fmt[32];
    char *p_value;
    int n_bytes_value, n_head_bytes, seek_byte, first_byte;
    int first_byte_p1,key_index_p1;
    int i_values, ierr;



    /* Get tag_header and check for returning errors */

    ierr=sep_get_header_format_tag(tag_history, tag_header);

    if(ierr!=0) {
/*        free(*tag_header);*/
        return ierr;
    }

    /* get format for key */
    ierr=sep_get_key_fmt(tag_history, key_index, key_fmt);
    if(ierr!=0) {
        free(*tag_header);
        return ierr;
    }

    /* Read number of bytes per header record */
    ierr=sep_get_header_bytes(tag_history, &n_head_bytes);

    if(ierr!=0) {
        free(*tag_header);
        return ierr;
    }

    /* find position of first byte */

    ierr= sep_get_key_first_byte(tag_history, key_index, &first_byte);
    key_index_p1=*key_index +1;
    ierr= sep_get_key_first_byte(tag_history, &key_index_p1, &first_byte_p1);

    if(ierr!=0) {
        free(*tag_header);
        return ierr;
    }

    /* compute number of bytes to write */
    n_bytes_value=first_byte_p1-first_byte;

    /* init pointer to data */
    p_value=values;

    
    for (i_values =0; i_values < *n_values ; i_values++) {

        seek_byte =(*record_number +i_values -1)*n_head_bytes + first_byte ;

        if(sseek(*tag_header,seek_byte,0) == -1) {
	    seperr("sseek error seek_byte=%d\n",seek_byte);
	}
        if(sreed2(*tag_header,p_value,n_bytes_value,key_fmt) != n_bytes_value){
	    seperr("sseek error seek_byte=%d\n",seek_byte);
	}
        
	/* increment pointer to data */
        p_value = p_value + n_bytes_value;

    }

    free(*tag_header);
    
    return 0;
}






