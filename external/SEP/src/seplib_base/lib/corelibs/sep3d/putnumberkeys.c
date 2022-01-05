/*$

=head1 NAME

   sep_put_number_keys  - put number of keys

=head1 SYNOPSIS

   int sep_put_number_keys(tag_history,n_keys)

=head1 INPUT PARAMETER

=over 4

=item   char*- tag_history 

        tag of History File

=item   int* - n_keys       

        number of Header Keys 

=back

=head1 RETURN VALUE

 -1 = if not Sep history file (or fails for other reasons)

 0 = if successful

 +2 = if tag_history is Sep77 History File

=head1 DESCRIPTION

auxputch nh1 to Header Format File


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

REVISIONS
	 Bob - Modified so that n_keys would be set in the header. Also allocate
   key_bytes and allow number_keys to change without destroying key_bytes
   array

*/
#include <sep3d.h>
#include "streamlist.h"


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_put_number_keys(const char *tag_history, int *n_keys)
_XFUNCPROTOEND
#else 
int sep_put_number_keys(tag_history, n_keys)
const char *tag_history;
int *n_keys;
#endif
{
    char *tag_header[1];
    int  ierr,*temp_i,i1;
		streaminf *info;


    /* Get tag for header format file */
    ierr = sep_get_header_format_tag(tag_history, tag_header);

    /* Check for errors on header format file tag */
    if(ierr!=0) return ierr;

	  info = tag_info(tag_history, TAG_INQUIRE);  /* get info on this tag */

		if(info->n_key!=-1) { /* this isn't the first time put_number_keys has been
    called */
			temp_i=(int *) malloc (sizeof(int) * info->n_key);
	
			for(i1=0;i1< info->n_key; i1++) temp_i[i1]=info->key_bytes[i1];
			free(info->key_bytes);
			info->key_bytes=(int *) malloc (sizeof(int) * (*n_keys));
			for(i1=0;i1< MIN(info->n_key,*n_keys);i1++){ info->key_bytes[i1]=temp_i[i1];}
			free(temp_i);	
		}
		else info->key_bytes=(int *) malloc (sizeof(int) * (*n_keys));
		info->n_key=*n_keys;	


    /* Putch number of keys into header format file */
    auxputch("n1", "d", n_keys, *tag_header);

    free(*tag_header);

    return 0;
}






