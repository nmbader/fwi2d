/*$

=head1 NAME

sep_get_key_type - get the key type associated with an index

=head1 SYNOPSIS

C<ierr= sep_get_key_type(tag_history,key_index,key_type)>

=head1 INPUT PARAMETER

=over 4

=item    char*- tag_history     

         tag of History File

=item    int* - key_index       

         Key Index of Header Key

=back

=head1 OUTPUT PARAMETER

=over 4

=item   char*-key_type         

         Key Type of Header Key

=back

=head1 RETURN VALUE

 -1= if fails for other reasons

 0= if successful

 +1= if tag_history is an Sep77 History File

 +2= if tag_history is a Sep3d History File but
      no matching Key Index is in Header Format File

=head1 DESCRIPTION

   auxpar from the Header Format File for the parameter hdrtype"#"
   for a given key_index

=head1 SEE ALSO

L<sep_get_key_type>

=head1 LIBRARY

B<sep3d>

=cut


>*/
/*
KEYWORDS
   header 

SEE ALSO
   sep3d

AUTHOR
   SEP , July ... 1995

*/
#include <string.h>
#include "sep3d.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_get_key_type(const char *tag_history, int *key_index, char *key_type)
_XFUNCPROTOEND
#else 
int sep_get_key_type(tag_history, key_index, key_type)
char *tag_history; 
int *key_index;
char *key_type;
#endif
{
    char *tag_header[1];
    char hdrtype[32], type[32];
    int  ierr;

    ierr = sep_get_header_format_tag(tag_history, tag_header);

    /* Check for returning errors */
    if(ierr!=0) {
        return ierr;
    }

    /* Get header type for key index */
    sprintf(hdrtype, "hdrtype%d", *key_index);
    if(!auxpar(hdrtype, "s", type, *tag_header)) {
	return -1;
    }

    strcpy(key_type, type);

    free(*tag_header);

    return 0;	
	
}
