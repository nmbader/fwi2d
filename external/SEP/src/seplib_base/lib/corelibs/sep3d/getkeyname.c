/*$

=head1 NAME

sep_get_key_name  -get the key name associated with a key index

=head1 SYNOPSIS

C<ierr= sep_get_key_name(tag_history,key_index,key_name)>

=head1 INPUT PARAMETER

=over 4

=item   char*- tag_history     

        tag of History File

=item   int* - key_index       

        Key Index of Header Key

=back

=head1 OUTPUT PARAMETER

=over 4

=item   char*-key_name         

        Key Name of Header Key

=back


=head1 RETURN VALUE

 -1= if fails for other reasons

 0= if successful

 +1= if tag_history is an Sep77 History File

 +2= if tag_history is a Sep3d History File but
      no matching Key Index is in Header Format File

=head1 DESCRIPTION

   auxpar from the Header Format File for the parameter hdrkey"#"
   for a given key_index

=head1 SEE ALSO

L<sep_get_key_index>


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
int sep_get_key_name(const char *tag_history, int *key_index, char *key_name)
_XFUNCPROTOEND
#else 
int sep_get_key_name(tag_history, key_index, key_name)
char *tag_history; 
int *key_index;
char *key_name;
#endif
{

    char *tag_header[1];
    char hdrkey[255], name[255];
    int ierr;

    ierr = sep_get_header_format_tag(tag_history, tag_header);

    if(ierr!=0) {
        return ierr;
    }

    /* Find hdrkey for hdrindex passed */
    sprintf(hdrkey, "hdrkey%d", *key_index);
    if(!auxpar(hdrkey, "s", name, *tag_header)) {
	return -1;
	}

    strcpy(key_name, name);

    free(*tag_header);

    return 0;	
	
}
