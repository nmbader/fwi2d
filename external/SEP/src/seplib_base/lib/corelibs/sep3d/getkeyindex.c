/*$

=head1 NAME
   sep_get_key_index  - get the index of a key

=head1 SYNOPSIS

C<int sep_get_key_index(tag_history,key_name,key_index)>

=head1 INPUT PARAMETER

=over 4

=item   char* - tag_history      

        tag of History File

=item   char* - key_name         

        Key Name of Header Key

=back

=head1 OUTPUT PARAMETER

=over 4

=item   int*  - key_index        

        Key Index of Header Key

=back

=head1 RETURN VALUE

 -1= if fails for other reasons

 0= if successful

 +1= if tag_history is an Sep77 History File

 +2= if tag_history is a Sep3d History File but
      no matching Key Name is in Header Format File

=head1 DESCRIPTION

   auxpar from the Header Format File for the parameter hdrindex"#"
   for a given key_name

=head1 SEE ALSO

L<sep_get_key_name>

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

Bob 5/99 - Replaced strncmp with strcmp

*/
#include "sep3d.h"
#include <string.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_get_key_index(const char *tag_history, char *key_name, int *key_index)
_XFUNCPROTOEND
#else 
int sep_get_key_index(tag_history, key_name, key_index)
const char *tag_history; 
char *key_name;
int *key_index;
#endif
{
    register int i_index;
    char *tag_header[1];
    char hdrkey[255], name[255];
    int nkeys, ierr, not_found = 1;

    ierr = sep_get_number_keys(tag_history, &nkeys);

    /* Check for returning errors */
    if(ierr!=0) {
        return ierr;
    }
    *key_index = 0;

    ierr=sep_get_header_format_tag(tag_history, tag_header);
    if(ierr!=0) {
        return ierr;
    }

    /* Read from header format file different hdrkeys names
       and compare with key_name passed.  */
    i_index=0;
    while(i_index++<=nkeys && not_found) {
	sprintf(hdrkey, "hdrkey%d", i_index);
	auxpar(hdrkey, "s", name, *tag_header);
	not_found = strcmp(name, key_name);
    } 
    
    if(not_found) {
	*key_index = 0;
	return 3;
    }

    *key_index = i_index-1;

    free(*tag_header);

    return 0;	
	
}
