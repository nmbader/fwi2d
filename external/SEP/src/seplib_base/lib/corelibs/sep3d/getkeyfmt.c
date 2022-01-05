/*$

=head1 NAME

sep_get_key_fmt -get the format of a key

=head1 SYNOPSIS

C<int sep_get_key_fmt(tag_history,key_index,key_fmt)>

=head1 INPUT PARAMETER

=over 4

=item   char*-tag_history     

        tag of History File

=item   int* -key_index       

        Key Index of Header Key

=back


=head1 OUTPUT PARAMETER

=over 4

=item   char* - key_fmt       

        Key Name of Header Key

=back

=head1 RETURN VALUE

 -1= if fails for other reasons

 0= if successful

 +1= if tag_history is an Sep77 History File

 +2= if tag_history is a Sep3d History File but
      no matching Key Index is in Header Format File

=head1 DESCRIPTION

   auxpar from the Header Format File for the parameter hdrfmt"#"
   for a given key_index

=head1 SEE ALSO

L<sep_get_key_index>, L<sep_put_key>, L<sep_get_key_name>


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
int sep_get_key_fmt(char *tag_history, int *key_index, char *key_fmt)
_XFUNCPROTOEND
#else 
int sep_get_key_fmt(tag_history, key_index, key_fmt)
char *tag_history; 
int *key_index;
char *key_fmt;
#endif
{
    char *tag_header[1];
    char hdrfmt[32], fmt[32];
    int  ierr;

    ierr = sep_get_header_format_tag(tag_history, tag_header);

    /* Check for returning errors */
    if(ierr!=0) {
/*        free(*tag_header);*/
        return ierr;
    }

    sprintf(hdrfmt, "hdrfmt%d", *key_index);
    if(!auxpar(hdrfmt, "s", fmt, *tag_header)) {
        free(*tag_header);
	return -1;
	}

    strcpy(key_fmt, fmt);

    free(*tag_header);

    return 0;	

}
