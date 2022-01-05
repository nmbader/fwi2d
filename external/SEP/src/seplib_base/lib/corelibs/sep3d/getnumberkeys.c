/*$

=head1 NAME

   sep_get_number_keys  -get number of keys

=head1 SYNOPSIS

C<ierr= sep_get_number_keys(tag_history,n_keys)>

=head1 INPUT PARAMETER

=over 4

=item   char* - tag_history   

        tag of History File

=back

=head1 OUTPUT PARAMETER

=over 4

=item   int* - n_keys         

         number of Header Keys described in Header Format File

=back

=head1 RETURN VALUE

 -1 = if fails for other reasons

 0 = if successful

 +1 = if tag_history is an Sep77 History File


=head1 DESCRIPTION

    auxpar from Header Format File for parameter nh1


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

Modified:
Bob-Dec'97:PROTOTYPES, Changed first_bytes to key_bytes.  Error checking
so we don't auxpar from STREAMOUT.  Goal: Allow piping and coding style
that easy to transfer accross platforms.


*/
#include <string.h>
#include "sep3d.h"
#include "streamlist.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_get_number_keys(const char *tag_history, int *n_keys)
_XFUNCPROTOEND
#else 
int sep_get_number_keys(tag_history, n_keys)
const char *tag_history;
int *n_keys;
#endif
{
    streaminf *info;
    char *tag_header[1];
    char hdrtype[255], type[255];
    int nkeys, ierr, i_key,bytes=0;

    /* Get tag for header format file */

    ierr = sep_get_header_format_tag(tag_history, tag_header);

    /* Check for errors on header format file tag */
    if(ierr!=0) {
			*n_keys = 0;
			return ierr;
    }

    info = tag_info(tag_history, TAG_INQUIRE);  /* get info on this tag */



    if(info->n_key == -1){

		if( ((info->entrytype == STREAMOUT)))
			seperr("get_number_keys:Can not get number of keys from STREAMOUT \n");


    /* Get number of keys from header format file */
      if(!auxpar("n1", "d", &nkeys, *tag_header)){
	seperr("\n sep_get_number_keys: Must specify n1 in %s\n", *tag_header);  
      }

      info->n_key=nkeys;
      info->key_bytes= (int *) malloc((nkeys)*(int)sizeof(int));

    for(i_key=0; i_key < nkeys; i_key++) { 
			sprintf(hdrtype, "hdrtype%d", i_key+1);
			if(!auxpar(hdrtype, "s", type, *tag_header)) 
	    	seperr("\n sep_get_number_keys: Must specify %s in %s\n", hdrtype, *tag_header);  
   
			if(strcmp(type,"scalar_float") || strcmp(type,"scalar_int")) 
				info->key_bytes[i_key]=4;

	 		else {
	    	/* Not supported data type */
        free(*tag_header);
	    	return -1;
	    	/* Future releases will include other data types */
			}
    }


    }
    else {
      nkeys=info->n_key;
    }

    
    *n_keys = nkeys;

    free(*tag_header);

    return 0;
}

