/*<
   sep_insert_val_by_index 

USAGE
   int sep_insert_val_by_index(tag_history,key_index,n_values,values,all_values)

INPUT PARAMETER
   char*  - tag_history   Tag of History File
   int*   - key_index     Key Index of Header Key
   int*   - n_values      Number of Values to be written in contiguos records.  
   void*  - values        Header Values at for Key Index
OUTPUT PARAMETER
   void*  - all_values    vector of Header Values with mixed Header Types


RETURN VALUE
   -1=if fails for other reasons
    0=if successful
   +1=if tag_history is a Sep77 History File
   +2=if tag_history is a Sep3d History File but
      no matching Key Index is in Header Format File
DESCRIPTION
   return=sep_get_header_bytes(char *tag_history, int *n_bytes)
   i_byte=(record_number-1)*n_bytes+(key_index-1)*4
   sseek into the Header Format File at position i_byte
   sreed2 from Header Format File into values

CATEGORY
Lib:Sep3d:Header access

COMPILE LEVEL
DISTR
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
int sep_insert_val_by_index(char *tag_history, int *key_index, 
                    int *n_values, void *values, void *all_values)
_XFUNCPROTOEND
#else 
int sep_insert_val_by_index(tag_history,  key_index, n_values, values, all_values)
char *tag_history;
int  *key_index;
int *n_values;
void *values;
void *all_values;
#endif
{
    char *p_values, *p_all_values ;
    int n_bytes_value, n_head_bytes, first_byte, first_byte_p1;
    int i_values, i_byte, ierr,key_index_p1;


    ierr=sep_get_key_first_byte(tag_history, key_index, &first_byte);
    if(ierr!=0) {
        return ierr;
    }

    /* Read number of bytes per header record */
    ierr=sep_get_header_bytes(tag_history, &n_head_bytes);
    if(ierr!=0) {
        return ierr;
    }

    /* find position of first byte */
    /* for the moment this returns (key_index-1)*4 but it can be more general*/

    ierr= sep_get_key_first_byte(tag_history, key_index, &first_byte);
    key_index_p1=*key_index +1;
    ierr= sep_get_key_first_byte(tag_history, &key_index_p1, &first_byte_p1);
    if(ierr!=0) {
        return ierr;
    }


    /* compute number of bytes to write */
    n_bytes_value=first_byte_p1-first_byte;

    /* init pointer to data */
    p_values=values;
    p_all_values=all_values;
    p_all_values=p_all_values+first_byte;

    
    for (i_values =0; i_values < *n_values ; i_values++) {

      for (i_byte =0; i_byte < n_bytes_value ; i_byte++) {

/* copy byte */
	*p_all_values=*p_values;

/* increment pointer to data */
        p_values = p_values + 1;
        p_all_values = p_all_values + 1;
	
      }

        p_all_values = p_all_values + n_head_bytes - n_bytes_value;
    }

    
    return 0;
}

