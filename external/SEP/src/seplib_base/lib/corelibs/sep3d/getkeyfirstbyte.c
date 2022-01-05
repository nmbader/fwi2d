/*<
   sep_get_key_first_byte
USAGE
   int sep_get_key_first_byte(tag_history,key_index,key_byte)

INPUT PARAMETER
   char* -  tag_history   tag of History File
   int* - key_index       Key Index of Header Key
OUTPUT PARAMETER
   char* - key_byte       Key first byte of Header Key

RETURN VALUE
   -1=if fails for other reasons
    0=if successful
   +1=if tag_history is an Sep77 History File
   +2=if tag_history is a Sep3d History File but
      no matching Key Index is in Header Format File

DESCRIPTION
   auxpar from the Header Format File for the parameter hdrkey"#"
   for a given key_index

CATEGORY
Lib:Sep3d:Header access

COMPILE LEVEL
DISTR
>*/
/*

KEYWORDS
   header 

SEE ALSO
   sep3d

AUTHOR
   SEP , July ... 1995

REVISED
Bob-Dec'97-Changed key_first_byte to key_bytes
*/
#include "sep3d.h"
#include "streamlist.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_get_key_first_byte(char *tag_history, int *key_index, int *key_byte)
_XFUNCPROTOEND
#else 
int sep_get_key_first_byte(tag_history, key_index, key_byte)
char *tag_history; 
int *key_index;
int *key_byte;
#endif
{
    streaminf *info;
    int ierr,n_keys,i1;

    ierr=sep_get_number_keys(tag_history, &n_keys);
    if(ierr!=0) {
	return ierr;
    }
    if(*key_index > n_keys+1) seperr("\n sep_get_key_first_byte key_index > n_keys");

    info = tag_info(tag_history, TAG_INQUIRE);  /* get info on this tag */

		*key_byte=0;
		for(i1=0; i1 < *key_index-1; i1++) *key_byte += info->key_bytes[i1];


    return 0;	
	
}
