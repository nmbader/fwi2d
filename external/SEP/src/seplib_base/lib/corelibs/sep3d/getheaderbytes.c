/*<
sep_get_header_bytes


usage
   int sep_get_header_bytes(tag_history,n_bytes)

INPUT PARAMETER
   char* - tag_history   tag of History File
OUTPUT PARAMETER
   int* - n_bytes         number of bytes in Header Record

RETURN VALUE
    0=if successful
   +1=if tag_history is an Sep77 History File
   -1=if fails for other reasons

DESCRIPTION
	Get the number of bytes in a header

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

REVISED:
Bob-Dec'97: Changed key_first_byte to key_bytes

*/
#include "sep3d.h"
#include "streamlist.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_get_header_bytes(char *tag_history, int *n_bytes)
_XFUNCPROTOEND
#else 
int sep_get_header_bytes(tag_history, n_bytes)
char *tag_history;
int *n_bytes;
#endif
{
    streaminf *info;
    int n_keys, ierr,i1;

    ierr=sep_get_number_keys(tag_history, &n_keys);
    if(ierr!=0) {
	return ierr;
    }

    info = tag_info(tag_history, TAG_INQUIRE);  /* get info on this tag */
		
		*n_bytes=0;
		for(i1=0;i1<n_keys;i1++){
			 *n_bytes+=info->key_bytes[i1];
			
		}
		


    
    return 0;
}






