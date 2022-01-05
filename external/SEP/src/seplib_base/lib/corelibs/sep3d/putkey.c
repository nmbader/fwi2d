/*$

=head1 NAME

   sep_put_key  - write key info to a tag


=head1 SYNOPSIS

   int sep_put_key (tag_history,key_name,key_type,key_fmt,key_index)

=head1 INPUT PARAMETER

=over 4

=item   char* - tag_history      

        tag of History File

=item   char* - key_name         

        Key Name of Header Key

=item   char* - key_type         

        Key Type of Header Key

=item   char* - key_fmt          

        Key Format of Header Key

=item   int*  - key_index        

        Key Index of Header Key

=back

=head1 RETURN VALUE

 -1= if fails for other reasons

 0= if successful

 +1= if tag_history is an Sep77 History File

 +2= if tag_history is a Sep3d History File but

=head1 DESCRIPTION

   auxputch into the Header Format File the key key_index"#"
   hdrkey"#" with hdrtype"#" and hdrfmt"#"

=head1 SEE ALSO

L<sep_get_key_index>, L<sep_get_key_fmt>, L<sep_get_key_name>


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

Modified:
Bob-(Dec'97)Removed check for key duplication (auxpar on output not allowed 
    [necessary for piping]). If someone wants it back should include keynames 
    in sepstream structure....Not a bad idea anyway. Also had it set the
    number of bytes in the key in the sepstream structure for the same reason.

*/
#include <sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include "sep3d.h"
#include "streamlist.h"
#define KEY_STR_LEN 15

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_put_key(const char *tag_history, const char *key_name, const char *key_type,                             const char *key_fmt, int *key_index)
_XFUNCPROTOEND
#else 
int sep_put_key(tag_history, key_name, key_type, key_fmt, key_index)
const char *tag_history; 
const char *key_name;
const char *key_type;
const char *key_fmt;
int *key_index;
#endif
{
    char *tag_header[1];
    char hdrkey[255], hdrtype[32], hdrfmt[32],space1[KEY_STR_LEN],space2[KEY_STR_LEN],space3[KEY_STR_LEN];
    int ierr=0,i,nlen;
		streaminf *info;
   char temp_ch[2045];
	

		info = tag_info(tag_history, TAG_INQUIRE);  /* get info on this tag */
		
		/*CHECK TO MAKE SURE NUMBER OF KEYS HAS BEEN SET IN OUTPUT HEADER 
     AND THAT WE ARE NOT TRYING TO WRITE A KEYNUMBER OUTSIDE RANGE */
		if(info->n_key ==-1) 
     seperr("sep_put_key:must call sep_put_number_keys before sep_put_key\n");
		if(*key_index < 1 || *key_index > info->n_key) 
			seperr("sep_put_key:invalid key number (outside 1 < key_index[%d] < n_key[%d])\n",*key_index,info->n_key);

		/*SET THE NUMBER OF HEADER BYTES FOR THIS KEY IN THE key_bytes array */
		if(0==strcmp ("scalar_float",key_type)){
			 info->key_bytes[*key_index-1]=4;
		}	
		else if(0==strcmp ("scalar_int",key_type)){
			 info->key_bytes[*key_index-1]=4;
		}
		else seperr("put_key:Unrecognized key_type=%s keyname=%s \n",key_type,key_name);



   /* Check for returning errors */
   ierr=sep_get_header_format_tag(tag_history, tag_header);
   if(ierr!=0) return ierr;
    
  if(*key_index < 10) nlen=KEY_STR_LEN-1;
  else nlen=KEY_STR_LEN-2;
  strcpy(space1,"             ");
  strcpy(space2,"             ");
 
  space1[MAX(2,nlen-(int)strlen(key_name))]='\0';
  space2[MAX(nlen-(int)strlen(key_type),2)]='\0';


  /*Auxputch key to header format file */
/*
  auxputhead(tag_header[0],
  "\t hdrkey%d=\"%s\"%shdrtype%d=\"%s\"%shdrfmt%d=\"%s\"\n",
   *key_index,key_name,space1,*key_index,key_type,space2,*key_index,key_fmt);
 sprintf(temp_ch,
  "\t hdrkey%d=\"%s\"    hdrtype%d=\"%s\"   hdrfmt%d=\"%s\"\n",
   *key_index,key_name,*key_index,key_type,*key_index,key_fmt);
*/
 sprintf(temp_ch,
  "\t hdrkey%d=\"%s\"    hdrtype%d=\"%s\"   hdrfmt%d",
   *key_index,key_name,*key_index,key_type,*key_index);

   auxputch(temp_ch,"s",key_fmt,tag_header[0]);

if(0!=strcmp(tag_history,"outtag") && 0== strcmp("data_record_number",key_name))
  return(0);

  free(tag_header[0]);
  return 0;	
}



