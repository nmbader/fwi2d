#define SET_SDOC 1
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Sun Aug 16 15:04:02 PDT 1998

Purpose: 

*/	 

#include "superset_internal.h"
#include<string.h>
#include <superset.h> 

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_clear_headers

USAGE
ierr= sep3d_clear_headers(char *sep3dname)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure


RETURN VALUES
0   =  if it works correctly

DESCRIPTION
Deletes the headers in memory
(check superset.h for other return values and there meaning)

>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_clear_headers(char *sep3dname) 
_XFUNCPROTOEND
#else
int sep3d_clear_headers(sep3dname ) 
char *sep3dname;
#endif 
{ 
int i1;
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s \n",sep3dname));   /* Not a valid */

if(info->headers != SEPNULL) {
	free(info->headers);
	info->headers=SEPNULL;
}
if(info->drn != SEPNULL) {
	free(info->drn);
	info->drn=SEPNULL;
}
info->nh=0;
return(SUCCESS);



}
 

/*<
sep3d_set_nh

USAGE
ierr= sep3d_set_nh(char *sep3dname, int nh)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
nh       -  int     number of traces to allocate


RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Allocates a header array that we can write into when constructing an
output or scratch file
>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_nh(char *sep3dname, int nh)
_XFUNCPROTOEND
#else
int sep3d_set_nh(sep3dname,nh)
char *sep3dname; int nh;
#endif 
{
sep_3d *info;
int i,tot;
int ierr;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s \n",sep3dname));   /* Not a valid */

if(info->file_format==REGULAR) 
	return(sepwarn(NOT_MET,"tag %s can not be regular2\n",sep3dname));
if(info->ndims ==0) 
	return(sepwarn(NOT_MET,"tag %s: ndims not set \na,sep3dname"));
if(info->nkeys==0)
	return(sepwarn(NOT_MET,"tag %s: nkeys not set \n",sep3dname));


if(info->nh!=nh){
	ierr=sep3d_clear_headers(info->name);
	if(ierr!=SUCCESS) return(ierr);
         if(nh>0) {
	   info->headers=(int*) alloc(sizeof(int)*info->nkeys*nh);
	   info->drn=(int*) alloc(sizeof(int)*nh); info->drn[0]=-1;
         }
}

info->nh=nh;
return(SUCCESS);
}
/*<
sep3d_grab_nh

USAGE
ierr= sep3d_grab_nh(char *sep3dname, int *nh)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
nh        -  int    number of traces


RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Grabs the number of headers currently stored
>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_nh(char *sep3dname, int *nh)
_XFUNCPROTOEND
#else
int sep3d_grab_nh(sep3dname,nh)
char *sep3dname; int  *nh; 
#endif 
{
sep_3d *info;
int i,tot;
;
info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

*nh=info->nh;
return(SUCCESS);
}





/*<
sep3d_grab_header_pointer

USAGE
int*= sep3d_grab_header_pointer(char *tag)

INPUT PARAMETERS
tag       -  char*   tag to get the pointer for


RETURN VALUES
 int*    pointer to the header;

DESCRIPTION
Grabs the pointer to the currently stored headers
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int* sep3d_grab_header_pointer(char *tag)
/*, int * pt)*/
_XFUNCPROTOEND
#else
int* sep3d_grab_header_pointer(tag)
/*,pt)*/
char *tag;
/*int *pt;*/
#endif 
{
sep_3d *info;
int i,tot;
float one;

info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL){
  sepwarn(0,"tag:%s  invalid struc\n",tag);
	return(SEPNULL);
}

if(info->headers==SEPNULL ) {
  sepwarn(0,"tag:%s  headers not set\n",tag);
	return(SEPNULL);;
}
	
return(info->headers);

}


/*<
sep3d_set_header_pointer

USAGE
int*= sep3d_set_header_pointer(char *tag, int nh,int *pt)

INPUT PARAMETERS
tag       -  char*   tag to get the pointer for
nh        -  int     number of  headers
pt        -  int*    pointer to the header


RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Sets a header pointer for the dataset
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_header_pointer(char *tag, int nh, int *pt)
_XFUNCPROTOEND
#else
int sep3d_set_header_pointer(tag,nh,pt)
int *pt,nh;
char *tag;
#endif 
{
sep_3d *info;
int i,tot;
float one;

info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));

if(info->nh!=0)
  return (sepwarn(NOT_MET,"tag:%s  headers have already been allocated \n",
  tag));

info->headers=pt;
info->nh=nh;
info->drn=(int*)alloc(sizeof(int)*info->nh);
info->drn[0]=-1;
return(SUCCESS);

}


/*<
sep3d_copy_headers

USAGE
int*= sep3d_compy_headers(char *tag_in, char *tag_out)

INPUT PARAMETERS
tag_in       -  char*   tag to copy the header from
tag_out      -  char*   tag to copy the header to


RETURN VALUES
SUCCESS   =  if it works correctly
(see superset.h for other return values and there meanings)

DESCRIPTION
Copies header values from one tag to another
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_copy_headers(char *tag_in, char *tag_out)
_XFUNCPROTOEND
#else
int sep3d_copy_headers(tag_in,tag_out)
char *tag_in,*tag_out;
#endif 
{
sep_3d *info_in,*info_out;
int i,tot,ierr;
float one;

info_in = tag_info_sep3d(tag_in, INQUIRE);  /* get info on this tag */
if(info_in == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag_in));
info_out = tag_info_sep3d(tag_out, INQUIRE);  /* get info on this tag */
if(info_out == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag_out));



/*if(info_in->headers==SEPNULL || info_in->nh ==0 )*/
/*  return (sepwarn(NOT_MET,"tag:%s  headers have not been allocated \n",*/
/*    tag_in));*/

if(info_out->headers!=SEPNULL){
	ierr=sep3d_clear_headers(tag_out);
	if(ierr!=SUCCESS) return(ierr);
}
ierr=sep3d_set_nh(tag_out,info_in->nh);
if(ierr!=SUCCESS) return(ierr);
if(info_in->nh>0){
	
  memcpy((void *) info_out->headers, (const void *) info_in->headers,info_in->nh*info_in->nkeys*sizeof(int));
}
sep3d_drn_copy(tag_in,tag_out);
if(ierr!=SUCCESS) return(ierr);

return(SUCCESS);

}



/*<
sep3d_set_header_block

USAGE
int*= sep3d_set_header_block(char *tag, void *block)

INPUT PARAMETERS
tag       -  char*   tag to copy the header from

OUTPUT PARAMETERS
block      -  void*  block to copy headers too (must already be allocated)


RETURN VALUES
SUCCESS   =  if it works correctly
(see superset.h for other return values and there meanings)

DESCRIPTION
Copies the headers (int and float) to a single array
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_header_block(char *tag, void *block)
_XFUNCPROTOEND
#else
int sep3d_set_header_block(tag,block)
char *tag;
void *block;
#endif 
{
sep_3d *info;
int i,tot,ierr;
float one;

info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));



if(info->headers==SEPNULL || info->nh ==0 )
  return (sepwarn(NOT_MET,"tag:%s  headers have not been allocated \n",
    tag));
	
memcpy((void *) info->headers,(const void*) block, info->nkeys*info->nh*sizeof(int));

return(SUCCESS);
}
/*<
sep3d_grab_header_block

USAGE
int*= sep3d_grab_header_block(char *tag, void *block)

INPUT PARAMETERS
tag       -  char*   tag to copy the header from

OUTPUT PARAMETERS
block      -  void*  block to copy headers too (must already be allocated)


RETURN VALUES
SUCCESS   =  if it works correctly
(see superset.h for other return values and there meanings)

DESCRIPTION
Copies the headers (int and float) to a single array
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_header_block(char *tag, void *block)
_XFUNCPROTOEND
#else
int sep3d_grab_header_block(tag,block)
char *tag;
void *block;
#endif 
{
sep_3d *info;
int i,tot,ierr;
float one;

info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));



if(info->headers==SEPNULL || info->nh ==0 )
  return (sepwarn(NOT_MET,"tag:%s  headers have not been allocated \n",
    tag));
	
memcpy((void *) block, (const void *) info->headers,info->nkeys*info->nh*sizeof(int));

return(SUCCESS);
}




/*<
sep3d_nullify

USAGE
int*= sep3d_nullify_header_pointer(char *tag)

INPUT PARAMETERS
tag       -  char*   tag to get the pointer for


RETURN VALUES
0   =  if it works correctly

DESCRIPTION
Nullifys a header pointer
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_nullify_header_pointer(char *tag)
_XFUNCPROTOEND
#else
int sep3d_nullify_header_pointer(tag)
char *tag;
#endif 
{
sep_3d *info;
int i,tot;
float one;

info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));


if(info->headers==0)
  return (sepwarn(NOT_MET,"tag:%s  headers have not been associated \n",
    tag));


info->headers=SEPNULL;
if(info->drn != SEPNULL){
	 free(info->drn);
	info->drn=SEPNULL;
}
info->nh=0;
return(SUCCESS);

}

_XFUNCPROTOBEGIN
int sep3d_header_copy(sep_3d *input, sep_3d *output)
_XFUNCPROTOEND
{
int i;
                                                                                          
if(input->nkeys != output->nkeys)
  return(sepwarn(NOT_MET,"number of keys not the same  \n"));
                                                                                          
for(i=0; i< input->nkeys; i++){
  if(0!=strcmp(input->keyname[i],output->keyname[i]))
    return(sepwarn(NOT_MET,"key %d not the same key=%s key=%s  \n",
     i,input->keyname[i],output->keyname[i]));
  if(0!=strcmp(input->keyfmt[i],output->keyfmt[i]))
    return(sepwarn(NOT_MET,"key %d not the same keyfmt=%s keyfmt=%s  \n",
     i,input->keyfmt[i],output->keyfmt[i]));
  if(0!=strcmp(input->keytype[i],output->keytype[i]))
    return(sepwarn(NOT_MET,"key %d not the same keytype=%s keytype=%s  \n",
     i,input->keytype[i],output->keytype[i]));
}
if(input->nh >0){
  if(0!=sep3d_set_nh(output->name,input->nh))
    return(sepwarn(NOT_MET,"trouble setting nh for tag %s \n",output->name));
  memcpy((void*)output->headers,(const void*)input->headers,4*input->nkeys*input->nh);
  if(input->drn!=SEPNULL)
    memcpy((void*)output->drn,(const void*)input->drn,4*input->nh);
}
return(SUCCESS);
}
