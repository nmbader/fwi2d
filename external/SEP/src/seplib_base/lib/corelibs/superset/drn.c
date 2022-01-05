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
sep3d_set_drn

USAGE
ierr= sep3d_set_drn(char *sep3dname,  int *drn)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure

OUTPUT PARAMETERS
drn       -  int*    data record number for the dataset

RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Set the data record number for the headers currently in memory
>
*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_drn(char *sep3dname,int *drn)
_XFUNCPROTOEND
#else
int sep3d_set_drn(sep3dname, drn)
char *sep3dname;
int *drn;
#endif 
{
int i1,ierr;
sep_3d *info;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));


if(info->drn!=SEPNULL) free(info->drn);
if(info->nh==0) return(sepwarn(NOT_MET,
	"size of headers (sep3d_set_nh) must be called before sep3d_set_drn  (%s)\n",
	sep3dname));

info->drn=(int*)alloc(sizeof(int)*info->nh);
for(i1=0; i1 < info->nh; i1++)info->drn[i1]=drn[i1];
return(SUCCESS);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_drn

USAGE
ierr= sep3d_grab_drn(char *sep3dname,  int *drn)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure

OUTPUT PARAMETERS
drn       -  int*    data record number for the dataset

RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Find the data record number of the dataset
>
*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_inorder(const char *sep3dname,int *inorder)
_XFUNCPROTOEND
#else
int sep3d_grab_inorder(sep3dname, inorder)
const char *sep3dname;
int *inorder;
#endif 
{
int i1,ierr;
sep_3d *info;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

*inorder=info->in_order;
return(SUCCESS);
}




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_drn(char *sep3dname,int *drn)
_XFUNCPROTOEND
#else
int sep3d_grab_drn(sep3dname, drn)
char *sep3dname;
int *drn;
#endif 
{
int i1,ierr;
sep_3d *info;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));


if(info->drn==SEPNULL){
	sep3d_print(info);
	return(sepwarn(NOT_MET,"drn has not been set in %s\n",sep3dname));
}
for(i1=0; i1 < info->nh; i1++) drn[i1]=info->drn[i1];
return(SUCCESS);
}




/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_drn_set

USAGE
ierr= sep3d_drn_set(char *sep3dname)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure


RETURN VALUES
0   =  if drn is set to reasonable values

-1  =  drn not set to reasonable values


DESCRIPTION
Sets that the dataset is in order
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_drn_set(char *sep3dname)
_XFUNCPROTOEND
#else
int sep3d_drn_set(sep3dname)
char *sep3dname;
#endif 
{
int i1,ierr;
sep_3d *info;


info = tag_info_sep3d(sep3dname, INQUIRE); 
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));
if(info->in_order==1) return(1);

if(info->drn == SEPNULL) return(-1);
for(i1=0; i1 < info->nh; i1++){
   if(info->drn[i1] < 1) return(-1);
}
return(SUCCESS);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_set_inorder

USAGE
ierr= sep3d_set_inorder(char *sep3dname)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure


RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Sets that the dataset is in order
>
*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sep3d_set_inorder_flag(const char *sep3dname, int flag)
_XFUNCPROTOEND
#else
int sep3d_set_inorder(sep3dname)
const char *sep3dname;
#endif 
{
int i1,ierr;
sep_3d *info;


info = tag_info_sep3d(sep3dname, INQUIRE); 
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

info->in_order=flag;
if(flag != 0 && info->drn != SEPNULL){
	 free(info->drn);
	 info->drn=SEPNULL;
}

return(SUCCESS);
}
int sep3d_set_inorder(const char *sepname) {
    sep3d_set_inorder_flag(sepname, 1);
   return 0;
}
int sep3d_unset_inorder(const char *sepname) {
    sep3d_set_inorder_flag(sepname, 0);
   return 0;
}

_XFUNCPROTOBEGIN
int sep3d_drn_copy(char *sepin,char *sepout)
_XFUNCPROTOEND
{
int i1,ierr,i;
int *ndim,n123,idim[99];
sep_3d *in,*out;


in = tag_info_sep3d(sepin , INQUIRE); 
out= tag_info_sep3d(sepout, INQUIRE); 
if(in == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sepin));
if(out == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sepout));

if(in->nh!=out->nh) 
  return (sepwarn(INVALID_STRUC,"tag:%s and %s attempt to copy drn when headers not same length\n",sepin,sepout));

if(in->nh==0) return(0);

if(in->drn==SEPNULL) { /*we have to construct the drn list */
  if(in->file_format==GRID)
    return(sepwarn(NOT_MET,
      "attempted to copy drn when input in order and gridded dataset (must create drn first"));
  else{
    if(in->fwind[1]==-1) 
      return(sepwarn(NOT_MET,
      "attempted to copy drn when input in order and window has not ben set (must set window first"));
    for(i1=1,n123=1; i1 < in->ndims; i1++)n123=n123*in->nwind[i1];
    for (i1=0; i1 < n123; i1++){
      h2c(i1,&in->nwind[1], in->ndims-1, idim);
      for(i=1; i < in->ndims; i++) idim[i-1]=idim[i-1]*in->jwind[i]+in->fwind[i];
      c2h(&out->drn[i1],&in->n[1],in->ndims-1,idim);
    }
  }
}
else if (out->in_order!=1){
 memcpy((void*)out->drn,(const void*)in->drn,sizeof(int)*in->nh);
 }


  
  
   
in->in_order=1;

return(SUCCESS);
}


