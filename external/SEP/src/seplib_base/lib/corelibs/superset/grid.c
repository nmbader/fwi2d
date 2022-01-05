#define SET_SDOC 1
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Sun Aug 16 15:04:02 PDT 1998

Purpose: 
To handle the grid 

*/	 

#include "superset_internal.h"
#include <superset.h> 

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_clear_grid

USAGE
ierr= sep3d_clear_grid(char *sep3dname)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure


RETURN VALUES
SUCCESS   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Deletes the grid in memory

>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_clear_grid(char *sep3dname) 
_XFUNCPROTOEND
#else
int sep3d_clear_grid(sep3dname ) 
char *sep3dname;
#endif 
{ 
int i1;
sep_3d *info;


/*Dummy  for backward compart now that grid has been eliminated*/

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s \n",sep3dname));   /* Not a valid */

return(SUCCESS);



}
 

/*<
sep3d_set_ng

USAGE
ierr= sep3d_set_ng(char *sep3dname, int ng)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
ng       -  int     number of grid elements to allocate


RETURN VALUES
SUCCESS   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Allocates a grid array that we can write into when constructing an
output or scratch file
>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_ng(char *sep3dname, int ng)
_XFUNCPROTOEND
#else
int sep3d_set_ng(sep3dname,ng)
char *sep3dname; int ng;
#endif 
{
sep_3d *info;
int i,tot,ierr;
/*Dummy  for backward compart now that grid has been eliminated*/

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s \n",sep3dname));   /* Not a valid */

if(info->file_format==REGULAR) 
	return(sepwarn(NOT_MET,"tag %s can not be regular1\n",sep3dname));
if(info->ndims ==0) 
	return(sepwarn(NOT_MET,"tag %s: ndims not set \n",sep3dname));
if(info->nkeys==0)
	return(sepwarn(NOT_MET,"tag %s: nkeys not set \n",sep3dname));


return(SUCCESS);
}
/*<
sep3d_grab_ng

USAGE
ierr= sep3d_grab_ng(char *sep3dname, int *ng)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
ng        -  int    number of grid elements


RETURN VALUES
SUCCESS   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Grabs the number of grid elements currently stored
>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_ng(char *sep3dname, int *ng)
_XFUNCPROTOEND
#else
int sep3d_grab_ng(sep3dname,ng)
char *sep3dname; int  *ng; 
#endif 
{
sep_3d *info;
int i,tot;
;
info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

if(info->ndims<1) 
  return (sepwarn(INVALID_STRUC,"tag:%s  ndims not set\n",sep3dname));

*ng=0;
for(i=1; i< info->ndims; i++) *ng=*ng*info->nwind[i];

return(SUCCESS);
}





/*<
sep3d_grab_grid_pointer

USAGE
int*= sep3d_grab_grid_pointer(char *tag)

INPUT PARAMETERS
tag       -  char*   tag to get the pointer for


RETURN VALUES
 int*    pointer to the header;

DESCRIPTION
Grabs the pointer to the currently stored grid
>
*/

/*
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int* sep3d_grab_grid_pointer(char *tag)
_XFUNCPROTOEND
#else
int* sep3d_grab_grid_pointer(tag)
char *tag;
#endif 
{
sep_3d *info;
int i,tot;
float one;

info = tag_info_sep3d(tag, INQUIRE);  
if(info == SEPNULL){
  sepwarn(0,"tag:%s  invalid struc\n",tag);
	return(SEPNULL);
}

if(info->grid==SEPNULL ) {
  sepwarn(0,"tag:%s  grid not set\n",tag);
	return(SEPNULL);;
}
	
return(info->grid);

}
*/


/*<
sep3d_set_grid_vals

USAGE
int*= sep3d_set_grid_vals(char *tag, int *vals)

INPUT PARAMETERS
tag       -  char*   tag to get the pointer for
vals      -  int*    header values


RETURN VALUES
SUCCESS   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Store the grid values
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_grid_vals(char *tag, int *vals)
_XFUNCPROTOEND
#else
int sep3d_set_grid_vals(tag,vals)
int *vals;
char *tag;
#endif 
{
sep_3d *info;
int i,tot,nc,idim,index,ng;
float one;
long long *blockb,*blocks;
sep_coord_type a1,a2,a3;

info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));

for(i=1,ng=1; i < info->ndims; i++) ng=ng*info->nwind[i];

if(ng<0) 
  return(sepwarn(NOT_MET,"window in memory not set for tag %s \n",tag));
  
for(i=0,nc=0; i <ng;i++) if(vals[i]>0) nc++;



if(nc==0) {
 SEP3D_delete_coord(info);
 return(SUCCESS);
}

if(0!=SEP3D_alloc_coord(info,nc))
  return(sepwarn(NOT_MET,"trouble allocating coord  for tag %s\n",tag));

blockb=(long long*)malloc(sizeof(long long)*(info->ndims-1));
blocks=(long long*)malloc(sizeof(long long)*(info->ndims-1));
blockb[0]=1;
for(idim=1; idim<info->ndims-1;idim++){
  blockb[idim]=blockb[idim-1]*(long long)info->n[idim];
}
blocks[0]=1;
for(idim=1; idim<info->ndims-1;idim++)
  blocks[idim]=blocks[idim-1]*(long long)info->nwind[idim];

for(i=0,nc=0; i <ng;i++) {
  if(vals[i]>0){
    a3=0;
    for(idim=1; idim< info->ndims; idim++){
      index=SEPH2C(i,blocks[idim-1],info->nwind[idim]);
/*      info->coord[nc]+=(long long)*/
/*         SEPC2H(index*info->jwind[idim]+info->fwind[idim],blockb[idim-1]);*/
         a1=( sep_coord_type)(index*info->jwind[idim]+info->fwind[idim])*
         (sep_coord_type)blockb[idim-1];
         a2=a3;
         a3=a2+a1;
    }
/*         info->coord[nc]=(sep_coord_type)a3;*/
         info->coord[nc]=(sep_coord_type)a3;
    nc++;
  }
}
free(blocks);
free(blockb);
return(SUCCESS);

}


/*<
sep3d_set_grid_pointer

USAGE
int*= sep3d_set_grid_pointer(char *tag, int ng,int *pt)

INPUT PARAMETERS
tag       -  char*   tag to get the pointer for
ng        -  int     number of  headers
pt        -  int*    pointer to the header


RETURN VALUES
SUCCESS   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Sets a grid pointer for the dataset
>
*/

/*
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_grid_pointer(char *tag, int ng, int *pt)
_XFUNCPROTOEND
#else
int sep3d_set_grid_pointer(tag,ng,pt)
int *pt,ng;
char *tag;
#endif 
{
sep_3d *info;
int i,tot;
float one;

info = tag_info_sep3d(tag, INQUIRE);  
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));

if(info->ng!=0)
  return (sepwarn(NOT_MET,"tag:%s  grid have already been allocated \n",
  tag));

info->grid=pt;
info->ng=ng;
return(SUCCESS);

}

*/

/*<
sep3d_copy_grid

USAGE
int*= sep3d_copy_grid(char *tag_in, char *tag_out)

INPUT PARAMETERS
tag_in       -  char*   tag to copy the grid from
tag_out      -  char*   tag to copy the grid to


RETURN VALUES
SUCCESS   =  if it works correctly
(see superset.h for other return values and there meanings)

DESCRIPTION
Copies grid values from one tag to another
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_copy_grid(char *tag_in, char *tag_out)
_XFUNCPROTOEND
#else
int sep3d_copy_grid(tag_in,tag_out)
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

if(0!= SEP3D_coord_copy(info_in,info_out))
  return(sepwarn(NOT_MET,"trouble copying coordinates tagin=%s tagout=%s\n",
   tag_in,tag_out));
return(SUCCESS);

}



/*<
sep3d_grab_grid_block

USAGE
int*= sep3d_grab_grid_block(char *tag, void *block)

INPUT PARAMETERS
tag       -  char*   tag to copy the grid from

OUTPUT PARAMETERS
block      -  int*  block to copy grid too (must already be allocated)


RETURN VALUES
SUCCESS   =  if it works correctly
(see superset.h for other return values and there meanings)

DESCRIPTION
Copies the grid (int) to a single array
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_grid_block(char *tag, int *grid)
_XFUNCPROTOEND
#else
int sep3d_grab_grid_block(tag,grid)
char *tag;
int *grid;
#endif 
{
sep_3d *info;
int i,tot,ierr,ng,index,ig,idim,myc,ngbig;
float one;
long long *blocks,*blockb;

info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));

for(i=1,ng=1; i < info->ndims; i++) ng=ng*info->nwind[i];
for(i=1,ngbig=1; i < info->ndims; i++)  ngbig=ngbig*info->n[i];


/*if(ng<1 || info->ncoord==0)*/
/*  return(sepwarn(NOT_MET,"can't grab grid values for tag, ncoord or nwind not set tag=%s \n",tag));*/

blocks=(long long*)malloc(sizeof(long long)*(info->ndims-1));
blockb=(long long*)malloc(sizeof(long long)*(info->ndims-1));
blockb[0]=1;
for(idim=1; idim<info->ndims-1;idim++)
  blockb[idim]=blockb[idim-1]*info->n[idim];
blocks[0]=1;
for(idim=1; idim<info->ndims-1;idim++){
  blocks[idim]=blocks[idim-1]*info->nwind[idim];
}

for(i=0; i < ng; i++) grid[i]=-1;
     
for(i=0; i <info->ncoord;i++) {
    for(idim=1,ig=0; idim< info->ndims; idim++){
      index=SEPH2C(info->coord[i],blockb[idim-1],info->n[idim]);

      ig+=
         SEPC2H((index-info->fwind[idim])/info->jwind[idim],blocks[idim-1]);
    }
    if(ig < 0 || ig >= ng) {
     return(sepwarn(NOT_MET,
      "invalid coordinate given window tag=%s ig=%d ng=%d coord[%d]=%lld\n", 
         tag,ig,ng,i,info->coord[i]));
     }
    grid[ig]=i+1;
}

free(blocks);
free(blockb);


return(SUCCESS);
}




/*<
sep3d_nullify_grid_pointer

USAGE
int*= sep3d_nullify_grid_pointer(char *tag)

INPUT PARAMETERS
tag       -  char*   tag to get the pointer for


RETURN VALUES
0   =  if it works correctly

DESCRIPTION
Nullifys a grid pointer
>
*/
/*

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_nullify_grid_pointer(char *tag)
_XFUNCPROTOEND
#else
int sep3d_nullify_grid_pointer(tag)
char *tag;
#endif 
{
sep_3d *info;
int i,tot;
float one;

info = tag_info_sep3d(tag, INQUIRE);  
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));


if(info->grid==SEPNULL)
  return (sepwarn(NOT_MET,"tag:%s  grid have not been associated \n",
    tag));


info->grid=SEPNULL;
info->ng=0;
return(SUCCESS);

}
*/

/*
_XFUNCPROTOBEGIN
int sep3d_grid_copy(sep_3d *input, sep_3d *output)
_XFUNCPROTOEND
{
int i;



if(input->ng >0){
  if(0!=sep3d_set_ng(output->name,input->ng))
    return(sepwarn(NOT_MET,"trouble allocating grid  tag=%s\n",output->name));
   
  memcpy((void*)output->grid,(const void*) input->grid,4*input->ng);
}
  
return(SUCCESS);
}
*/
