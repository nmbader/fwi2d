#define SET_SDOC 1
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Thu Jan  8 08:20:53 PST 2004

Purpose: 

*/	 

#include "superset_internal.h"
#include<string.h>
#include <superset.h> 


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_set_coord_vals

USAGE
ierr= sep3d_set_coord_vals(char *sep3dname,  int *coords)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure

OUTPUT PARAMETERS
coords       -  int*    coordinate value (ndim-1 * ncoord)

RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Set the coordinate values for the current in block
>
*/


_XFUNCPROTOBEGIN
int sep3d_grab_ncoord(char *sep3dname,int *ncoord)
_XFUNCPROTOEND
{
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

*ncoord=info->ncoord;

return(SUCCESS);
}
_XFUNCPROTOBEGIN
int sep3d_set_coord_vals(char *sep3dname,int *coords)
_XFUNCPROTOEND
{
int i1,ierr,idim,i;
sep_3d *info;
long long *block;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

if(info->coord==SEPNULL) 
 return(sepwarn(NOT_MET, "coords have not been allocated (%s)\n", sep3dname));

block=(long long*) malloc(sizeof(long long)*(info->ndims-1));
for(idim=1,block[0]=1; idim < info->ndims-1; idim++) 
  block[idim]=block[idim-1]*info->n[idim];

for(i1=0,i=0; i1 < info->ncoord; i1++){
  info->coord[i1]=0;
  for(idim=0; idim < info->ndims-1; idim++,i++){
    info->coord[i1]+=(sep_coord_type)(coords[i]-1)*block[idim];
  }
}

free(block);
return(SUCCESS);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_coord_vals

USAGE
ierr= sep3d_grab_coord_vals(char *sep3dname,  int *coords)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure

OUTPUT PARAMETERS
coords       -  int*    coordinate (ndim-1*ncoord)

RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Find the data record number of the dataset
>
*/


_XFUNCPROTOBEGIN
int sep3d_grab_coord_vals(char *sep3dname,int *coords)
_XFUNCPROTOEND
{
int i1,ierr,i,idim,j,nsz;
sep_3d *info;
long long *block,*n;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));




if(info->coord==SEPNULL){
  if(info->file_format!=REGULAR || info->nwind[1] <1)
	return(sepwarn(NOT_MET,"coord has not been set in %s\n",sep3dname));
  if(0!=SEP3D_wind_coords(info,&nsz))
    return(sepwarn(NOT_MET,"trouble window to coords \n"));
}

block=(long long*) malloc(sizeof(long long)*(info->ndims-1));
n=(long long*) malloc(sizeof(long long)*(info->ndims-1));
for(idim=1,block[0]=1; idim < info->ndims-1; idim++) {
  block[idim]=block[idim-1]*info->n[idim];
  n[idim]=info->n[idim];
}

for(i1=0,i=0; i1 < info->ncoord; i1++){
  for(idim=0; idim < info->ndims-1; idim++,i++){
/*    j=(int)((info->coord[i1]/block[idim])/(long long)n[idim]);*/
/*    coords[i]=*/
/*     (int)((long long)(info->coord[i1]/block[idim])-(long long)j*(long long)n[idim]);*/
      coords[i]=((int)(info->coord[i1]/block[idim]))%info->n[idim+1]+1;
  }

}
free(block); free(n);
return(SUCCESS);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_set_coordh

USAGE
ierr= sep3d_set_coordh(char *sep3dname,  int *coords, int *coordb, int block)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure

OUTPUT PARAMETERS
coords       -  int*    coordinate value part1

coordb       -  int*    coordinate value part2

block        -  int     block (coord=coordb*block+coords)

RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Set the coordinate values for the current in block
>
*/


_XFUNCPROTOBEGIN
int sep3d_set_coordh(char *sep3dname,long long *cooordh)
_XFUNCPROTOEND
{
int i1,ierr;
sep_3d *info;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));


if(info->coord==SEPNULL) 
 return(sepwarn(NOT_MET, "coords have not been allocated (%s)\n", sep3dname));

memcpy((void*) info->coord,(const void*) cooordh,info->ncoord*sizeof(long long));

return(SUCCESS);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_coord

USAGE
ierr= sep3d_grab_coordh(char *sep3dname,  int *coords, int *coordb, int *block)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure

OUTPUT PARAMETERS
coords       -  int*    coord part 1 of dataset

coordb       -  int*    coord part 2 of dataset

block        -  int*    block (coord=coords+block*coordb)

RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Find the data record number of the dataset
>
*/


_XFUNCPROTOBEGIN
int sep3d_grab_coordh(char *sep3dname,long long *coords)
_XFUNCPROTOEND
{
int i1,ierr;
sep_3d *info;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));


if(info->coord==SEPNULL)
	return(sepwarn(NOT_MET,"coord has not been set in %s\n",sep3dname));

memcpy((void*)coords,(const void*) info->coord,info->ncoord*sizeof(long long));
return(SUCCESS);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
SEP3D_delete_coord

USAGE
ierr= SEP3D_delete_coord(sep_3d *info)

INPUT PARAMETERS
info -  sep_3d*   delete coord


RETURN VALUES
0   =  if deletion successful

-1  =  if deletion not successful


DESCRIPTION
Delete coordinate array of the dataset
>
*/

_XFUNCPROTOBEGIN
int SEP3D_delete_coord(sep_3d *info)
_XFUNCPROTOEND
{
int i1,ierr;

if(info->coord!=SEPNULL){
 free(info->coord); info->coord=SEPNULL;
}
info->ncoord=0;
return(SUCCESS);
}


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
SEP3D_alloc_coord

USAGE
ierr= SEP3D_alloc_coord(sep_3d *info, int ncoord)

INPUT PARAMETERS
info -  sep_3d*   structure associated with sep3d structure

ncoord    -  int     size to allocate array


RETURN VALUES
0   =  if allocation successful

-1  =  if allocation not successful


DESCRIPTION
Allocate coordinate array of the dataset
>
*/

_XFUNCPROTOBEGIN
int sep3d_alloc_coord(char *tag, int ncoord)
_XFUNCPROTOEND
{
sep_3d *info;
int i1,ierr;
info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
return(SEP3D_alloc_coord(info,ncoord));
}


_XFUNCPROTOBEGIN
int SEP3D_alloc_coord(sep_3d *info, int ncoord)
_XFUNCPROTOEND
{
int i1,ierr;
SEP3D_delete_coord(info);

info->ncoord=ncoord;
if(ncoord>0) info->coord=(long long*)malloc(sizeof(long long)*info->ncoord);

return(SUCCESS);
}

int SEP3D_coord_copy(sep_3d *input, sep_3d *output){
int i;

if(input->ncoord>0){
  if(0!=SEP3D_alloc_coord(output,input->ncoord))
    return(sepwarn(NOT_MET,"trouble allocating coord for tag %s \n",output->name));
  memcpy((void*)output->coord,(const void*)input->coord,sizeof(long long)*input->ncoord);



}
return(SUCCESS);
}
_XFUNCPROTOBEGIN
int sep3d_copy_coords(char *sep3din,char *sep3dout)
_XFUNCPROTOEND
{
int i1,ierr,idim,i;
sep_3d *in,*out;


in = tag_info_sep3d(sep3din, INQUIRE);  /* get info on this tag */
if(in == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3din));

out = tag_info_sep3d(sep3dout, INQUIRE);  /* get info on this tag */
if(out == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dout));

return(SEP3D_coord_copy(in,out));

}


_XFUNCPROTOBEGIN
int SEP3D_wind_coords(sep_3d *info,int *nsz)
_XFUNCPROTOEND
{
long long *blocks,*blockb,*tempd;
int idim,i,n123,ig,ic,nc,iax;
int *have,*mine;

for(i=1,n123=1; i < info->ndims; i++) n123=info->nwind[i]*n123;

if(n123<1){
  *nsz=0;
   SEP3D_delete_coord(info);
   return(0);
}
*nsz=n123;


blocks=(long long*)malloc(sizeof(long long)*(info->ndims-1));
blockb=(long long*)malloc(sizeof(long long)*(info->ndims-1));
blockb[0]=1;
for(idim=1; idim<info->ndims-1;idim++)
  blockb[idim]=blockb[idim-1]*info->n[idim];
blocks[0]=1;
                                                                                

for(idim=1; idim<info->ndims-1;idim++){
  blocks[idim]=blocks[idim-1]*info->nwind[idim];
}

 if(0!=SEP3D_alloc_coord(info,n123))
    return(sepwarn(NOT_MET,"trouble allocating coord for tag %s \n",
     info->name));
                                                                                
for(i=0; i <info->ncoord;i++) {
    info->coord[i]=0;
    for(idim=1,ig=0; idim< info->ndims; idim++){
      ic=SEPH2C(i,blocks[idim-1],info->nwind[idim]);

      info->coord[i]+=SEPC2H( (ic*info->jwind[idim]+info->fwind[idim]),blockb[idim-1]);
    }
}

free(blocks);free(blockb);

return(SUCCESS);
}
