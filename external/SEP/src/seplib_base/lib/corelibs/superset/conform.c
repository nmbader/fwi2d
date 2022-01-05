#define SET_SDOC 1
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Tue Aug 27 08:40:24 PDT 2002


Purpose: 
	Handle ntraces element of superstructure

*/	 
#include "superset_internal.h"
#include <superset.h>

/*--------------------------------------------------------------------------*/
/*
<
sep3d_conform

USAGE
ierr= sep3d_conform(char *sep3dname1, char *sep3dname2)

INPUT PARAMETERS
sep3dname1 -  char*   tag associated with sep3d structure

sep3dname2 -  char*   tag associated with sep3d structure

RETURN VALUES

0   =  if they conform

5   =  if they don't conform (NOT_MET)

DESCRIPTION
Check to see if the spaces (o,d, n , ndim) are the same
>
*/


int sep3d_conform(char *sep3dname1,char  *sep3dname2)
{
int i1,ierr;
sep_3d *info1,*info2;


info1 = tag_info_sep3d(sep3dname1, INQUIRE);  /* get info on this tag */
if(info1 == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname1));

info2 = tag_info_sep3d(sep3dname2, INQUIRE);  /* get info on this tag */
if(info2 == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname2));

for(i1=0; i1 <MIN(info1->ndims,info2->ndims); i1++){
  if(info1->n[i1] != info2->n[i1]) return(NOT_MET);
  if(info1->o[i1] != info2->o[i1]) return(NOT_MET);
  if(info1->d[i1] != info2->d[i1]) return(NOT_MET);
}
  
return(SUCCESS);
}




/*--------------------------------------------------------------------------*/
/*
<
sep3d_ge_space

USAGE
ierr= sep3d_ge_space(char *sep3dname1, char *sep3dname2)

INPUT PARAMETERS
sep3dname1 -  char*   tag associated with sep3d structure

sep3dname2 -  char*   tag associated with sep3d structure

RETURN VALUES

0   =  if they conform

5   =  if they don't conform (NOT_MET)

DESCRIPTION
Check to see if the first space is greater or equal in all dimensions
>
*/


int sep3d_ge_space(char *sep3dname1,char  *sep3dname2)
{
int i1,ierr;
sep_3d *info1,*info2;
float beg1,end1,beg2,end2;


info1 = tag_info_sep3d(sep3dname1, INQUIRE);  /* get info on this tag */
if(info1 == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname1));

info2 = tag_info_sep3d(sep3dname2, INQUIRE);  /* get info on this tag */
if(info2 == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname2));

if(info1->ndims!= info2->ndims)  return(NOT_MET);
for(i1=0; i1 <info1->ndims; i1++){
  beg1=info1->o[i1]; end1=beg1+(info1->n[i1]-1)*info1->d[i1];
  beg2=info2->o[i1]; end2=beg2+(info2->n[i1]-1)*info2->d[i1];
  if(beg1 > beg2) return(NOT_MET);
  if(end2 > end1) return(NOT_MET);
}
  
return(SUCCESS);
}





