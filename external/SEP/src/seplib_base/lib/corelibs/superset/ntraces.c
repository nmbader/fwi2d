#define SET_SDOC 1
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Sun Aug 16 15:04:02 PDT 1998

Purpose: 
	Handle ntraces element of superstructure

*/	 
#include "superset_internal.h"
#include <superset.h>


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
a
sep3d_grab_ntraces

USAGE
ierr= sep3d_grab_ntraces(char *sep3dname, int *ntraces)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep_3d structure

OUTPUT PARAMETERS
ntraces   -  int*     number of traces in the dataset



RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Grab the number of traces from the dataset
>
*/






_XFUNCPROTOBEGIN
int sep3d_grab_long_ntraces(char *sep3dname,long long *ntraces)
_XFUNCPROTOEND
{
int i1,ierr;
sep_3d *info;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

  return(SEP3D_grab_ntraces(info,ntraces));
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_ntraces(char *sep3dname,int *ntraces)
_XFUNCPROTOEND
#else
int sep3d_grab_ntraces(sep3dname, ntraces)
char *sep3dname;
int *ntraces;
#endif 
{
int i1,ierr;
sep_3d *info;
long long ntr;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

  ierr=SEP3D_grab_ntraces(info,&ntr);
  *ntraces=(int)ntr;
  return(ierr);
}


_XFUNCPROTOBEGIN
int SEP3D_grab_ntraces(sep_3d *info,long long *ntraces)
_XFUNCPROTOEND
{

/*if(0==info->ntraces) return(sepwarn(NOT_MET,*/
/*	"tag:%s number of traces not set\n",sep3dname));*/

*ntraces=info->ntraces;

return(SUCCESS);
}



/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_set_ntraces

USAGE
ierr= sep3d_set_ntraces(char *sep3dname, int ntraces)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep_3d structure
ntraces   -  int      number of traces in the dataset



RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Set the number of traces into the dataset
>
*/







_XFUNCPROTOBEGIN
int sep3d_set_long_ntraces(char *sep3dname,long long ntraces )
_XFUNCPROTOEND
{
int i1,ierr;
sep_3d *info;
long long big;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

 return(SEP3D_set_ntraces(info,ntraces));
}
_XFUNCPROTOBEGIN
int sep3d_set_big_ntraces(char *sep3dname,int nblock, int blocksize, int nextra )
_XFUNCPROTOEND
{
int i1,ierr;
sep_3d *info;
long long big;

big=(long long)nblock*(long long)blocksize+nextra;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

 return(SEP3D_set_ntraces(info,big));
}

_XFUNCPROTOBEGIN
int sep3d_set_ntraces(char *sep3dname,int ntraces)
_XFUNCPROTOEND
{
int i1,ierr;
sep_3d *info;
long long ntr;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

 ntr=ntraces;
 return(SEP3D_set_ntraces(info,ntr));
}

_XFUNCPROTOBEGIN
int SEP3D_set_ntraces(sep_3d *info,long long ntraces)
_XFUNCPROTOEND
{

info->ntraces=ntraces;
return(SUCCESS);
}

_XFUNCPROTOBEGIN
long long SEP3D_count_ntraces(sep_3d *info)
_XFUNCPROTOEND
{
int i,nt,ntot;
sep_3d *info_temp;

info->ntraces=info->ntraces_wrote;
return(SUCCESS);
}
_XFUNCPROTOBEGIN
int sep3d_count_ntraces(char *tag)
_XFUNCPROTOEND
{
sep_3d *info;
info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));

 return((int)SEP3D_count_ntraces(info));
}


/*  $Id: ntraces.c,v 1.2 2004/04/08 22:32:28 bob Exp $ */
