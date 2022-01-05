/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Mon Aug 17 14:19:08 PDT 1998

Purpose: 

*/	 

#include <sepConfig.h>

#ifdef HAVE_STRING_H
#include <string.h>
#endif
#define SET_SDOC 3
#include <superset.h> 
#include <superset_internal.h> 
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/* 
<
sep3d_grab_header_vals_s

USAGE
ierr=sep3d_grab_header_vals_s(sep3dname,key,values)

INPUT PARAMETERS
sep3dname   -  char*     pointer to sep3d structure
keyname     -  char*     name of the key to get


OUTPUT PARAMETERS
values      -  int*    values for the key

RETURN VALUES
0   = if it works

DESCRIPTION
Gets key values form sep3d structure

>
*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_header_vals_s(char *sep3dname, char *key, int *values) 
_XFUNCPROTOEND
#else
int sep3d_grab_header_vals_s(sep3dname,key,values) 
char *sep3dname, *key;
int *values;
#endif 
{ 
sep_3d *info;
int i1,kindex,ierr,logic;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == 0) return (INVALID_STRUC);   /* Not a valid */

if(info->nkeys==0) return(sepwarn(NOT_MET,
  "key array has not been set up for tag %s \n",sep3dname));

kindex=0;logic=0;i1=0;
while(logic==0 && i1 < info->nkeys){
	if(0==strcmp(key,info->keyname[i1])){
		logic=1;
		return(sep3d_grab_header_vals_i(sep3dname,i1+1,values));
	}
	i1++;
}
if(logic==0) return(sepwarn(NOT_MET,"key %s not initialized for tag %s \n",key,
  sep3dname));
return(SUCCESS);
}


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/

/* 
<
sep3d_grab_header_vals_i

USAGE
ierr=sep3d_grab_header_vals_i(sep3dname,kindex,values)

INPUT PARAMETERS
sep3dname   -  char*     pointer to sep3d structure
kindex      -  int       index of the key to get


OUTPUT PARAMETERS
values      -  int*    values for the key

RETURN VALUES
0   = if it works

DESCRIPTION
Gets key values form sep3d structure

>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_header_vals_i(char *sep3dname, int kindex, int *values) 
_XFUNCPROTOEND
#else
int sep3d_grab_header_vals_i(sep3dname,kindex,values) 
char *sep3dname;
int *values;
int kindex;
#endif 
{ 
sep_3d *info;
int i1;
int n;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == 0) return (INVALID_STRUC);   /* Not a valid */

if(info->nkeys==0) return(sepwarn(NOT_MET,
  "key array has not been set up for tag  %s \n",sep3dname));
if(info->headers==SEPNULL) return(sepwarn(NOT_MET,
	"header array has not been set up for tag %s \n",sep3dname));
if(info->nkeys < kindex || kindex <1) 
	return(sepwarn(NOT_MET,
	"key index out of range 1 < %d < %d :%s \n",kindex,info->nkeys,sep3dname));


for(i1=0; i1< info->nh; i1++){
	 values[i1]=info->headers[kindex-1+i1*info->nkeys];
}
return (SUCCESS);
} 


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/* 
<
sep3d_set_header_vals_s

USAGE
ierr=sep3d_set_header_vals_s(sep3dname,key,values)

INPUT PARAMETERS
sep3dname   -  char*     pointer to sep3d structure
keyname     -  char*     name of the key to get
values      -  int*    values for the key

RETURN VALUES
0   = if it works

DESCRIPTION
Puts key values into sep3d structure

>
*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_header_vals_s(char *sep3dname, char *key, int *values) 
_XFUNCPROTOEND
#else
int sep3d_et_header_vals_s(sep3dname,key,values) 
char *sep3dname, *key;
int *values;
#endif 
{ 
sep_3d *info;
int i1,kindex,ierr,logic;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL) return (sepwarn(INVALID_STRUC,
	"%s: struct has not been set up \n",sep3dname));   /* Not a valid */

if(info->nkeys==0)return(sepwarn(NOT_MET,"key array has not been set up %s \n",
   sep3dname));


kindex=0;logic=0;i1=0;
while(logic==0 && i1 < info->nkeys){
	if(0==strcmp(key,info->keyname[i1])){
		logic=1;
		return(sep3d_set_header_vals_i(sep3dname,i1+1,values));
	}
	i1++;
}
if(logic==0) return(sepwarn(NOT_MET,"key %s not initialized in %s \n",
  key,sep3dname));

return(SUCCESS);
}



/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/* 
<
sep3d_set_header_vals_i

USAGE
ierr=sep3d_set_header_vals_i(sep3dname,kindex,values)

INPUT PARAMETERS
sep3dname   -  char*     pointer to sep3d structure
kindex      -  int       index of the key to get
values      -  int*    values for the key

RETURN VALUES
0   = if it works

DESCRIPTION
Puts key values into sep3d structure

>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_header_vals_i(char *sep3dname, int kindex, int *values) 
_XFUNCPROTOEND
#else
int sep3d_set_header_vals_i(sep3dname,kindex,values) 
char *sep3dname;
int *values;
int kindex;
#endif 
{ 
sep_3d *info;
int i1;
int n;
float *fval;
fprintf(stderr,"In set header \n");

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL) /* Not a valid */
	return (sepwarn(INVALID_STRUC,
  "Tag has not been initialized %s \n",sep3dname));   

if(info->nkeys==0) return(sepwarn(NOT_MET,
  "key array has not been set up for tag  %s \n",sep3dname));
if(info->headers==0) return(sepwarn(NOT_MET,
	"header array has not been set up for tag %s \n",sep3dname));
if(info->nkeys < kindex || kindex <1) 
	return(sepwarn(NOT_MET,
	"key index out of range 1 < %d < %d :%s \n",kindex,info->nkeys,sep3dname));

for(i1=0; i1< info->nh; i1++) info->headers[kindex-1+i1*info->nkeys]=values[i1];
return (SUCCESS);
} 



/*  $Id: headval.c,v 1.2 2004/04/08 22:32:28 bob Exp $ */
