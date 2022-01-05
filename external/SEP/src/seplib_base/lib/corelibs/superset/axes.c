#define SET_SDOC 1
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Sun Aug 16 15:04:02 PDT 1998

Purpose: To handle axis information in the superset data

*/	 

#include<string.h>
#include "superset_internal.h"
#include <superset.h> 

#include <ctype.h>
static int fstrlen(const char *s, int maxlen)
{
  int i;
  for(i=(maxlen-1); i>=0; i--) {
     if((!isspace(s[i])) && s[i] != '\0') break;
  }
  return i+1;
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*<
sep3d_set_ndims

USAGE
ierr= sep3d_set_ndims(char *sep3dname, int ndims)


INPUT PARAMETERS
sep3dname    -  char*  tag of sep3d structure
ndims        -  int    number of dimension in the dataset

RETURN VALUES
0  =  works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Sets the number of dimesions in the dataset into the sep3d structure

>*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_ndims(char *sep3dname, int ndims) 
_XFUNCPROTOEND
#else
int sep3d_set_ndims(sep3dname ,ndims) 
char *sep3dname;
int ndims;
#endif 
{ 
int i1;
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));
return(SEP3D_set_ndims(info,ndims));
}

_XFUNCPROTOBEGIN
int SEP3D_set_ndims(sep_3d *info, int ndims)
_XFUNCPROTOEND
{
int i1;


if(ndims==info->ndims)  return(0);

 if(0!=SEP3D_del_axes(info))
   return (sepwarn(NOT_MET,"tag:%s  trouble deleting axis info \n",info->name));


info->ndims=ndims;

info->nwind=(int *) malloc( info->ndims * sizeof(int));
info->fwind=(int *) malloc( info->ndims * sizeof(int));
info->jwind=(int *) malloc( info->ndims * sizeof(int));

info->n=(int *) malloc( info->ndims * sizeof(int));
info->o=(float *) malloc( info->ndims * sizeof(float));
info->d=(float *) malloc( info->ndims * sizeof(float));
info->label=(char **) malloc( info->ndims * sizeof(char*));
info->unit=(char **) malloc( info->ndims * sizeof(char*));

for(i1=0; i1 < info->ndims; i1++){
	info->n[i1]=1;
	info->o[i1]=0;
	info->d[i1]=1;
  info->nwind[i1]=1;
  info->fwind[i1]=0;
  info->jwind[i1]=1;
	info->label[i1]=(char *) malloc(sizeof(char)* 256);
	info->unit[i1]=(char *) malloc(sizeof(char)* 256);
	strcpy(info->label[i1],"Undefined");
	strcpy(info->unit[i1],"Undefined");
}
return(SUCCESS);
}

int  SEP3D_del_axes(sep_3d *info){
int i1;


if(info->ndims!=0){ /*Deallocate old axis stuff */
	free(info->n); free(info->o); free(info->d); 
	free(info->nwind); free(info->fwind); free(info->jwind); 
	for(i1=0; i1 < info->ndims; i1++){
		if(info->label[i1]!=SEPNULL){  free(info->label[i1]);}
		if(info->unit[i1]!=SEPNULL) { free(info->unit[i1]);}
	}
  if(info->label!=SEPNULL){ free(info->label); }
  if(info->unit!=SEPNULL){ free(info->unit); }
	
}
info->ndims=0;


return(SUCCESS);
}



/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_ndims

USAGE
ierr= sep3d_grab_ndims(char *sep3dname, int *ndims)


INPUT PARAMETERS
sep3dname    -  char*  tag of sep3d structure


OUTPUT PARAMETERS
ndims        -  int*    number of dimension in the dataset

RETURN VALUES
0 = if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Gets the number of axes in the dataset from the associated structure

>
*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_ndims(char *sep3dname, int *ndims) 
_XFUNCPROTOEND
#else
int sep3d_grab_ndims(sep3dname ,ndims) 
char *sep3dname;
int *ndims;
#endif 
{ 
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));
return(SEP3D_grab_ndims(info,ndims));
}

_XFUNCPROTOBEGIN
int SEP3D_grab_ndims(sep_3d *info, int *ndims) 
_XFUNCPROTOEND
{
*ndims=info->ndims;
return(SUCCESS);
}

 
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_set_axis

USAGE
ierr= sep3d_set_axis(char *sep3dname, int axis,  int n, float o, float d, char *label, char *unit) 

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
axis      -  int     axis to set
n         -  int     length of given axis
o         -  float   first sample of given axis
d         -  float   sampling of given axis
label     -  char*   label of given axis
unit      -  char*   unit of given axis



RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Sets axis properties in the sep3d structure

>
*/


 
 
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_axis(char *sep3dname, int axis, int n, float o, float d,char *label, char *unit) 
_XFUNCPROTOEND
#else
int sep3d_set_axis(sep3dname ,axis,n,o,d,label,unit)
char *sep3dname,*label,*unit;
int axis,n;
float o,d;
#endif 
{ 
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

return(SEP3D_set_axis(info,axis,n,o,d,label,unit));
}

_XFUNCPROTOBEGIN
int SEP3D_set_axis(sep_3d *info, int axis, int n, float o, float d,char *label, char *unit) 
_XFUNCPROTOEND
{
int isect,i;
sep_3d *info_l;

if(info->ndims==0) return(sepwarn(NOT_MET,
 "set_axis:Must set number of dimensions :%s \n",info->name));

if(axis <1 || axis > info->ndims)
	return(sepwarn(NOT_MET,
	"%s: Axis is out of range 1 < index < %d : \n",info->name,info->ndims,axis));


if((int)strlen(unit)>(int)strlen(info->unit[axis-1])){
	free(info->unit[axis-1]);
/*	info->unit[axis-1]=(char*) malloc((strlen(unit)+1)*sizeof(char)); */
	info->unit[axis-1]=(char*) malloc((SEP_3D_STRING_LEN+1)*sizeof(char));
	strncpy(info->unit[axis-1],unit,SEP_3D_STRING_LEN);
        info->unit[axis-1][SEP_3D_STRING_LEN] = '\0';
}
else strncpy(info->unit[axis-1],unit,SEP_3D_STRING_LEN);
if(NULL!=label){
if((int)strlen(label)>(int)strlen(info->label[axis-1])){
	free(info->label[axis-1]);
/*	info->label[axis-1]=(char*) malloc((strlen(label)+1)*sizeof(char)); */
	info->label[axis-1]=(char*) malloc((SEP_3D_STRING_LEN+1)*sizeof(char));
	strncpy(info->label[axis-1],label,SEP_3D_STRING_LEN);
        info->label[axis-1][SEP_3D_STRING_LEN] = '\0';
}
else strncpy(info->label[axis-1],label,SEP_3D_STRING_LEN);
}


info->n[axis-1]=n;
for(i=1,info->ntraces=1; i < info->ndims; i++)  info->ntraces=info->ntraces*info->n[i];
info->nwind[axis-1]=n;
info->o[axis-1]=o;
info->d[axis-1]=d;

return(SUCCESS);
}


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_axis

USAGE
ierr= sep3d_grab_axis(char *sep3dname, int axis,  int *n, float *o, float *d, char *label, char *unit) 

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
axis      -  int     axis to set

OUTPUT PARAMETERS
n         -  int     length of given axis
o         -  float   first sample of given axis
d         -  float   sampling of given axis
label     -  char*   label of given axis
unit      -  char*   unit of given axis



RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Grab axis properties from the sep3d structure

>
*/




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_axis(char *sep3dname, int axis, int *n, float *o, float *d,char *label, char *unit) 
_XFUNCPROTOEND
#else
int sep3d_grab_axis(sep3dname ,axis,n,o,d,label,unit)
char *sep3dname,*label,*unit;
int axis,*n;
float *o,*d;
#endif 
{ 
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

return(SEP3D_grab_axis(info,axis,n,o,d,label,unit));
}

_XFUNCPROTOBEGIN
int SEP3D_grab_axis(sep_3d *info, int axis, int *n, float *o, float *d,char *label, char *unit) 
_XFUNCPROTOEND
{

if(info->ndims==0) return(sepwarn(NOT_MET,
 "%s:grab_axis:Must set number of dimensions \n",info->name));
if(axis <1 || axis > info->ndims) return(sepwarn(NOT_MET,
	"%s: Axis is out of range 1 < index < %d:%d \n",info->name,info->ndims,axis));

if(info->unit[axis-1]==SEPNULL) strcpy(unit,"Unspecified");
else strncpy(unit,info->unit[axis-1],SEP_3D_STRING_LEN);
if(info->label[axis-1]==SEPNULL) strcpy(label,"Unspecified");
else strncpy(label,info->label[axis-1],SEP_3D_STRING_LEN);


*n=info->n[axis-1];
*o=info->o[axis-1];
*d=info->d[axis-1];

return(SUCCESS);
}


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_nod

USAGE
ierr= sep3d_grab_nod(char *sep3dname,   int *n, float *o, float *d)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure

OUTPUT PARAMETERS
n         -  int*    length of given axis
o         -  float*  first sample of given axis
d         -  float*  sampling of given axis


RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Grab axes properties from the sep3d structure

>
*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_nod(char *sep3dname, int *n, float *o, float *d)
_XFUNCPROTOEND
#else
int sep3d_grab_axes(sep3dname ,axis,n,o,d)
char *sep3dname;
int axis,*n;
float *o,*d;
#endif 
{
int i1,ierr;
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

  return(SEP3D_grab_nod(info,n,o,d));
}

_XFUNCPROTOBEGIN
int SEP3D_grab_nod(sep_3d *info, int *n, float *o, float *d)
_XFUNCPROTOEND
{
int i1;

if(info->ndims==0)
	return(sepwarn(NOT_MET,
	"%s:grab_axes:Must set number of dimensions \n",info->name));


for(i1=0;i1 <  info->ndims; i1++){
  n[i1]=info->n[i1];
  o[i1]=info->o[i1];
  d[i1]=info->d[i1];
}

return(SUCCESS);
}
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_axes

USAGE
ierr= sep3d_grab_axes(char *sep3dname,   int *n, float *o, float *d, char **label, char **unit) 

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure

OUTPUT PARAMETERS
n         -  int*    length of given axis
o         -  float*  first sample of given axis
d         -  float*  sampling of given axis
label     -  char**  label of given axis
unit      -  char**  unit of given axis



RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Grab axes properties from the sep3d structure

>
*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_axes(char *sep3dname, int *n, float *o, float *d,char **label, char **unit) 
_XFUNCPROTOEND
#else
int sep3d_grab_axes(sep3dname ,axis,n,o,d,label,unit)
char *sep3dname,**label,**unit;
int axis,*n;
float *o,*d;
#endif 
{
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

  return(SEP3D_grab_axes(info,n,o,d,label,unit));
}
_XFUNCPROTOBEGIN
int SEP3D_grab_axes(sep_3d *info, int *n, float *o, float *d,char **label, char **unit) 
_XFUNCPROTOEND
{
int i1,ierr;

if(info->ndims==0)
	return(sepwarn(NOT_MET,
	"%s:grab_axes:Must set number of dimensions \n",info->name));


for(i1=1;i1 <=  info->ndims; i1++){
	ierr=SEP3D_grab_axis(info,i1,&n[i1],&o[i1],&d[i1],label[i1],unit[i1]);
	if(ierr!=0) return(ierr);
}

return(SUCCESS);
}



/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_set_axes

USAGE
ierr= sep3d_set_axes(char *sep3dname,   int *n, float *o, float *d, char ***label, char **unit) 

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
n         -  int*    length of given axis
o         -  float*  first sample of given axis
d         -  float*  sampling of given axis
label     -  char**  label of given axis
unit      -  char**  unit of given axis



RETURN VALUES
0   =  if it works correctly

DESCRIPTION
Sets axes properties into the sep3d structure

>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_axes(char *sep3dname, int *n, float *o, float *d,char **label, char **unit) 
_XFUNCPROTOEND
#else
int sep3d_set_axes(sep3dname ,axis,n,o,d,label,unit)
char *sep3dname,**label,**unit;
int axis,*n;
float *o,*d;
#endif 
{
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

  return(SEP3D_set_axes(info,n,o,d,label,unit));

}
_XFUNCPROTOBEGIN
int SEP3D_set_axes(sep_3d *info, int *n, float *o, float *d,char **label, char **unit) 
_XFUNCPROTOEND
{

int i1,ierr;
if(info->ndims==0) return(sepwarn(NOT_MET,
	"%s:set_axes:Must set number of dimensions \n",info->name));

for(i1=0;i1 <  info->ndims; i1++){
	ierr=SEP3D_set_axis(info,i1+1,n[i1],o[i1],d[i1],label[i1],unit[i1]);
	if(ierr!=0) return(ierr);
}

return(SUCCESS);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_rite_axis

USAGE
ierr= sep3d_rite_axis(char *sep3dname, int axis, int n, float o, float d char *label, char *unit)

INPUT PARAMETERS
tag       -  char*   tag of history, header, or grid
axis      -  int     axis 
n         -  int     length of axis
o         -  float   first sample of axis
d         -  float   sampling of axis
label     -  char*   label for axis
unit      -  char*   unit for axis


RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Write axis information 
>
*/





#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_rite_axis(char *tag, int axis, int n, float o, float d, char *label, char *unit)
_XFUNCPROTOEND
#else
int sep3d_rite_axis(tag,axis,n, o,d, label,unit)
char *tag;int axis; int n; float o; float d; char *label, *unit;
#endif 
{
int ierr;
char temp_ch[2048],temp2_ch[1024];


sprintf(temp_ch,"\tn%d=%d  o%d=%f  d%d=%f",
axis,n,axis,o,axis,d);
if((int)strlen(label)>0){ sprintf(temp2_ch,"%s   label%d=\"%.*s\"",temp_ch,
  axis, fstrlen(label,SEP_3D_STRING_LEN), label);  strcpy(temp_ch,temp2_ch);}
if((int)strlen(unit)>0){ sprintf(temp2_ch,"%s   unit%d=\"%.*s\"",temp_ch,
  axis, fstrlen(unit,SEP_3D_STRING_LEN),unit);  strcpy(temp_ch,temp2_ch);}
if(0!=auxputhead(tag,temp_ch))
  return(sepwarn(FAIL_OTHER,"could not puthead to tag %s \n",tag));
return(SUCCESS);
}



/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_reed_axis

USAGE
ierr= sep3d_reed_axis(char *sep3dname, int axis, int *n, float *o, float *d char *label, char *unit)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
axis      -  int     axis 

OUTPUT PARAMETERS
n         -  int     length of axis
o         -  float*  first sample of axis
d         -  float*  sampling of axis
label     -  char*   label for axis
unit      -  char*   unit for axis


RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Read axis information 
>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_reed_axis(char *tag, int axis, int *n, float *o, float *d, char *label, char *unit)

_XFUNCPROTOEND
#else
int sep3d_reed_axis(tag,axis,n, o,d, label,unit)
char *tag;int axis; int *n; float *o; float *d; char *label, *unit;
#endif 
{
int ierr;
char temp_ch[512];

sprintf(temp_ch,"n%d",axis);
ierr=auxpar(temp_ch,"d",n,tag);
if(ierr<0)return(sepwarn(ierr,"auxpar failed on tag %s for %s \n",temp_ch,tag));
if(ierr==0) *n=1;

sprintf(temp_ch,"o%d",axis);
ierr=auxpar(temp_ch,"f",o,tag);
if(ierr<0)return(sepwarn(ierr,"auxpar failed on tag %s for %s \n",temp_ch,tag));
if(ierr==0) *o=0;

sprintf(temp_ch,"d%d",axis);
ierr=auxpar(temp_ch,"f",d,tag);
if(ierr<0)return(sepwarn(ierr,"auxpar failed on tag %s for %s \n",temp_ch,tag));
if(ierr==0) *d=1;

sprintf(temp_ch,"label%d",axis);
ierr=auxpar(temp_ch,"s",label,tag);
if(ierr<0)return(sepwarn(ierr,"auxpar failed on tag %s for %s \n",temp_ch,tag));
if(ierr==0) strcpy(label,"Undefined");

sprintf(temp_ch,"unit%d",axis);
ierr=auxpar(temp_ch,"s",unit,tag);
if(ierr<0)return(sepwarn(ierr,"auxpar failed on tag %s for %s \n",temp_ch,tag));
if(ierr==0) strcpy(unit,"Undefined");
return(SUCCESS);
}

/*<
sep3d_grab_axis_index

USAGE
int*= sep3d_grab_axis_index(char *tag, char axisname, int *axisindex)

INPUT PARAMETERS
tag       -  char*   tag of sep3d structure
axisname   -  char*   name of key
axisindex  -  int*    index of key


RETURN VALUES
0   =  if axis is found
(check superset.h for other return values and there meaning)

DESCRIPTION
Find the index of a given axis (if it exists)
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_wind(char *tag, int *nwind, int *fwind, int *jwind)
_XFUNCPROTOEND
#endif
{
sep_3d *info;
int i,tot,ret;
float one;

info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));

if(info->ndims==0) return(sepwarn(NOT_MET,"ndims not set (%s)\n", tag));

for(i=0; i < info->ndims; i++){
  info->nwind[i]=nwind[i];
  info->fwind[i]=fwind[i];
  info->jwind[i]=jwind[i];
}

return(SUCCESS);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_wind(char *tag, int *nwind, int *fwind, int *jwind)
_XFUNCPROTOEND
#endif
{
sep_3d *info;
int i,tot,ret;
float one;

info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));

if(info->ndims==0) return(sepwarn(NOT_MET,"ndims not set (%s)\n", tag));

for(i=0; i < info->ndims; i++){
  nwind[i]=info->nwind[i];
  fwind[i]=info->fwind[i];
  jwind[i]=info->jwind[i];
}

return(SUCCESS);
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_axis_index(char *tag,char *axisname, int *axisindex)
_XFUNCPROTOEND
#else
int sep3d_grab_axis_index(tag,axisname,axisindex)
char *tag,*axisname;
int *axisindex;
#endif
{
sep_3d *info;
int i,tot,ret;
float one;

info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));

if(info->ndims==0) return(sepwarn(NOT_MET,"ndims not set (%s)\n", tag));

i=0; ret=1;
while(i< info->ndims && ret !=0){
	if(info->label[i]!=SEPNULL)
  	ret=strcmp(info->label[i],axisname);
  	i++;
}

*axisindex=i;
if(ret!=0 ) return(NOT_MET);

return(SUCCESS);
}



/*<
sep3d_change_axes_index

USAGE
int*= sep3d_change_axes_index(char *tag, int naxes, int *axisindex)

INPUT PARAMETERS
tag       -  char*   tag of sep3d structure
naxes     -  int      number of axes in new dataset
axisindex  -  int*    index of key


RETURN VALUES
0   =  if axis is found
(check superset.h for other return values and there meaning)

DESCRIPTION
Modifies the number of axes in sep3d dataset eg

3-D can become 2-D,etc

>
*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_change_dims(char *tag,int naxes, int *axisindex)
_XFUNCPROTOEND
#endif
{
sep_3d *info;

info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));
return(SEP3D_change_dims(info,naxes,axisindex));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int SEP3D_change_dims(sep_3d *info,int naxes, int *axisindex)
_XFUNCPROTOEND
#endif
{
int i,ret,ierr,old_nax,nax,reshape[41];
float one;
int n[40];
float o[40], d[40];
char label[40][256],unit[40][256];
int nout,last_axes,tot,iout,ndim;
sep_3d *info_l;


if(info->ndims==0) return(sepwarn(NOT_MET,"ndims not set (%s)\n", info->name));
if(info->ndims>40) return(sepwarn(NOT_MET,"ndims to large \n"));

for(i=0; i < naxes; i++) reshape[i]=axisindex[i];
nax=naxes;
if(axisindex[naxes-1]!=info->ndims){
  for(i=0; i < info->ndims-axisindex[naxes-1]; i++,nax++) {
    reshape[i+naxes]=axisindex[naxes-1]+i+1;
  }
}

ndim=info->ndims;
for(i=0; i <40; i++) n[i]=1;

/*first read in the old axes */
for(i=1; i <= info->ndims; i++){
	ierr=sep3d_grab_axis(info->name,i,&n[i-1],&o[i-1],&d[i-1],label[i-1],unit[i-1]);
	if(ierr!=SUCCESS) return(ierr);
}

ierr=sep3d_set_ndims(info->name,nax);
if(ierr!=SUCCESS) return(ierr);

last_axes=0;
tot=1;
for(i=0; i<nax; i++){
	if(reshape[i]< last_axes)
		return(sepwarn(INVALID_DATA,"sep3d_change_dims:axes must be increasing \n"));
/*		return(sepwarn(INVALID_DATA,"sep3d_change_dims:axis(%d) larger than input ndim(%d) tag=%s \n",reshape[i],ndim,tag));*/
	if(reshape[i]==last_axes){ /* this axis just defaults */
   ierr=sep3d_set_axis(info->name,i+1,1,0.,1.,"none","none");
  }
	else if(reshape[i]==0){ /*we are creating a fake first axis */
  }
	else {
		nout=1;
		for(iout=last_axes;iout<reshape[i];iout++) nout=nout*n[iout];
ierr=sep3d_set_axis(info->name,i+1,nout,o[reshape[i]-1],d[reshape[i]-1],
     label[reshape[i]-1], unit[reshape[i]-1]);
		if(ierr!=SUCCESS) return(ierr);
		last_axes=reshape[i];
	}
}
		
		

if(reshape[nax-1]<ndim) return(sepwarn(INVALID_DATA,
 "output axes (%d) does not span input dimensions(%d) tag=%s \n",reshape[naxes-1],ndim,info->name));

return(SUCCESS);

}




/*  $Id: axes.c,v 1.3 2004/07/08 18:15:32 bob Exp $ */
