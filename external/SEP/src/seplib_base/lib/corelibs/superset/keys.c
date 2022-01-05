#define SET_SDOC 1
#include "superset_internal.h"
#include<string.h>
#include <superset.h>

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_set_nkeys

USAGE
ierr= sep3d_set_nkeys(char *sep3dname, int nkeys) 

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
nkeys     -  int     number of keys in the dataset


RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Sets the number of keys in the dataset

>
*/


/* this forward declaration belongs in an include file, not here. */
int sep3d_set_nkeys_fake(int,char *, int);

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_nkeys_fake(int x,char *sep3dname, int nkeys)
_XFUNCPROTOEND
#else
int sep3d_set_nkeys_fake(x,sep3dname ,nkeys)
char *sep3dname;
int nkeys,x;
#endif
{

return sep3d_set_nkeys(sep3dname,nkeys);

}
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_nkeys(char *sep3dname, int nkeys) 
_XFUNCPROTOEND
#else
int sep3d_set_nkeys(sep3dname ,nkeys) 
char *sep3dname;
int nkeys;
#endif 
{ 
int i1,isect;
sep_3d *info,*info_l;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));



if(nkeys==info->nkeys) return(SUCCESS);
if(info->nkeys!=0){ /*Deallocate old axis stuff */
	free(info->keyname); free(info->keyfmt); free(info->keytype); 
}

if(info->headers!=SEPNULL) 
	return(sepwarn(NOT_MET,
	"can not set nkeys without first clearing the header array  %s \n",sep3dname));

info->nkeys=nkeys;


info->keyname=(char **) malloc( nkeys * sizeof(char*));
info->keytype=(char **) malloc( nkeys * sizeof(char*));
info->keyfmt=(char  **) malloc( nkeys * sizeof(char*));

for(i1=0; i1 < nkeys; i1++){
	info->keyname[i1]=SEPNULL;
	info->keyfmt[i1]=SEPNULL;
	info->keytype[i1]=SEPNULL;
}
return(SUCCESS);

}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_nkeys

USAGE
ierr= sep3d_grab_nkeys(char *sep3dname, int *nkeys) 

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure

OUTPUT PARAMETERS
nkeys     -  int*     number of keys in the dataset


RETURN VALUES
0   =  if it works correctly

DESCRIPTION
Grabs the number of keys in the dataset 

>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_nkeys(char *sep3dname, int *nkeys) 
_XFUNCPROTOEND
#else
int sep3d_grab_nkeys(sep3dname ,nkeys) 
char *sep3dname;
int *nkeys;
#endif 
{ 
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));


*nkeys=info->nkeys;
return(SUCCESS);
}


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_set_key

USAGE
ierr= sep3d_set_key(char *sep3dname, int keyindex, char *keyname, char *keytype, char *keyformat) 

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
keyindex  -  int     key to set
keyname   -  char*   name of the key
keytype   -  char*   type of the key
keyfmt    -  char*   format of the key



RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and thier meaning)

DESCRIPTION
Sets the property of given key into the sep3d structure

>
*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_key(char *sep3dname, int keyindex, char *keyname, char *keytype, char *keyfmt) 
_XFUNCPROTOEND
#else
int sep3d_set_key(sep3dname ,keyindex,keyname,keytype,keyfmt)
char *sep3dname,*keyname, *keytype, *keyfmt;
int keyindex;
#endif 
{ 
sep_3d *info,*info_l;
int isect;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));




if(info->nkeys==0) return(sepwarn(NOT_MET,"nkeys not set (%s)\n", sep3dname));

if(keyindex <1 || keyindex > info->nkeys)
	 return(sepwarn(NOT_MET,
	"key is out of range 1 < index < %d (%s) \n",info->nkeys,sep3dname));



if(info->keyname[keyindex-1]!=SEPNULL) free(info->keyname[keyindex-1]);
if(info->keytype[keyindex-1]!=SEPNULL) free(info->keytype[keyindex-1]);
if(info->keyfmt[keyindex-1]!=SEPNULL) free(info->keyfmt[keyindex-1]);
/*
info->keyname[keyindex-1]=(char*) malloc((strlen(keyname)+1)*sizeof(char));
info->keytype[keyindex-1]=(char*) malloc((strlen(keytype)+1)*sizeof(char));
info->keyfmt[keyindex-1]=(char*) malloc((strlen(keyfmt)+1)*sizeof(char));
*/
info->keyname[keyindex-1]=(char*) malloc((SEP_3D_STRING_LEN+1)*sizeof(char));
info->keytype[keyindex-1]=(char*) malloc((SEP_3D_STRING_LEN+1)*sizeof(char));
info->keyfmt[keyindex-1]=(char*) malloc((SEP_3D_STRING_LEN+1)*sizeof(char));

strncpy(info->keyname[keyindex-1],keyname,SEP_3D_STRING_LEN);
info->keyname[keyindex-1][SEP_3D_STRING_LEN] = '\0';
strncpy(info->keyfmt[keyindex-1],keyfmt,SEP_3D_STRING_LEN);
info->keyfmt[keyindex-1][SEP_3D_STRING_LEN] = '\0';
strncpy(info->keytype[keyindex-1],keytype,SEP_3D_STRING_LEN);
info->keytype[keyindex-1][SEP_3D_STRING_LEN] = '\0';

	

return(SUCCESS);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_key

USAGE
ierr= sep3d_grab_key(char *sep3dname, int keyindex, char *keyname, char *keytype, char *keyformat) 

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
keyindex  -  int     key to set

OUTPUT PARAMETERS
keyname   -  char*   name of the key
keytype   -  char*   type of the key
keyfmt    -  char*   format of the key



RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Grabs the property of given key from the sep3d structure

>
*/





#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_key(char *sep3dname, int keyindex, char *keyname, char *keytype, char *keyfmt) 
_XFUNCPROTOEND
#else
int sep3d_grab_key(sep3dname ,keyindex,keyname,keytype,keyfmt)
char *sep3dname,*keyname, *keytype, *keyfmt;
int keyindex;
#endif 
{ 
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));



if(info->nkeys==0) return(sepwarn(NOT_MET,"nkeys not set (%s)\n", sep3dname));

if(keyindex <1 || keyindex > info->nkeys)
	 return(sepwarn(NOT_MET,
	"key is out of range 1 < index < %d (%s) \n",info->nkeys,sep3dname));


if(info->keyname[keyindex-1]==SEPNULL || info->keytype[keyindex-1]==SEPNULL ||
	info->keyfmt[keyindex-1]==SEPNULL) return(sepwarn(NOT_MET,
  "key %d has not been set in %s \n",keyindex,sep3dname));


strncpy(keyname,info->keyname[keyindex-1],SEP_3D_STRING_LEN);
strncpy(keytype,info->keytype[keyindex-1],SEP_3D_STRING_LEN);
strncpy(keyfmt,info->keyfmt[keyindex-1],SEP_3D_STRING_LEN);

return(SUCCESS);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_set_keys

USAGE
ierr= sep3d_set_keys(char *sep3dname,   char **keyname, char **keytype, char **keyfmt) 

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
keyname   -  char**  name of  keys
keytype   -  char**  types of keys
keyfmt    -  char**  units of keys



RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Sets keys properties into the sep3d structure

>
*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_keys(char *sep3dname,char **keyname,char **keytype, char **keyfmt) 
_XFUNCPROTOEND
#else
int sep3d_set_keys(sep3dname ,keyname,keytype,keyfmt)
char *sep3dname,**keyname,**keytype,**keyfmt;
#endif 
{
int i1,ierr;
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

if(info->nkeys==0) return(sepwarn(NOT_MET,"nkeys not set (%s)\n", sep3dname));

for(i1=0;i1 <  info->nkeys; i1++){
	ierr=sep3d_set_key(sep3dname,i1+1,keyname[i1],keytype[i1],keyfmt[i1]);
	if(ierr!=SUCCESS) return(ierr);
}

return(SUCCESS);
}


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_keys

USAGE
ierr= sep3d_grab_keys(char *sep3dname,   char **keyname, char **keytype, char **keyfmt) 

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure

OUTPUT PARAMETERS
keyname   -  char**  name of  keys
keytype   -  char**  types of keys
keyfmt    -  char**  units of keys



RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Grab keys properties from the sep3d structure

>
*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_keys(char *sep3dname,char **keyname,char **keytype, char **keyfmt) 
_XFUNCPROTOEND
#else
int sep3d_grab_keys(sep3dname ,keyname,keytype,keyfmt)
char *sep3dname,**keyname,**keytype,**keyfmt;
#endif 
{
int i1,ierr;
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

if(info->nkeys==0) return(sepwarn(NOT_MET,"nkeys not set (%s)\n", sep3dname));

for(i1=0;i1 <  info->nkeys; i1++){
	ierr=sep3d_grab_key(sep3dname,i1+1,keyname[i1],keytype[i1],keyfmt[i1]);
	if(ierr!=SUCCESS) return(ierr);
}

return(SUCCESS);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_add_key

USAGE
ierr= sep3d_add_key(char *sep3dname, char *keyname, char *keytype)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure
keyname   -  char*   keyname to add
keytype   -  char*   keytype to add



RETURN VALUES
0   =  if it works correctly

DESCRIPTION
Add key to sep3d structure
>
*/








#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_add_key(char *sep3dname,char *keyname, char *keytype)
_XFUNCPROTOEND
#else
int sep3d_add_key(sep3dname, keyname, keytype)
char *sep3dname,*keyname, *keytype;
#endif 
{
int i2,i1,ierr,nsz,isect;
sep_3d *info,*info_l;
char **tname,**tfmt,**ttype;
float *buffer;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));


if(info->nkeys>0){
	tname=(char**)malloc(info->nkeys*sizeof(char*));
	tfmt=(char**)malloc(info->nkeys*sizeof(char*));
	ttype=(char**)malloc(info->nkeys*sizeof(char*));
	for(i1=0; i1< info->nkeys; i1++){
		tname[i1]=info->keyname[i1];
		tfmt[i1]=info->keyfmt[i1];
		ttype[i1]=info->keytype[i1];
	}
	if(0!=info->headers){
		buffer=(float *)malloc(info->nh*info->nkeys*sizeof(float));
		for(i1=0; i1 < info->nkeys*info->nh;i1++) buffer[i1]=info->headers[i1];;
	}
}
	ierr=sep3d_set_nkeys(sep3dname,info->nkeys+1);
	if(ierr!=SUCCESS) return(ierr);

if(info->nkeys>1){
	for(i1=1; i1<info->nkeys;i1++){
		ierr=sep3d_set_key(sep3dname,i1,tname[i1-1],ttype[i1-1],tfmt[i1-1]);
		if(ierr!=SUCCESS)		return(ierr);
	}
	free(tname);free(ttype);free(tfmt);
	if(info->headers!=0){
		free(info->headers);
		for(i2=0; i2< info->nh; i2++){
			for(i1=0; i1 < info->nkeys-1; i1++){
				info->headers[i1+i2*info->nkeys]=buffer[i1+i2*(info->nkeys-1)];
			}
				info->headers[(i2+1)*info->nkeys-1]=0.;
		}
	 free(buffer);
	}
}
if(0==strcmp("scalar_float",keytype)){
		ierr=sep3d_set_key(sep3dname,info->nkeys,keyname,keytype, "xdr_float");
		if(ierr!=0) return(ierr);
}else if(0==strcmp("scalar_int",keytype)){
		ierr=sep3d_set_key(sep3dname,info->nkeys,keyname,keytype, "xdr_int");
		if(ierr!=0) return(ierr);
}
else
	return(sepwarn(INVALID_DATA,"Keytype specified unrecognized \n"));



return(SUCCESS);
}


/*<
sep3d_grab_key_index

USAGE
int*= sep3d_grab_key_index(char *tag, char keyname, int *keyindex)

INPUT PARAMETERS
tag       -  char*   tag of sep3d structure
keyname   -  char*   name of key
keyindex  -  int*    index of key


RETURN VALUES
0   =  if key is found
(check superset.h for other return values and there meaning)

Nullifys a header pointer
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_key_index(char *tag,char *keyname, int *keyindex)
_XFUNCPROTOEND
#else
int sep3d_grab_key_index(tag,keyname,keyindex)
char *tag,*keyname;
int *keyindex;
#endif 
{
sep_3d *info;
int i,tot,ret;
float one;

info = tag_info_sep3d(tag, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",tag));

if(info->nkeys==0) return(sepwarn(NOT_MET,"nkeys not set (%s)\n", tag));

i=0; ret=1;
while(i< info->nkeys && ret !=0){
	ret=strcmp(info->keyname[i],keyname);
	i++;
}

*keyindex=i;
if(ret!=0 ) return(NOT_MET);

return(SUCCESS);
}
