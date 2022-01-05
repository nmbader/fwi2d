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

extern int sep_thread_num(void);
extern int SEP3D_grab_file_type(sep_3d *info,char *file_type);

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_usage

USAGE
ierr= sep3d_grab_usage(char *sep3dname,   char *usage)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep_3d structure

OUTPUT PARAMETERS
usage     -  char*   usage for the dataset (INPUT, OUTPUT, or SCRATCH)



RETURN VALUES
0   =  if it works correctly

DESCRIPTION
(check superset.h for other return values and there meaning)

>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_usage(char *sep3dname,char *usage)
_XFUNCPROTOEND
#else
int sep3d_grab_usage(sep3dname ,usage)
char *sep3dname,*usage;
#endif 
{
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc thread=%d\n",sep3dname,sep_thread_num()));

 return(SEP3D_grab_usage(info,usage));
}
_XFUNCPROTOBEGIN
int SEP3D_grab_usage(sep_3d *info,char *usage)
_XFUNCPROTOEND

{

int i1,ierr;


switch(info->usage){
	case(INPUT): strcpy(usage,"INPUT"); break;
	case(OUTPUT): strcpy(usage,"OUTPUT"); break;
	case(SCRATCH): strcpy(usage,"SCRATCH"); break;
}
return(SUCCESS);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_file_type

USAGE
ierr= sep3d_grab_file_type(char *sep3dname,   char *file_type)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep_3d structure

OUTPUT PARAMETERS
file_type -  char*   file type for the dataset (REGULAR, HEADER, or GRID)



RETURN VALUES
0   =  if it works correctly

DESCRIPTION
Grab the file type of the dataset
>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_file_type(char *sep3dname,char *file_type)
_XFUNCPROTOEND
#else
int sep3d_grab_file_type(sep3dname ,file_type)
char *sep3dname,*file_type;
#endif 
{
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

return(SEP3D_grab_file_type(info,file_type));
}
_XFUNCPROTOBEGIN
int SEP3D_grab_file_type(sep_3d *info,char *file_type)
_XFUNCPROTOEND
{

int i1,ierr;
switch(info->file_format){
	case(REGULAR): strcpy(file_type,"REGULAR"); break;
	case(HEADER): strcpy(file_type,"HEADER"); break;
	case(GRID): strcpy(file_type,"GRID"); break;
	case(UNSPECIFIED): strcpy(file_type,"UNSPECIFIED"); break;
}
return(SUCCESS);
}




/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_set_file_type

USAGE
ierr= sep3d_set_file_type(char *sep3dname,   char *file_type)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep_3d structure
file_type -  char*   file type for the dataset (REGULAR, HEADER, or GRID)



RETURN VALUES
0   =  if it works correctly

DESCRIPTION
Sets the file type into the dataset
>
*/




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_file_type(char *sep3dname,char *file_type)
_XFUNCPROTOEND
#else
int sep3d_set_file_type(sep3dname ,file_type)
char *sep3dname,*file_type;
#endif 
{
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

  return(SEP3D_set_file_type(info,file_type));
}

_XFUNCPROTOBEGIN
int SEP3D_set_file_type(sep_3d *info,char *file_type)
_XFUNCPROTOEND
{
	if(0==strcmp("REGULAR",file_type)) info->file_format=REGULAR;
	else if(0==strcmp("HEADER",file_type)) info->file_format=HEADER;
	else if(0==strcmp("GRID",file_type)) info->file_format=GRID;
	else seperr("Unknown file type %s for sep3dtag=%s \n",file_type,info->name);
  return(SUCCESS);
}





/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_grab_data_type

USAGE
ierr= sep3d_grab_data_type(char *sep3dname,   char *data_type)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep_3d structure
data_type -  char*   file type for the dataset (FLOAT, INTEGER, COMPLEX, or BYTE)



RETURN VALUES
0   =  if it works correctly

DESCRIPTION
Sets the data type into the dataset
>
*/






#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_grab_data_type(char *sep3dname,char *data_type)
_XFUNCPROTOEND
#else
int sep3d_grab_data_type(sep3dname ,data_type)
char *sep3dname,*data_type;
#endif 
{
int i1,ierr;
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

return(SEP3D_grab_data_type(info,data_type));
}

_XFUNCPROTOBEGIN
int SEP3D_grab_data_type(sep_3d *info,char *data_type)
_XFUNCPROTOEND
{

switch(info->data_format){
	case(FLOAT): strcpy(data_type,"FLOAT"); break;
	case(INTEGER): strcpy(data_type,"INTEGER"); break;
	case(COMPLEX): strcpy(data_type,"COMPLEX"); break;
	case(BYTE): strcpy(data_type,"BYTE"); break;
	case(UNKNOWN): strcpy(data_type,"UNKNOWN"); break;
}
return(SUCCESS);
}



/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_set_data_type

USAGE
ierr= sep3d_set_data_type(char *sep3dname,   char *data_type)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep_3d structure
data_type -  char*   file type for the dataset (FLOAT, INTEGER, COMPLEX, or BYTE)



RETURN VALUES
0   =  if it works correctly

DESCRIPTION
Sets the data type into the dataset
>
*/





#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_data_type(char *sep3dname,char *data_type)
_XFUNCPROTOEND
#else
int sep3d_set_data_type(sep3dname ,data_type)
char *sep3dname,*data_type;
#endif 
{
int i1,ierr;
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));
return(SEP3D_set_data_type(info,data_type));


}

_XFUNCPROTOBEGIN
int SEP3D_set_data_type(sep_3d *info,char *data_type)
_XFUNCPROTOEND
{
	if(0==strcmp("FLOAT",data_type)) info->data_format=FLOAT;
	else if(0==strcmp("INTEGER",data_type)) info->data_format=INTEGER;
	else if(0==strcmp("COMPLEX",data_type)) info->data_format=COMPLEX;
	else if(0==strcmp("BYTE",data_type)) info->data_format=BYTE;
	else  return(sepwarn(INVALID_DATA,"Unknown file type %s for sep3dtag=%s \n",
   data_type,info->name));
return(SUCCESS);
}
/*  $Id: type_usage.c,v 1.2 2004/04/08 22:32:28 bob Exp $ */
