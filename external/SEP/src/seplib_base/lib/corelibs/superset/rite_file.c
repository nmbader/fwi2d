#define SET_SDOC 1
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Sun Aug 16 15:04:02 PDT 1998

Purpose: 

*/	 

#include "superset_internal.h"
#include<string.h>
#ifndef YES
#define YES 1
#endif
#ifndef NO
#define NO 0
#endif
#include <superset.h> 


/*<
sep3d_rite_format

USAGE
ierr= sep3d_rite_format(char *tag, char *sep3dname)

INPUT PARAMETERS
tag       -  char*   tag  to write out to
sep3dname -  char*   tag associated with sep_3d structure


RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Writes out structure to disk
>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_rite_format(char *tag,char *sep3dname)
_XFUNCPROTOEND
#else
int sep3d_rite_format(tag,sep3dname)
char *tag,*sep3dname; 
#endif 
{

return(sep3d_rite_format_ith(tag,sep3dname,0));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_rite_format_ith(char *tag,char *sep3dname, int ith)
_XFUNCPROTOEND
#else
int sep3d_rite_format_ith(tag,sep3dname,int ith)
char *tag,*sep3dname; 
int ith;
#endif 
{


sep_3d *info;
int i,tot,ierr,in_order;
float one;
char temp_ch[1024],temp_name[2049],temp2_ch[1024];


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));





/*if(info->usage == INPUT) return(sepwarn(NOT_MET,*/
/*		"can not write out to tag in (%s) \n",sep3dname));*/
if(info->file_format == UNSPECIFIED) return(sepwarn(NOT_MET,
    "can not write to a file whose file_format is unspecified1 (%s) \n",
	sep3dname));

i=0;
auxputch("junkME","d",&i,tag);

if(info->data_format == UNKNOWN)  return(sepwarn(NOT_MET,
    "can not write to a file whose data_format is unspecified2 (%s) \n",
     sep3dname));
else if(info->data_format==INTEGER) set_format(tag,"xdr_int");
else if (info->data_format==BYTE)   set_format(tag,"native_byte");
if(info->ndims==0) return(sepwarn(NOT_MET,
    "can not write to a file whose dimensions are  unspecified3 (%s) \n",
  sep3dname));

one=1.;


if(info->file_format==GRID){
	if(SUCCESS!=fget_grid_format_tag(tag,temp_ch))
		return(sepwarn(FAIL_OTHER,
		"rite_format:trouble obtaining grid format tag for %s \n",tag));

	for(i=2; i<=info->ndims;i++){
		if(SUCCESS!=sep_put_grid_axis_par(tag,&i,&(info->n[i-1]),&(info->o[i-1]),
    &(info->d[i-1]), info->label[i-1])) return(sepwarn(FAIL_OTHER,
    "trouble writing out grid par for tag %s \n",tag));
		ierr=sep3d_rite_axis(temp_ch,i,info->n[i-1],info->o[i-1],info->d[i-1],
     info->label[i-1],info->unit[i-1]);
		if(ierr!=SUCCESS) return(ierr);
	}

}
if(info->file_format!=REGULAR){
	if(info->nkeys==0) return(sepwarn(NOT_MET,
    "number of keys has not been set for tag %s \n",tag));

	if(info->file_format==HEADER){
		if(SUCCESS!=fget_header_format_tag(tag,temp_ch))
			return(sepwarn(FAIL_OTHER,
    	"trouble obtaining header format tag for tag %s \n",tag));

		for(i=2; i<=info->ndims;i++){
			ierr=sep3d_rite_axis(temp_ch,i,info->n[i-1],info->o[i-1],info->d[i-1],
     info->label[i-1],info->unit[i-1]);
		 if(ierr!=0) return(ierr);
		}
	}

  if(sep3d_grab_inorder(info->name,&in_order)!=0)
    return(sepwarn(NOT_MET,"trouble getting in_order \n"));



  if(in_order!=1) i=info->nkeys+1;
	else i=info->nkeys;

	if(SUCCESS!=sep_put_number_keys(tag,&i))
			return(sepwarn(FAIL_OTHER,
    	"trouble putting number of keys for for tag %s \n",tag));

	for(i=1; i<= info->nkeys; i++){
		if(SUCCESS!=sep_put_key(tag,info->keyname[i-1],info->keytype[i-1],
     info->keyfmt[i-1],&i)) return(sepwarn(FAIL_OTHER,
      "trouble putting  key %d for for tag %s \n",i,tag));
	}

  if(in_order!=1){
		if(SUCCESS!=sep_put_key(tag,"data_record_number","scalar_int","xdr_int",&i))
			return(sepwarn(FAIL_OTHER, 
      	"trouble writing drn  key for for tag %s \n",tag));
		i=0;
  }
  else i=1;
	if(SUCCESS!=auxputch("same_record_number","d",&i,tag))
			return(sepwarn(FAIL_OTHER, 
      	"trouble writing same_record_number for for tag %s \n",tag));

	if(SUCCESS!=sep3d_rite_ntraces(tag,sep3dname))
   return(sepwarn(FAIL_OTHER,
        "trouble writing out ntraces for tag %s \n",tag));

   
   

}
else{ /* REGULAR SEP3D FILE */
/*if(0!=strcmp("/scr1/bob/stemp.H",tag)){}*/
if(1==1){
		for(i=2; i<=info->ndims;i++){
			ierr=sep3d_rite_axis(tag,i,info->n[i-1],info->o[i-1],info->d[i-1],
     info->label[i-1],info->unit[i-1]);
			if(SUCCESS!=ierr) return(ierr);
		}
	if(SUCCESS!=sep_set_no_grid(tag)) return(sepwarn(FAIL_OTHER,
    "trouble setting no grid for tag %s \n",tag));
	if(SUCCESS!=sep_set_no_headers(tag))return(sepwarn(FAIL_OTHER, 
    "trouble setting no headers for tag %s \n",tag));

   strcpy(temp_ch,"hff=-1 gff=-1 ");
   for(i=info->ndims+1;i < 10; i++) sprintf(temp2_ch,"%s n%d=1 ",temp_ch,i);
   auxputhead(tag,"%s",temp2_ch);

}
}

/*if(0!=strcmp("/scr1/bob/stemp.H",tag)){}*/
i=1;

ierr=sep3d_rite_axis(tag,1,info->n[0],info->o[0],info->d[0],
     info->label[0],info->unit[0]);
if(ierr!=SUCCESS) return(ierr);

if(info->file_format==REGULAR) auxputch("hff","s","-1",tag);
if(info->file_format!=GRID) auxputch("gff","s","-1",tag);


switch(info->data_format){
	case   (COMPLEX): i=8; break;
	case   (BYTE): i=1; break;
	default: i=4; break;
}
if(SUCCESS!=auxputch("esize","d",&i,tag)) return(sepwarn(FAIL_OTHER,
	"trouble writing out esize to tag %s \n",tag));




return(SUCCESS);
}


/*<
sep3d_rite_ntraces

USAGE
ierr= sep3d_rite_ntraces(char *tag, char *sep3dname)

INPUT PARAMETERS
tag       -  char*   tag  to write out to
sep3dname -  char*   tag associated with sep_3d structure


RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Writes out number of traces to disk
>
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_set_rite_status(char *sep3dname,int data,int  header)
_XFUNCPROTOEND
#else
int sep3d_set_rite_status(sep3dname,data,header)
char *sep3dname;
int data, header;
#endif
{
sep_3d *info;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));
/*if(info->usage == INPUT) return(sepwarn(NOT_MET,*/
/*		"can not write out to tag in (%s) \n",sep3dname));*/
if(info->file_format == UNSPECIFIED) return(sepwarn(NOT_MET,
    "can not write to a file whose file_format is unspecified (%s) \n",
	sep3dname));

info->wrote_data=data;
info->wrote_headers=header;
return(SUCCESS);

}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_rite_ntraces(char *tag,char *sep3dname)
_XFUNCPROTOEND
#else
int sep3d_rite_ntraces(tag,sep3dname)
char *tag,*sep3dname; 
#endif 
{
sep_3d *info;
int i,ik,nt;
float one=1.;
char temp_name[2049];


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

/*if(info->usage == INPUT) return(sepwarn(NOT_MET,*/
/*		"can not write out to tag in (%s) \n",sep3dname));*/
if(info->file_format == UNSPECIFIED) return(sepwarn(NOT_MET,
    "can not write to a file whose file_format is unspecified (%s) \n",
	sep3dname));

			

if(info->wrote_data==YES && info->file_format!=REGULAR){
      nt=info->ntraces;
	i=2;if(SUCCESS!=sep_put_data_axis_par(tag,&i,&nt,&one,&one,
 	 "trace number")) return(sepwarn(FAIL_OTHER,
       "trouble writing out data axis for tag %s \n",tag));
}


if(info->file_format==GRID || (info->file_format==HEADER && info->ndims<3 && info->wrote_headers==YES)){
      nt=info->ntraces;
	i=2;if(SUCCESS!=sep_put_header_axis_par(tag,&i,&nt,&one,&one,
   "trace number"))  return(sepwarn(FAIL_OTHER,
    "trouble writing out header axis par for tag %s \n",tag));
	auxputhead(tag,"n3=1 n4=1 n5=1 n6=1 n7=1 n8=1 n9=1 \n");
}
return(SUCCESS);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_rite_file_stat(char *sep3dname,int data_stat,int header_stat)
_XFUNCPROTOEND
#else
int sep3d_rite_file_stat(sep3dname,data_stat,header_stat)
char *sep3dname; 
int data_stat,header_stat;
#endif 
{
sep_3d *info;
info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */

if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

info->wrote_data=data_stat;
info->wrote_headers=header_stat;

return(SUCCESS);
}

/*  $Id: rite_file.c,v 1.4 2004/07/08 18:15:32 bob Exp $ */
