#define SET_SDOC 1
/*
Author: Robert Clapp, ESMB 463, 7230253

Date Created:Wed Nov 7 15:18:13 PDT 2001

Purpose: 

*/	 

#include <seplib.h> 
#include<string.h>
#include <superset.h>
#include "superset_internal.h"

#ifndef YES
#define YES 1
#define NO 0
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int search_key_string(int i2,char *title_out,char *key_string);
int grab_head_par(char *name,double *value);
int quick_d_int_array(sep_3d *info,int,int *array,double *value);
int quick_d_float_array(sep_3d *info,int,float *array,double *value);
int quick_int_d_array(sep_3d *info,int *array,double *value);
int quick_float_d_array(sep_3d *info,float *array,double *value);

int sep3d_add_key(char *sep3dname,char *keyname, char *keytype);

_XFUNCPROTOEND
#else
int search_key_string();
int grab_head_par();
int quick_d_float_array();
int quick_d_int_array();
int quick_float_d_array();
int quick_int_d_array();

int sep3d_add_key();
#endif

int *ipt_math,*nd_math,ndim_math;
float *od_math,*dd_math;
sep_3d *info_math;
int ikey_math,nelem_math;
int *temp_math_head;


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_check_add_keys(char  *sep3dname)
_XFUNCPROTOEND
#else
int sep3d_check_add_keys(sep3dname)
char *sep3dname;
#endif
{
sep_3d *info;
char *temp_ch;
char big_buf[10000];
char key_name[1000];
char key_type[1000];
char temp2[1000];
char key_eqn[10000];
int nextra,i;
char *ptr;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */

if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

info->nkeys_in=info->nkeys;


if(info->usage==OUTPUT)return(SUCCESS);

temp_ch=(char*)malloc(sizeof(char)*(12+strlen(info->name)));
sprintf(temp_ch,"%s-extra_keys",info->name);
if(0==getch(temp_ch,"s",big_buf)){
free(temp_ch);
return(SUCCESS);
}
else { /*we are going to do additional mapping */
free(temp_ch);
nextra=1;
for( ptr=big_buf; *ptr!='\0'; ptr++ ){ if(*ptr == ':') nextra++;}

info->exp=(char **) alloc( nextra * sizeof(char*));

  for(i=0; i < nextra;i++){
    if(SUCCESS!=search_key_string(i,key_name,big_buf)){
      return(sepwarn(NOT_MET,"Internal error searching extra_keys string \n"));
    }
    sprintf(temp2,"%s-%s-type",info->name,key_name);
    if(0==getch(temp2,"s",key_type)) strcpy(key_type,"scalar_float");
    sprintf(temp2,"%s-%s-eqn",info->name,key_name);
    if(0==getch(temp2,"s",key_eqn)){
      return(sepwarn(NOT_MET,"Key equation not specified for extra key %s \n",
         key_name));
    } 

  if(0==strncmp(key_type,"scalar_float",12)||
       0==strncmp(key_type,"scalar_int",10)){
      if(SUCCESS!=sep3d_add_key(sep3dname,key_name,key_type)){
        return(sepwarn(NOT_MET,"trouble adding key %s \n",key_type));
      }
    }
    else return(sepwarn(NOT_MET,"invalid key type %s: acceptable types : scalar_float and scalar float \n",key_type));
 
      info->exp[i]=(char*) alloc((strlen(key_eqn)+1)*sizeof(char));
      strcpy(info->exp[i],key_eqn);
  }
 info->nextra_keys=nextra;
}
if(info->file_format==REGULAR) info->file_format=HEADER;


return(SUCCESS);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int calc_additional_headers(sep_3d  *info,int *nwind,int *fwind, int *jwind)
_XFUNCPROTOEND
#else
int calc_additional_headers(info,nwind,fwind,jwind)
sep_3d *info;
int *nwind,*fwind, *jwind;
#endif
{
int *idim,wind_elem,i1;
int i,ikey,i2,temp2i;
double *array;



  /* first make  a list of the axis elements we are reading in */
  idim=(int *) alloc(sizeof(int) * (info->ndims-1));
  nelem_math=1;
  for(i1=0; i1 < info->ndims-1; i1++) nelem_math=nelem_math*nwind[i1];
  ipt_math=(int *) alloc(sizeof(int) * nelem_math);
  for(i2=0; i2< nelem_math; i2++){
    h2c(i2,nwind,info->ndims-1,idim); /* convert to helical coordinate sys*/
    for(i1=0; i1< info->ndims-1; i1++)idim[i1]=fwind[i1]+jwind[i1]*idim[i1];
    c2h(&temp2i,&(info->n[1]),info->ndims-1,idim); /*convert back to reg*/
    ipt_math[i2]=temp2i;
  }
  free(idim);
info_math=info;
#ifdef SU_SUPPORT
if(info->su_input){
  nd_math=&info->n_su[1];
  od_math=&info->o_su[1];
  dd_math=&info->d_su[1];
  ndim_math=info->su_ndims-1;
}
else{
#endif
  nd_math=&info->n[1];
  od_math=&info->o[1];
  dd_math=&info->d[1];
  ndim_math=info->ndims-1;
#ifdef SU_SUPPORT
}
#endif



  temp_math_head=(int *) alloc(sizeof(int) * nelem_math);
  array=(double *) alloc(sizeof(double) * nelem_math);

  /* loop over the expressions */
  for(ikey_math=info->nkeys_in,ikey=0;ikey_math<info->nkeys;ikey++,ikey_math++)
{

    if(0!=evaluate_expression(info->exp[ikey],
       grab_head_par,nelem_math,array))
        return(sepwarn(NOT_MET,"trouble evaluating expression  %s for key=%s\n",
         info->exp[ikey],info->keyname[ikey_math]));


    if(0==strncmp(info->keytype[ikey_math],"scalar_float",12)){
      if(0!=quick_d_float_array(info,ikey_math,(float*)info->headers,array)){
        return(sepwarn(NOT_MET,"trouble converting headers for key=%s\n",
          info->keyname[ikey_math]));
       }
    }
    else if(0==strncmp(info->keytype[ikey_math],"scalar_int",10)){
      if(0!=quick_d_float_array(info,ikey_math,(float*)info->headers,array)){
        return(sepwarn(NOT_MET,"trouble converting headers for key=%s\n",
         info->keyname[ikey_math]));
       }
    }
    else return(sepwarn(NOT_MET,"unknown key type(%s) for key(%s) \n",
          info->keytype[ikey_math], info->keyname[ikey_math]));
  }
  free(ipt_math);
  free(array);
  free(temp_math_head);
  return(SUCCESS);

}
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int quick_d_float_array(sep_3d *info,int iloc,float *array,double *value)
_XFUNCPROTOEND
#else
int quick_d_float_array(info,iloc,array,value)
sep_3d *info;
int iloc;
float *array;
double *value;
#endif
{
int ih,ia;

for(ih=0,ia=iloc;ih < info->nh;ih++,ia=ia+info->nkeys) {
   array[ia]=(float)value[ih];
}
return(SUCCESS);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int quick_d_int_array(sep_3d *info,int iloc,int *array,double *value)
_XFUNCPROTOEND
#else
int quick_d_int_array(info,iloc,array,value)
sep_3d *info;
int iloc;
int *array;
double *value;
#endif
{
int ih,ia;

for(ih=0,ia=iloc;ih < info->nh;ih++,ia=ia+info->nkeys) array[ia]=(int)value[ih];
return(SUCCESS);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int quick_float_d_array(sep_3d *info,float *array,double *value)
_XFUNCPROTOEND
#else
int quick_float_d_array(info,array,value)
sep_3d *info;
float *array;
double *value;
#endif
{
int ih;

for(ih=0;ih < info->nh;ih++) value[ih]=(double)array[ih];
return(SUCCESS);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int quick_int_d_array(sep_3d *info,int *array,double *value)
_XFUNCPROTOEND
#else
int quick_int_d_array(info,array,value)
sep_3d *info;
int *array;
double *value;
#endif
{
int ih;

for(ih=0;ih < info->nh;ih++) value[ih]=(double)array[ih];
return(SUCCESS);
}





#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int grab_head_par(char *name,double *value)
_XFUNCPROTOEND
#else
int grab_head_par(char *name,double *value)
char *name;
double *value;
#endif
{
int nsz;
int an_axis;
int naxis;
int iaxis[9];
int logic,i;
nsz=strlen(name);
if(nsz==5){ /* could be an axis */
 if( name[0]=='a' && name[1]=='x' && name[2]=='i' && name[3]=='s' ){
   naxis=atoi(&name[4])-2;
   if(naxis<0 || naxis >= ndim_math)
    seperr("invalid request for axis %d. Only axis2 to  axis%d valid for this dataset \n",naxis+1,ndim_math);
   
   for(i=0; i <  nelem_math; i++){
    h2c(ipt_math[i],nd_math,ndim_math,iaxis); 
    value[i]=iaxis[naxis]*dd_math[naxis]+od_math[naxis];
   }
    
  return(0);
 }
}
for(i=0;i < ikey_math; i++){
  if(SUCCESS==strcmp(name,info_math->keyname[i])){
    if(0!=sep3d_grab_header_vals_i(info_math->name,i+1,temp_math_head)){
      seperr("trouble grabbing key values for key %s \n",name);
    }
    if(0==strncmp(info_math->keytype[i],"scalar_float",12)){
       if(0!=quick_float_d_array(info_math,(float*)temp_math_head,value))
        sepwarn(NOT_MET,"trouble converting from float to double \n");
    }
    else if(0==strncmp(info_math->keytype[i],"scalar_int",10)){
       if(0!=quick_int_d_array(info_math,(int*)temp_math_head,value))
        sepwarn(NOT_MET,"trouble converting from int to double \n");
    }
    else {
      seperr("unknown keytype %s \n",info_math->keytype);
    }
   return(0);
  }
}
seperr("can't interpret key conversion, trouble with %s \n",name);
    
return(0);
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int search_key_string(int i2,char *title_out,char *key_string)
_XFUNCPROTOEND
#else
int search_key_string(i2,char *key_string)
/*  find the i3rd membr of labels=first:second:third*/
int i2;
char *title_out;
#endif
{
  char *ptr;
  int i, colon,junk;
  colon = 0;
  title_out[0] = '\0';
  junk=1;
  for( ptr=key_string; *ptr!='\0'; ptr++ ) {
    if(*ptr == ':') {
      colon++;
      }
    else if( colon == i2 ) {
      for( i=0; *ptr!='\0' && *ptr !=':'; ptr++) {
        title_out[i++] = *ptr;
        }
      title_out[i] = '\0';
      junk=0;
      break;
      }
    }
  return(junk);
}

