#define SET_SDOC 1
/*
Author: Robert Clapp, ESMB 463, 7230253

Date Created:Sun Aug 16 13:18:13 PDT 1998

Purpose: 

*/	 

#include <seplib.h> 
#include <superset.h>
#include "superset_internal.h"
#include "sepstream.h"
#include "sep_main_external.h"
#include<string.h>

#ifndef YES
#define YES 1
#define NO 0
#endif

extern int have_header_format_tag(char *tag, char *tag_temp);
extern int have_grid_format_tag(char *tag, char *tag_temp);
extern int sep3d_tag_init_thread(char  *tag,char *sep3dout,char *usage,int ithread);
extern int sep_thread_num(void);
extern int sep_num_thread(void);
extern int sep3d_set_inorder(const char *tag);
extern int sep3d_unset_inorder(const char *tag);
extern int tag_exists(const char *tag);

/*<
sep3d_struct_init

Usage
ierr=sep3d_struct_init(char  *sep3din, char *stuct_out,char *usage)

Return Values
SUCESS        =  if works correctly
EXISTS        =  if tag has already been declared as something else
INVALID_USAGE =  if the usage passed is unacceptable
FAIL_OTHER    =  if fails for other reasons

Input Paramters
sep3din     -  char*   pointer to structure to copy from
usage       -  char*   usage for new structure

Output Paramters
sep3dout    -  char*   pointer to structure to copy to

Description
Initialize sep3d tag  from another tag
*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_struct_init(char  *sep3din, char *sep3dout,char *usage) 
_XFUNCPROTOEND
#else
int sep3d_struct_init(sep3din,sep3dout,usage) 
char *sep3din,*sep3dout,*usage; 
#endif 
{
int ndim,icreate;
int  i1,master_data;
sep_3d *input,*output;
char temp_ch[1024];

input = tag_info_sep3d(sep3din,INQUIRE);
if(input==SEPNULL) {
	return(sepwarn(FAIL_OTHER,"Invalid SEP3d structure to copy from (%s) \n",sep3din));
}

/*FIRST INITIALIZE THE TAG */
if(0==strcmp(usage,"INPUT")) {
	output=tag_info_sep3d(sep3dout,INPUT);
	if(output->usage!=INPUT){
		return(sepwarn(EXISTS,"Tag already declared as a diffent type\n"));
	}
/*	if(0!=strcmp("in",sep3dout) && NULL== auxin(sep3dout))*/
/*		return(sepwarn(INVALID_DATA,"trouble opening %s as type scratch \n"));*/
}
else if(0==strcmp(usage,"OUTPUT")) {
	output=tag_info_sep3d(sep3dout,OUTPUT);
	if(output->usage!=OUTPUT){
		return(sepwarn(EXISTS,"Tag already declared as a diffent type\n"));
	}
/*	if(0!=strcmp("out",sep3dout) && NULL== auxout(sep3dout))*/
/*		return(sepwarn(INVALID_DATA,"trouble opening %s as type scratch \n"));*/
}
else if(0==strcmp(usage,"SCRATCH")) {
	output=tag_info_sep3d(sep3dout,SCRATCH);
	if(output->usage!=SCRATCH) {
		return(sepwarn(EXISTS,"Tag already declared as a diffent type\n"));
	}
/*	if(NULL== auxinout(sep3dout))*/
/*		return(sepwarn(INVALID_DATA,"trouble opening %s as type scratch \n"));*/
}
else{
	return(sepwarn(INVALID_USAGE,
     "Invalid usage (INPUT, OUTPUT, or SCRATCH) %s \n",usage));
}

output = tag_info_sep3d(sep3dout,INQUIRE);
if(output==SEPNULL) {
	return(sepwarn(FAIL_OTHER,"Unable to create SEP3d structure \n"));
}
	


/*Next copy the miscellaneous information */
output->data_format=input->data_format;
output->file_format=input->file_format;
output->ntraces=input->ntraces;
output->drn=input->drn;


/*Next copy the axis information */
if(input->ndims>0){
	if(SUCCESS!=sep3d_set_ndims(sep3dout,input->ndims)){
		return(sepwarn(FAIL_OTHER,
     "trouble setting the number of dimensions in the output \n"));
	}
	for(i1=0; i1< input->ndims; i1++){
		if(input->n[i1]!=0){
			if(SUCCESS!=sep3d_set_axis(sep3dout,i1+1, input->n[i1],
				input->o[i1],input->d[i1],input->label[i1],input->unit[i1])){
					return(sepwarn(FAIL_OTHER,"trouble copying axis %d \n",i1+1));
			}
		}
	}
}
			
/*Next copy the key information */
if(input->nkeys>0){
	if(SUCCESS!=sep3d_set_nkeys(sep3dout,input->nkeys)){
		return(sepwarn(FAIL_OTHER,
     "trouble setting the number of dimensions in the output \n"));
	}
	for(i1=0; i1< input->nkeys; i1++){
			if(SUCCESS!=sep3d_set_key(sep3dout,i1+1, input->keyname[i1],
				input->keytype[i1],input->keyfmt[i1])){
					return(sepwarn(FAIL_OTHER,"trouble copying key %d \n",i1+1));
			}
	}
}

/*next add any user specified math keys */
if(0!=sep3d_check_add_keys(sep3dout)){
					return(sepwarn(FAIL_OTHER,"trouble setting up additional keys \n"));
}
output->nkeys_in=input->nkeys_in;

/*  mkrandom_string(sep3dout,temp_ch);*/
master_data=1;



if(input->nextra_keys >0){
  output->nextra_keys=input->nextra_keys;
  output->nkeys_in=input->nkeys_in;
  output->exp=(char**)malloc(sizeof(char*)*output->nextra_keys);
  for(i1=0; i1 < output->nextra_keys; i1++){
    output->exp[i1]=(char*)malloc(sizeof(char)*(strlen(input->exp[i1])+1));
    strcpy(output->exp[i1],input->exp[i1]);
  }
}
/*  sep3d_grab_file_type(sep3din,temp_ch);*/
/*  sep3d_grab_file_type(sep3dout,temp_ch);*/

 
return (SUCCESS) ;
}

/*<
sep3d_par_init

Usage
ierr=sep3d_par_init(char  *sep3din, char *usage)

Return Values
SUCESS        =  if works correctly
EXISTS        =  if tag has already been declared as something else
INVALID_USAGE =  if the usage passed is unacceptable
FAIL_OTHER    =  if fails for other reasons

Input Paramters
usage       -  char*   usage for new structure

Output Paramters
sep3dout    -  char*   pointer to structure to create

Description
Initialize sep3d tag

*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_par_init(char  *sep3dout,char *usage) 
_XFUNCPROTOEND
#else
int sep3d_par_init(sep3dout,usage) 
char *sep3dout,*usage; 
#endif 
{ 
int ndim;
int  i1;
sep_3d *input,*output;
char temp_ch[1024];


/*FIRST INITIALIZE THE TAG */
if(0==strcmp(usage,"INPUT")) {
	output=tag_info_sep3d(sep3dout,INPUT);
	if(output->usage!=INPUT){
		return(sepwarn(EXISTS,"Tag %s already declared as a diffent type\n",sep3dout));
	}
}
else if(0==strcmp(usage,"OUTPUT")) {
	output=tag_info_sep3d(sep3dout,OUTPUT);
	if(output->usage!=OUTPUT){
		return(sepwarn(EXISTS,"Tag already declared as a diffent type\n"));
	}
}

else if(0==strcmp(usage,"SCRATCH")) {
	output=tag_info_sep3d(sep3dout,SCRATCH);
	if(output->usage!=SCRATCH){
		return(sepwarn(EXISTS,"Tag already declared as a diffent type\n"));
	}
}
else{
	return(sepwarn(INVALID_DATA,
	"Invalid usage (INPUT, OUTPUT, or SCRATCH) %s \n",usage));
}
                                                                                 
	return(SUCCESS);
}

/*<
sep3d_tag_init

Usage
ierr=sep3d_tag_init(char  *tag, char *sep3dout, char *usage)

Return Values
SUCESS        =  if works correctly
EXISTS        =  if tag has already been declared as something else
INVALID_USAGE =  if the usage passed is unacceptable
FAIL_OTHER    =  if fails for other reasons

Input Paramters
tag         -  char*   tag to initialize the tag from
usage       -  char*   usage for new structure

Output Paramters
sep3dout    -  char*   pointer to structure to create

Description
Initialize sep3d tag from tag
*/




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_tag_init(char  *tag,char *sep3dout,char *usage) 
_XFUNCPROTOEND
#else
int sep3d_tag_init(tag,sep3dout,usage) 
char *tag,*sep3dout,*usage; 
#endif 
{


return(sep3d_tag_init_thread(tag,sep3dout,usage,-1));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
#endif 
int sep3d_tag_init_thread(char  *tag,char *sep3dout,char *usage,int ithread) 
_XFUNCPROTOEND
{

int ndim,master_data;
int  i1, ntr;
sep_3d *input,*output;
char temp_ch1[1024],temp_ch2[1024],temp_ch3[1024],tag_temp[1024],temp_ch[1204];
int n,tempi,temp2i,esize,n1,drn,nk;
int ignore_gff,ignore_hff;
int inorder;
int nbuf;
char *buf;
float o,d;
int has_headers;
long long bigi;



if(ithread==sep_thread_num() || ithread==-1){


if(0==strcmp(usage,"SCRATCH")){
	if(NULL== auxinout(sep3dout))
		return(sepwarn(INVALID_DATA,"trouble opening %s as type scratch \n"));
}
if(0==auxpar("esize","d",&esize,tag)) esize=4;


/*FIRST INITIALIZE THE TAG */
if(0==strcmp(usage,"INPUT")) {
	output=tag_info_sep3d(sep3dout,INPUT);
	if(output->usage!=INPUT){
		return(sepwarn(EXISTS,"Tag already declared as a diffent type\n"));
	}
}
else if(0==strcmp(usage,"OUTPUT")) {
	output=tag_info_sep3d(sep3dout,OUTPUT);
	if(output->usage!=OUTPUT){
		return(sepwarn(EXISTS,"Tag already declared as a diffent type\n"));
	}
}
else if(0==strcmp(usage,"SCRATCH")) {
	output=tag_info_sep3d(sep3dout,SCRATCH);
	if(output->usage!=SCRATCH){
		return(sepwarn(EXISTS,"Tag already declared as a diffent type\n"));
	}
}
else{
	return(sepwarn(INVALID_DATA,
   "Invalid usage (INPUT, OUTPUT, or SCRATCH):%s \n",usage));
}

sprintf(temp_ch1,"%s.%s",tag,"ignore_hff");
if(0==getch(temp_ch1,"d",&ignore_hff)) {
  sprintf(temp_ch1,"%s-%s",tag,"ignore_hff");
  if(0==getch(temp_ch1,"d",&ignore_hff)) ignore_hff=0;
}
sprintf(temp_ch1,"%s.%s",tag,"ignore_gff");
if(0==getch(temp_ch1,"d",&ignore_gff)) {
  sprintf(temp_ch1,"%s-%s",tag,"ignore_gff");
  if(0==getch(temp_ch1,"d",&ignore_gff)) ignore_gff=0;;
}

/* NOW DETERMINE THE TYPE OF DATASET */
has_headers=NO;
if(0== ignore_gff &&
   SUCCESS==have_grid_format_tag(tag,tag_temp)){ /*DATA HAS A GRID */
	has_headers=YES;
	if(SUCCESS!=sep_get_non_default_axes(tag_temp,&ndim,2)){
		return(sepwarn(FAIL_OTHER,"trouble reading number of grid axes\n"));
	}
		if(SUCCESS!=sep3d_set_ndims(sep3dout, ndim)){
			return(sepwarn(FAIL_OTHER,"trouble setting number of dimensions \n"));
		}
		if(SUCCESS!=sep3d_set_file_type(sep3dout, "GRID")){
			return(sepwarn(FAIL_OTHER,"trouble setting file_type \n"));
		}
	for(i1=2; i1<= ndim; i1++){
		if(SUCCESS!=sep_get_grid_axis_par(tag,&i1,&n,&o,&d,temp_ch1)){
			return(sepwarn(FAIL_OTHER,"trouble getting grid axis %d \n",i1));
		}
		sprintf(temp_ch2,"%s%d","unit",i1);
		temp2i=auxpar(temp_ch2,"s",temp_ch3,tag_temp);	
		if(SUCCESS!=sep3d_set_axis(sep3dout,i1,n,o,d,temp_ch1,temp_ch2)){
			return(sepwarn(FAIL_OTHER,"trouble setting axis %d \n",i1));
		}
	}
	if(SUCCESS!=sep_get_number_header_axes(tag,&tempi)){
		return(sepwarn(FAIL_OTHER,"trouble getting number of header axes \n"));
	}
	for(i1=2,bigi=1; i1<= tempi; i1++){
		if(SUCCESS!=sep_get_header_axis_par(tag,&i1,&n,&o,&d,temp_ch1)){
			return(sepwarn(FAIL_OTHER,"trouble getting header axis %d \n",i1));
		}
		bigi=bigi*(long long)n;
	}
	if(SUCCESS!=SEP3D_set_ntraces(output,bigi)){
		return(sepwarn(FAIL_OTHER,"trouble setting number of traces \n"));
	}
}
else if(ignore_hff==0 &&
   SUCCESS==have_header_format_tag(tag,tag_temp)){ /*DATA HAS A HEADER */
	has_headers=YES;
	if(SUCCESS!=sep_get_non_default_axes(tag_temp,&ndim,1))
		return(sepwarn(FAIL_OTHER,
     "trouble reading number of header axes for tag  %s \n",tag));
		if(SUCCESS!=sep3d_set_ndims(sep3dout, ndim))
			return(sepwarn(FAIL_OTHER,"trouble setting number of dimensions (%s) \n",
      sep3dout));
		if(SUCCESS!=sep3d_set_file_type(sep3dout, "HEADER"))
			return(sepwarn(FAIL_OTHER,"trouble setting file_type %s \n",sep3dout));
	for(i1=2,bigi=0; i1<= ndim; i1++){
		if(SUCCESS!=sep3d_reed_axis(tag_temp,i1,&n,&o,&d,temp_ch1,temp_ch2))
			return(sepwarn(FAIL_OTHER,"trouble getting header axis %d (%s)\n",i1,tag));

		if(SUCCESS!=sep3d_set_axis(sep3dout,i1,n,o,d,temp_ch1,temp_ch2))
			return(sepwarn(FAIL_OTHER,"trouble setting axis %d  (%s)\n",i1,sep3dout));
		bigi=bigi*(long long)n;
	}
	if(SUCCESS!=SEP3D_set_ntraces(output,bigi))
		return(sepwarn(FAIL_OTHER,"trouble setting number of traces (%s) \n",
     sep3dout));
}
else if(0==sep_get_non_default_axes(tag,&ndim,1)){ /*DATA IS REGULAR */
		if(SUCCESS!=sep3d_set_ndims(sep3dout, ndim))
			return(sepwarn(FAIL_OTHER,"trouble setting number of dimensions (%s) \n",
        sep3dout));
		if(SUCCESS!=sep3d_set_file_type(sep3dout, "REGULAR"))
			return(sepwarn(FAIL_OTHER,"trouble setting file_type (%s) \n",sep3dout));
	ntr=1;
	for(i1=2,bigi=1; i1<= ndim; i1++){
		if(SUCCESS!=sep3d_reed_axis(tag,i1,&n,&o,&d,temp_ch1,temp_ch2))
			return(sepwarn(FAIL_OTHER,
       "trouble getting header axis %d (%s)\n",i1,tag));
		 ntr=ntr*n;
		if(SUCCESS!=sep3d_set_axis(sep3dout,i1,n,o,d,temp_ch1,temp_ch2))
			return(sepwarn(FAIL_OTHER,"trouble setting axis %d (%s) \n",i1,sep3dout));
		bigi=bigi*(long long)n;
	}
	if(SUCCESS!=SEP3D_set_ntraces(output,bigi))
		return(sepwarn(FAIL_OTHER,"trouble setting number of traces (%s) \n",
     sep3dout));

	has_headers=NO;

	if(SUCCESS!=SEP3D_set_ntraces(output,bigi))
		return(sepwarn(FAIL_OTHER,"trouble writing number of traces (%s) \n",
     sep3dout));

	
}
else{
		 return(sepwarn(INVALID_STRUC,"invalid sep dataset (%s) \n",tag));
}



/*DEAL WITH THE FIRST AXIS IN THE DATASET */
tempi=1;
if(SUCCESS!=sep3d_reed_axis(tag,tempi,&n,&o,&d,temp_ch1,temp_ch2))
	return(sepwarn(FAIL_OTHER,"trouble obtaining data axis 1 (%s) \n",tag));
if(SUCCESS!=sep3d_set_axis(sep3dout,1,n,o,d,temp_ch1,temp_ch2))
	return(sepwarn(FAIL_OTHER,"trouble setting axis 1 (%s) \n",sep3dout));

	
/*SET THE DATA TYPE */
switch(esize){
	case(1):
		if(SUCCESS!=sep3d_set_data_type(sep3dout,"BYTE"))
			return(sepwarn(FAIL_OTHER,"trouble setting datatype (%s) \n",sep3dout));
		break;
	case(8):
		if(SUCCESS!=sep3d_set_data_type(sep3dout,"COMPLEX"))
			return(sepwarn(FAIL_OTHER,"trouble setting datatype (%s) \n",sep3dout));
		break;
	case(4):
		if(0==auxpar("data_format","s",temp_ch1,tag)) strcpy(temp_ch1,"xdr_float");
		if(0==strcmp(temp_ch1,"xdr_float") || 0==strcmp(temp_ch1,"native_float")) {
			if(SUCCESS!=sep3d_set_data_type(sep3dout,"FLOAT"))
				return(sepwarn(FAIL_OTHER,"trouble setting datatype (%s) \n",sep3dout));
		}
		else if(0==strcmp(temp_ch1,"xdr_int") || 0==strcmp(temp_ch1,"native_int")) {
			if(SUCCESS!=sep3d_set_data_type(sep3dout,"INTEGER"))
				return(sepwarn(FAIL_OTHER,"trouble setting datatype (%s) \n",sep3dout));
		}
		else return(sepwarn(INVALID_DATA,"unrecognized data type: %s (%s) \n",
     temp_ch1,sep3dout));
		break;
	default:
		return(sepwarn(INVALID_DATA,"Unacceptable esize:%d (%s)\n",tempi,sep3dout));
}
	

if(has_headers==YES){

	if(SUCCESS!=sep_get_number_keys(tag,&nk))
		return(sepwarn(FAIL_OTHER,"trouble getting number of keys (%s)  \n",tag));


	if(SUCCESS==sep_get_key_index(tag,"data_record_number",&drn)) tempi=nk-1;
	else tempi=nk;
	
	
	if(SUCCESS!=sep3d_set_nkeys(sep3dout,tempi))
		return(sepwarn(FAIL_OTHER,
     "trouble setting number of keys for  (%s)\n",sep3dout));
	
	n1=0;
	for(i1=1; i1 <= nk; i1++){
		if(i1!=drn){
			n1++;
			if(SUCCESS!=sep_get_key_name(tag,&i1,temp_ch1) ||
			SUCCESS!=sep_get_key_type(tag,&i1,temp_ch2) ||
			SUCCESS!=sep_get_key_fmt(tag,&i1,temp_ch3) )
			return(sepwarn(FAIL_OTHER,
       "trouble getting key name, type, or format for key %d  (%s) \n",i1,tag));
			if(SUCCESS!=sep3d_set_key(sep3dout,n1,temp_ch1,temp_ch2,temp_ch3))
				return(sepwarn(FAIL_OTHER,"trouble setting key %d  (%s) \n",i1,tag));
		}
	}
        if(SUCCESS==sep3d_grab_inorder(tag,&inorder)){;
          if(inorder==-1){/*We haven't determined if the dataset is inorder yet*/
          auxpar("same_record_number","s",&inorder,tag);
          if(inorder==1 || drn <1) sep3d_set_inorder(tag);
          else{
             inorder=0;
           sep3d_unset_inorder(tag);
            } 
           }
        }
	if(SUCCESS!=sep_get_number_header_axes(tag,&ndim))
		return(sepwarn(FAIL_OTHER,"trouble reading number of header axes (%s)\n",
     tag));
	for(i1=2,bigi=1; i1 <= ndim; i1++){
		if(SUCCESS!=sep_get_header_axis_par(tag,&i1,&n,&o,&d,temp_ch1))
			return(sepwarn(FAIL_OTHER,"trouble obtaining header axis par (%s) \n",
       tag));
		bigi=bigi*(long long)n;
	}
	if(SUCCESS!=SEP3D_set_ntraces(output,bigi))
		return(sepwarn(FAIL_OTHER,"trouble writing number of traces (%s) \n",tag));
}
if(0!=sep3d_check_add_keys(sep3dout)){
					return(sepwarn(FAIL_OTHER,"trouble setting up additional keys \n"));
}
 

                                                                                
}

if(sep_num_thread()>1 && ithread>=0){
if(0!=sep3d_broadcast_headers(tag, ithread))
  return(sepwarn(NOT_MET,"trouble broadcasting headers \n"));
}
master_data=1;
getch("master_data","d",&master_data);
sprintf(temp_ch,"%s.master_data",sep3dout);
getch(temp_ch,"d",&master_data);


return(SUCCESS);
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_close_tag(char *tag,char  *sep3dname) 
_XFUNCPROTOEND
#else
int sep3d_close_tag(tag,sep3dname) 
char *sep3dname,*tag;
#endif 
{ 
sep_3d *info;
int i;
char temp_ch[1024];
info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag(%s) invalid struc\n",sep3dname));


if(info->file_format==GRID){ 
	if(SUCCESS==fget_grid_format_tag(tag,temp_ch)) auxclose(temp_ch);
	if(SUCCESS==fget_header_format_tag(tag,temp_ch)) auxclose(temp_ch);
}
if(info->file_format==HEADER){
	if(SUCCESS==fget_header_format_tag(tag,temp_ch)) auxclose(temp_ch);
}

auxclose(tag);

return(SUCCESS);
}
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_delete(char  *sep3dname) 
_XFUNCPROTOEND
#else
int sep3d_delete(sep3dname) 
char *sep3dname;
#endif 
{
sep_3d *info;
int i;
char temp_ch[1024];
info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag(%s) invalid struc\n",sep3dname));

SEP3D_clean(info);

free(info->name);
info->name=SEPNULL;

SEP3D_del(info);



return(0);


}

_XFUNCPROTOBEGIN
int SEP3D_clean(sep_3d  *info) 
_XFUNCPROTOEND
{
int i;


if(0!=SEP3D_del_axes(info))
  return(sepwarn(NOT_MET,"trouble deieting axes info  tag=%s \n",info->name));

if(info->nkeys>0){
	for(i=0; i < info->nkeys; i++){
		if(SEPNULL!=info->keyname[i]) free(info->keyname[i]);
		if(SEPNULL!=info->keytype[i]) free(info->keytype[i]);
		if(SEPNULL!=info->keyfmt[i]) free(info->keyfmt[i]);
	}
	free(info->keyname); free(info->keyfmt); free(info->keytype);
}

if(info->nextra_keys>0){
	for(i=0; i < info->nextra_keys; i++){
		if(SEPNULL!=info->exp[i]) free(info->exp[i]);
}

  free(info->name);
  free(info->exp);
}

info->nh=0;
if(info->drn!=SEPNULL) free(info->drn); 
info->drn=SEPNULL;
if(info->headers!=SEPNULL) free(info->headers);
info->headers=SEPNULL;
SEP3D_delete_coord(info);



return(SUCCESS);
}


int sep_get_non_default_axes(char *tag,int *naxes,int istart){

  char ni[8];
    int  exist=1, i_axis=0, j_axis=0;
    int n;
    float o,d;
    char extra[26];


    i_axis=istart-1;
    if(istart==2) strcpy(extra,"_grid");
    else strcpy(extra,"");
    sprintf(ni, "o%d%s", i_axis+1,extra);
    exist=auxpar(ni, "f", &o, tag);
    sprintf(ni, "d%d%s", i_axis+1,extra);
    exist=auxpar(ni, "f", &d, tag);
    sprintf(ni, "n%d%s", i_axis+1,extra);
    exist=auxpar(ni, "d", &d, tag);
    if(1!=auxpar(ni, "d", &n, tag))
      seperr("invalid seplib3d dataset,%s, n1 must exist (%d) and by a single integer thread=%d \n",tag,exist,sep_thread_num());


     o=0.;d=1.;
    /* Auxpar i_axis from History File  */
    while (exist==1){
        if(n >1 || o !=0. || d!=1.) j_axis=i_axis;
        o=0.;d=1.;
        i_axis++;
        sprintf(ni, "o%d%s", i_axis,extra);
        exist=auxpar(ni, "f", &o, tag);
        sprintf(ni, "d%d%s", i_axis,extra);
        exist=auxpar(ni, "f", &d, tag);
        sprintf(ni, "n%d%s", i_axis,extra);
        exist=auxpar(ni, "d", &n, tag);
    }
    if(j_axis==0) j_axis=1;
    *naxes=j_axis;



    return 0;


}
int have_header_format_tag(char *tag, char *tag_temp){
int i;
i=auxpar("hff","s",tag_temp,tag);
if(i==1 && 0!=strcmp(tag_temp,"-1")){
   fget_header_format_tag(tag,tag_temp);
   return(0);
}
else return(1);
}


int have_grid_format_tag(char *tag, char *tag_temp){
int i;
i=auxpar("gff","s",tag_temp,tag);
if(i==1 && 0!=strcmp(tag_temp,"-1")){
   fget_grid_format_tag(tag,tag_temp);
    return(0);
}
else return(1);
}

int SEP3D_get_esize(sep_3d *info){
int esize;

switch(info->data_format){ /* set the appropriate esize */
  case(BYTE): esize=1; break;
  case(COMPLEX): esize=8; break;
  default:  esize=4; break;
}
return(esize);

}

int sep3d_copy_struct(char *sep3din, char *sep3dout){
sep_3d *input,*output;
int i,off,nsize;
char temp_ch[256];
char *buf;

	input=tag_info_sep3d(sep3din,OUTPUT);
	if(input==SEPNULL)
		return(sepwarn(NOT_MET,"Tag not allocated  %s\n",sep3din));
	output=tag_info_sep3d(sep3dout,INQUIRE);
	if(output==SEPNULL)
		return(sepwarn(NOT_MET,"Tag not allocated  %s\n",sep3dout));

  if(input->ndims != output->ndims)
		return(sepwarn(NOT_MET,"Copy call ndim not the same  %s-%s\n",sep3din,sep3dout));


  if(input->file_format!=output->file_format)
     return(sepwarn(NOT_MET,"Can't copy  because file formats different \n"));


  if(input->file_format!=REGULAR){
    if(0!=sep3d_header_copy(input,output))
       return(sepwarn(NOT_MET,"Trouble copying headers from %s to %s \n",sep3din,sep3dout));
    if(input->file_format==GRID){
/*      if(0!=sep3d_grid_copy(input,output))*/
/*       return(sepwarn(NOT_MET,"Trouble copying grid from %s to %s \n",sep3din,sep3dout));*/
     }
  }

  if(0!=SEP3D_coord_copy(input,output))
     return(sepwarn(NOT_MET,"Trouble copying coords from %s to %s \n",sep3din,sep3dout));

   
  memcpy((void*)output->nwind,(const void*)input->nwind,sizeof(int)*input->ndims);
  memcpy((void*)output->fwind,(const void*)input->fwind,sizeof(int)*input->ndims);
  memcpy((void*)output->jwind,(const void*)input->jwind,sizeof(int)*input->ndims);

return(SUCCESS);
}

int sep3d_fileexist(char *tag, int ilocal){
int i;

if(sep_thread_num()==0 || ilocal==1) i=tag_exists(tag);

if(ilocal==0) sep3d_broadcast_ints(&i,1,0);
   


return(i);

}




