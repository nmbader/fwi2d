#define SET_SDOC 1
/*

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Sun Jun 16 14:02:50 PDT 1999

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
/* initial values for these pointers is zero (no streams initialized) */
static sep_3d *sep3dlist=SEPNULL;
static sep_3d *sep3dtail=SEPNULL;



/*<
sep3d_head


Usage
sep_3d *sep3d_head(void)

Return Values
sep_3d   -  first sep_3d in list 

Description

Return the begining of a sep_3d list
>*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
sep_3d * sep3d_head(void) 
_XFUNCPROTOEND
#else
sep_3d * sep3d_head() 
#endif 
{ return(sep3dlist); }

/*<
SEP3D_del


Usage 
sep_3d *SEP3D_del(sep_3d *cur)

Input Paramters
cur  -  sep3d*    sep_3d structure to delete

Output Paramters

Return Values
Returns a pointer to the created structure

Description

Delete structure from list
>*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void SEP3D_del( sep_3d *curr )
_XFUNCPROTOEND
#else
void SEP3D_del( curr )
     sep_3d *curr;
#endif
{
    sep_3d *tmp1,*tmp2;

    tmp1=curr->prev;
    tmp2=curr->next;

    if( curr == sep3dlist ) sep3dlist = tmp2;
    if( curr == sep3dtail ) sep3dtail = tmp1;

    if( tmp1 != 0 ) tmp1->next = tmp2;
    if( tmp2 != 0 ) tmp2->prev = tmp1;


}

/*<
sep3d_start


Usage 
sep_3d *sep3d_addstart(curr)

Input Paramters
curr  -  sep3d*  sep_3d structure to add

Output Paramters

Return Values
Returns a pointer to the first structure

Description
Start a sep_3d list
>*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void sep3d_addstart( sep_3d *curr )
_XFUNCPROTOEND
#else
void sep3d_addstart( curr )
sep_3d *curr;
#endif
{
    if( curr == sep3dlist )  return;

    curr->prev = 0;
    if( sep3dlist != 0 ) sep3dlist->prev=curr;
    curr->next=sep3dlist;
    sep3dlist=curr;
    if( sep3dtail == 0 ) sep3dtail=sep3dlist; /* the first entry */
}


/*<
sep3d_addend


Usage 
sep_3d *sep3d_addend(name)

Input Paramters
curr  -  sep3d*  sep_3d structure to add

Output Paramters

Return Values
Pointer to the structure

Description

Add entry to end of sep_3d list
>*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void sep3d_addend( sep_3d *curr )
_XFUNCPROTOEND
#else
void sep3d_addend( curr )
     sep_3d *curr;
#endif
{
    if( curr == sep3dtail )  return;

    curr->prev=sep3dtail;
    if( sep3dtail != 0 ) sep3dtail->next = curr;
    curr-> next = 0;
    sep3dtail=curr;
    if( sep3dlist == 0 ) sep3dlist = curr;
}

 
/*<
sep3d_new


Usage 
sep_3d *sep3d_new(name, usage)

Input Paramters
curr  -  sep3d*       sep_3d structure to add
usage -  usage_type   purpose (INPUT,OUTPUT,SCRATCH) for type

Output Paramters

Return Values
pointer to the structure

Description

Add entry to sep_3d list
>*/
 
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
sep_3d *sep3d_new( const char *tag, enum usage_type usage )
_XFUNCPROTOEND
#else
sep_3d * sep3d_new( tag, usage )
const  char *tag; 
enum usage_type usage;
#endif
{
  sep_3d *newinf;
	int i;

  /* temporary buffer for strings from par line, max length=1024 */
  char tmppar[1024];

  /* if this is just a query, generate an error */
  if( usage == INQUIRE ){ 
     seperr("Unknown tag \"%s\" \n",tag ); 
   }

  newinf = (sep_3d*) calloc(1,sizeof(sep_3d));

  newinf->name = (char*)calloc(1+(int)strlen(tag),sizeof(char));
  strcpy(newinf->name,tag);


	newinf->usage=usage;
	newinf->data_format=UNKNOWN;
	newinf->file_format=UNSPECIFIED;
	newinf->n=SEPNULL;
	newinf->o=SEPNULL;
	newinf->d=SEPNULL;
	newinf->label=SEPNULL;
	newinf->unit=SEPNULL;
	newinf->nkeys=0;
	newinf->nextra_keys=0;
	newinf->nkeys_in=0;
	newinf->drn=SEPNULL;
	newinf->ntraces=0;
	newinf->ndims=0;
	newinf->keyname=SEPNULL;
	newinf->keytype=SEPNULL;
	newinf->keyfmt=SEPNULL;
	newinf->headers=SEPNULL;
	newinf->in_order=-1;
	newinf->nh=0;
	newinf->ntraces_wrote=0;
	newinf->wrote_data=NO;
	newinf->wrote_headers=NO;
  newinf->su_input=0;
#ifdef SU_SUPPORT
 	newinf->beg_trace=BEG_TR_DEFAULT;
  newinf->end_trace=END_TR_DEFAULT;
  newinf->last_trace=BEG_TR_DEFAULT;
  newinf->su_tr_block=SEPNULL;
  newinf->su_hdr_block=SEPNULL;
  newinf->extra_keys=EXTRA_KEY_DEFAULT;
  newinf->nmem=0;
  for(i=0; i< SU_NKEYS+MAX_EXTRA_KEYS; i++){
      newinf->key_index[i]=0;
      newinf->float_fract[i]=0;
  }
  for(i=0; i < MAX_EXTRA_KEYS; i++){
    newinf->extra_offset[i]=0;
  }
#endif
  newinf->coord=SEPNULL;
  newinf->ncoord=0;

 
return(newinf);
} 


/*<
sep3d_print


Usage 
void sep3d_print(sep3d)

Input Paramters
info   -  sep3d*  sep_3d structure to print information about

Output Paramters

Return Values

Description

Print information in sep3d
>*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void sep3d_print( sep_3d *info )
_XFUNCPROTOEND
#else
void sep3d_print( info )
sep3d* info;
#endif
{
int i1,i;


sepwarn(0," name: %s \n",info->name);

switch( info->usage ){
    case( INPUT ): sepwarn(0," USAGE: INPUT \n"); break;
    case( OUTPUT ): sepwarn(0," USAGE: OUTPUT \n"); break;
    case( SCRATCH ): sepwarn(0," USAGE: SCRATCH \n"); break;
}

switch( info->data_format ){
    case( FLOAT ): sepwarn(0," DATA_TYPE: FLOAT \n"); break;
    case( INTEGER ): sepwarn(0," DATA_TYPE: INTEGER \n"); break;
    case( COMPLEX ): sepwarn(0," DATA_TYPE: COMPLEX \n"); break;
    case( BYTE ): sepwarn(0," DATA_TYPE: BYTE \n"); break;
    case( UNKNOWN ): sepwarn(0," DATA_TYPE: UNKNOWN \n"); break;
}


switch( info->file_format ){
    case( REGULAR ): sepwarn(0," FILE_TYPE: REGULAR \n"); break;
    case( HEADER ): sepwarn(0," FILE_TYPE: HEADER \n"); break;
    case( GRID ): sepwarn(0," FILE_TYPE: GRID \n"); break;
    case( UNSPECIFIED ): sepwarn(0," FILE_TYPE: UNSPECIFIED \n"); break;
}

if(info->ndims >0){
	sepwarn(0,"Number of dimensions %d \n \n",info->ndims);
	sepwarn(0,"Axis N\t O\t D\t LABEL\t UNIT\n");
	sepwarn(0,"---------------------------------------------------------\n");
	for(i1=1;i1<=info->ndims;i1++)
		sepwarn(0," %d \t %d \t %f \t %f \t %s \t %s \n",i1,info->n[i1-1],info->o[i1-1],info->d[i1-1],info->label[i1-1],info->unit[i1-1]);  
}

if(info->nkeys>0){
	sepwarn(0,"Number of keys %d \n \n",info->nkeys);
	sepwarn(0,"Key \t NAME \t TYPE\t FORMAT\t \n");
	sepwarn(0,"---------------------------------------------------------\n");
	for(i1=1;i1<=info->nkeys;i1++)
		sepwarn(0,"%d \t %s \t %s \t %s \n",i1,info->keyname[i1-1],info->keytype[i1-1],info->keyfmt[i1-1]);
}

if(info->su_input==1) sepwarn(0,"This is an SU input file \n");


/*SU */

#ifdef SU_SUPPORT
if(info->su_input==1){
	for(i1=0;i1<=info->su_ndims;i1++)
		sepwarn(0,"Axis%d \t n=%d \t o=%f \t d=%f \n ",i1+1,info->n_su[i1],info->o_su[i1-1],info->d_su[i1-1]);
}



sepwarn(0,"last_trace=%d \n",info->last_trace);
sepwarn(0,"beg_trace=%d end_trace=%d nmem=%d \n",info->beg_trace,info->end_trace,info->nmem);
if(info->su_tr_block!=SEPNULL) 
	sepwarn(0,"trace block has been allocated \n");
if(info->su_hdr_block!=SEPNULL) 
  sepwarn(0,"header block has been allocated \n");
if(info->extra_keys!=-1)sepwarn(0,
  "we have %d non-su keys (copied) for this tag \n",info->extra_keys);
else sepwarn(0,"we have not initialzed the keys  \n");
  for(i=0; i< SU_NKEYS; i++){
    if(info->key_index[i]!=0)
      sepwarn(0,"key %d is assigned to sepkey %d with %f fract \n",
        i,info->key_index[i],info->float_fract[i]);
  }
#endif
sepwarn(0,"Other Information \n --------------------------------------------\n");
if(info->ntraces>0) sepwarn(0,"Number of traces : %d \n",info->ntraces);
else sepwarn(0,"Number of traces not set \n");

if(info->drn==SEPNULL) sepwarn(0,"data_record_number has not been set\n");
else sepwarn(0,"data_record_number exists in the data \n");
if(info->headers!=SEPNULL) {
	sepwarn(0,"Storing portion of headers in structure \n");
	sepwarn(0,"total number of headers %d \n",info->nh);
}
if(info->wrote_data==YES) sepwarn(0,"Have written data \n");
if(info->wrote_headers==YES) sepwarn(0,"Have written headers \n");

}
/*  $Id: list.c,v 1.2 2004/04/08 22:32:28 bob Exp $ */
