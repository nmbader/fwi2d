#define SET_SDOC 1
#define NO_DRN -100
#define NBUF 100000
#define YES 1
#define NO 0

#include "superset_internal.h"
#include<string.h>
#include <superset.h> 



/*
sep3d_convert
<
USAGE
ierr=sep3d_nconvert(char *tag)


INPUT PARAMETERS
tag          -   char*    tag to read

OUTPUT PARAMETERS
nsz    -   int*   convert tag size

>
*/

_XFUNCPROTOBEGIN
int SEP3D_convert(sep_3d  *info,int *offset,char *buf)
_XFUNCPROTOEND
{
int nbuf,tempi,i;

nbuf=*offset;

tempi=(strlen(info->name)+1)*sizeof(char);
memcpy((void*)  (buf+nbuf), (const void*) &tempi, sizeof(int));
nbuf+=sizeof(int);                                            
memcpy((void*)  (buf+nbuf), (const void*) info->name, tempi);
nbuf+=tempi;

switch(info->data_format){
  case(FLOAT): tempi=0; break;
  case(INTEGER): tempi=1; break;
  case(COMPLEX):    tempi=2; break;
  case(BYTE):    tempi=3; break;
  case(UNKNOWN): tempi=4; break;
}
memcpy((void*)  (buf+nbuf), (const void*) &tempi, sizeof(int));
nbuf+=sizeof(int);                                      /*data_format*/

switch(info->file_format){
  case(REGULAR): tempi=0; break;
  case(HEADER): tempi=1; break;
  case(GRID):    tempi=2; break;
  case(UNSPECIFIED): tempi=3; break;
}
memcpy((void*)  (buf+nbuf), (const void*) &tempi, sizeof(int));
nbuf+=sizeof(int);                                      /*file_format*/

switch(info->usage){
  case(INPUT): tempi=0; break;
  case(OUTPUT): tempi=1; break;
  case(SCRATCH):    tempi=2; break;
  case(INQUIRE): tempi=3; break;
}
memcpy((void*)  (buf+nbuf), (const void*) &tempi, sizeof(int));
nbuf+=sizeof(int);                                      /*usage*/


memcpy((void*)  (buf+nbuf), (const void*) &(info->ndims), sizeof(int));
nbuf+=sizeof(int);                                      /*ndims*/

memcpy((void*)  (buf+nbuf), (const void*) (info->n), sizeof(int)*info->ndims);
nbuf+=sizeof(int)*info->ndims;                          /*n*/

memcpy((void*)  (buf+nbuf), (const void*) (info->o), sizeof(float)*info->ndims);
nbuf+=sizeof(int)*info->ndims;                          /*o*/

memcpy((void*)  (buf+nbuf), (const void*) (info->d), sizeof(float)*info->ndims);
nbuf+=sizeof(int)*info->ndims;                          /*d*/

for(i=0; i< info->ndims; i++){
  tempi=(strlen(info->label[i])+1)*sizeof(char);
  memcpy((void*)  (buf+nbuf), (const void*) &tempi, sizeof(int));
  nbuf+=sizeof(int);
  memcpy((void*)  (buf+nbuf),(const void*)info->label[i],tempi);
  nbuf+=tempi;  /*label*/
}
for(i=0; i< info->ndims; i++){
  tempi=(strlen(info->unit[i])+1)*sizeof(char);
  memcpy((void*)  (buf+nbuf), (const void*) &tempi, sizeof(int));
  nbuf+=sizeof(int);
  memcpy((void*)  (buf+nbuf),(const void*)info->unit[i],tempi);
  nbuf+=tempi;  /*unit*/
}

memcpy((void*)  (buf+nbuf), (const void*) (info->nwind), sizeof(int)*info->ndims);
nbuf+=sizeof(int)*info->ndims;                          /*nwind*/
memcpy((void*)  (buf+nbuf), (const void*) (info->fwind), sizeof(int)*info->ndims);
nbuf+=sizeof(int)*info->ndims;                          /*fwind*/
memcpy((void*)  (buf+nbuf), (const void*) (info->jwind), sizeof(int)*info->ndims);
nbuf+=sizeof(int)*info->ndims;                          /*jwind*/
  
memcpy((void*)  (buf+nbuf), (const void*) &(info->nkeys), sizeof(int));
nbuf+=sizeof(int);                                      /*nkeys*/
if(info->nkeys >0){
for(i=0; i< info->nkeys; i++){
  tempi=(strlen(info->keyname[i])+1)*sizeof(char);
  memcpy((void*)  (buf+nbuf), (const void*) &tempi, sizeof(int));
  nbuf+=sizeof(int);
  memcpy((void*)  (buf+nbuf),(const void*)info->keyname[i],tempi);
  nbuf+=tempi;  /*keyname*/
}
for(i=0; i< info->nkeys; i++){
  tempi=(strlen(info->keytype[i])+1)*sizeof(char);
  memcpy((void*)  (buf+nbuf), (const void*) &tempi, sizeof(int));
  nbuf+=sizeof(int);
  memcpy((void*)  (buf+nbuf),(const void*)info->keytype[i],tempi);
  nbuf+=tempi;  /*keytype*/
}
for(i=0; i< info->nkeys; i++){
  tempi=(strlen(info->keyfmt[i])+1)*sizeof(char);
  memcpy((void*)  (buf+nbuf), (const void*) &tempi, sizeof(int));
  nbuf+=sizeof(int);
  memcpy((void*)  (buf+nbuf),(const void*)info->keyfmt[i],tempi);
  nbuf+=tempi;  /*keyfmt*/
}

memcpy((void*)  (buf+nbuf), (const void*) &(info->nh), sizeof(int));
nbuf+=sizeof(int);                                      /*nh*/
memcpy((void*)  (buf+nbuf), (const void*) (info->headers), sizeof(int)*info->nkeys*info->nh);
nbuf+=sizeof(int)*info->nkeys*info->nh;                 /*headers*/

memcpy((void*)  (buf+nbuf), (const void*) (info->drn), sizeof(int)*info->nh);
nbuf+=sizeof(int)*info->nh;                             /*drn*/
}

memcpy((void*)  (buf+nbuf), (const void*) &(info->ncoord), sizeof(int));
nbuf+=sizeof(int);                                      /*ncoord*/
memcpy((void*)  (buf+nbuf), (const void*) (info->coord), sizeof(long long)*info->ncoord);
nbuf+=sizeof(long long)*info->ncoord;                         /*coord*/

memcpy((void*)  (buf+nbuf), (const void*) &(info->ntraces), sizeof(int));
nbuf+=sizeof(int);                                      /*ntraces*/
memcpy((void*)  (buf+nbuf), (const void*) &(info->nextra_keys), sizeof(int));
nbuf+=sizeof(int);                                      /*nextra_keys*/
memcpy((void*)  (buf+nbuf), (const void*) &(info->nkeys_in), sizeof(int));
nbuf+=sizeof(int);                                      /*nkeys_in*/

for(i=0; i< info->nextra_keys; i++){
  tempi=(strlen(info->exp[i])+1)*sizeof(char);
  memcpy((void*)  (buf+nbuf), (const void*) &tempi, sizeof(int));
  nbuf+=sizeof(int);
  memcpy((void*)  (buf+nbuf),(const void*)info->exp[i],tempi);
  nbuf+=tempi;  /*exp*/
}

memcpy((void*)  (buf+nbuf), (const void*) &(info->wrote_data), sizeof(int));
nbuf+=sizeof(int);                                      /*wrote_data*/
memcpy((void*)  (buf+nbuf), (const void*) &(info->wrote_headers), sizeof(int));
nbuf+=sizeof(int);                                      /*wrote_headers*/
memcpy((void*)  (buf+nbuf), (const void*) &(info->su_input), sizeof(int));
nbuf+=sizeof(int);                                      /*su_input*/
                                                                                
if(info->su_input==1){
#ifdef SU_SUPPORT

memcpy((void*)  (buf+nbuf), (const void*) (info->key_index), sizeof(int)*(SU_NKEYS+MAX_EXTRA_KEYS));
nbuf+=sizeof(int)*(SU_NKEYS+MAX_EXTRA_KEYS);                         /*key_index*/
memcpy((void*)  (buf+nbuf), (const void*) (info->float_fract), sizeof(float)*(SU_NKEYS+MAX_EXTRA_KEYS));
nbuf+=sizeof(float)*(SU_NKEYS+MAX_EXTRA_KEYS);                       /*float_fract*/

memcpy((void*)  (buf+nbuf), (const void*) &(info->nmem), sizeof(int));
nbuf+=sizeof(int);                                       /*nmem*/

memcpy((void*)  (buf+nbuf), (const void*) (info->su_tr_block), sizeof(float)*(info->n[0]*info->nmem));
nbuf+=sizeof(float)*info->n[0]*info->nmem;                     /*su_tr_block*/
memcpy((void*)  (buf+nbuf), (const void*) (info->su_hdr_block), HDRBYTES*info->nmem);
nbuf+=HDRBYTES*info->nmem;                                     /*su_hdr_block*/

memcpy((void*)  (buf+nbuf), (const void*) &(info->extra_keys), sizeof(int));
nbuf+=sizeof(int);                                       /*extra_keys*/
memcpy((void*)  (buf+nbuf), (const void*) &(info->beg_trace), sizeof(int));
nbuf+=sizeof(int);                                       /*beg_trace*/
memcpy((void*)  (buf+nbuf), (const void*) &(info->last_trace), sizeof(int));
nbuf+=sizeof(int);                                       /*last_trace*/
memcpy((void*)  (buf+nbuf), (const void*) &(info->nextra_mapping), sizeof(int));
nbuf+=sizeof(int);                                       /*nextra_mapping*/

for(i=0; i < MAX_EXTRA_KEYS;i++){
   memcpy((void*)  (buf+nbuf), (const void*) &(info->extra_name[i]), sizeof(char)*257);
    nbuf+=sizeof(char)*257;                               /*extra name*/
}
memcpy((void*)  (buf+nbuf), (const void*) (info->extra_offset), sizeof(int)*MAX_EXTRA_KEYS);
nbuf+=sizeof(int)*MAX_EXTRA_KEYS;                                    /*extra offset*/
for(i=0; i < MAX_EXTRA_KEYS;i++){
  memcpy((void*)  (buf+nbuf), (const void*) &(info->extra_type[i]), sizeof(char)*5);
  nbuf+=sizeof(char)*5;                                 /*extra type*/
}

memcpy((void*)  (buf+nbuf), (const void*) &(info->su_ndims), sizeof(int));
nbuf+=sizeof(int);                                       /*su_ndims*/
memcpy((void*)  (buf+nbuf), (const void*) (info->n_su), sizeof(int)*info->su_ndims);
nbuf+=sizeof(int)*info->su_ndims;                         /*n_su*/
memcpy((void*)  (buf+nbuf), (const void*) (info->o_su), sizeof(float)*info->su_ndims);
nbuf+=sizeof(float)*info->su_ndims;                       /*o_su*/
memcpy((void*)  (buf+nbuf), (const void*) (info->d_su), sizeof(float)*info->su_ndims);
nbuf+=sizeof(float)*info->su_ndims;                       /*d_su*/

#endif
}


   *offset=nbuf ;

return(0);
}




_XFUNCPROTOBEGIN
int SEP3D_unconvert(sep_3d  *info,int *offset,char *buf)
_XFUNCPROTOEND
{
int nbuf,tempi,i;
nbuf=*offset;

SEP3D_clean(info);
memcpy((void*) &tempi, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                            
info->name=(char*)malloc(sizeof(char)*tempi);
memcpy((void*)  (info->name), (const void*) (buf+nbuf), tempi);
nbuf+=tempi;



memcpy((void*) &tempi, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                            
switch(tempi){
  case(0): info->data_format=FLOAT; break;
  case(1): info->data_format=INTEGER; break;
  case(2):    info->data_format=COMPLEX; break;
  case(3):    info->data_format=BYTE; break;
  case(4): info->data_format=UNKNOWN; break;
}

memcpy((void*) &tempi, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                            
switch(tempi){
  case(0): info->file_format=REGULAR; break;
  case(1): info->file_format=HEADER; break;
  case(2):    info->file_format=GRID; break;
  case(4): info->file_format=UNSPECIFIED; break;
}


memcpy((void*) &tempi, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                            
switch(tempi){
  case(0): info->usage=INPUT; break;
  case(1): info->usage=OUTPUT; break;
  case(2):    info->usage=SCRATCH; break;
  case(4): info->usage=INQUIRE; break;
}


memcpy((void*) &info->ndims, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                      /*ndims*/


if(info->ndims>0){
info->n=(int*)malloc(sizeof(int)*info->ndims);
memcpy((void*) info->n, (const void*) (buf+nbuf), sizeof(int)*info->ndims);
nbuf+=sizeof(int)*info->ndims;                          /*n*/

info->o=(float*)malloc(sizeof(float)*info->ndims);
memcpy((void*) info->o, (const void*) (buf+nbuf), sizeof(float)*info->ndims);
nbuf+=sizeof(float)*info->ndims;                          /*o*/

info->d=(float*)malloc(sizeof(float)*info->ndims);
memcpy((void*) info->d, (const void*) (buf+nbuf), sizeof(float)*info->ndims);
nbuf+=sizeof(float)*info->ndims;                          /*d*/


info->label=(char**)malloc(sizeof(char*)*info->ndims);
for(i=0; i< info->ndims; i++){
	memcpy((void*) &tempi, (const void*) (buf+nbuf), sizeof(int));
	nbuf+=sizeof(int);                                            
  info->label[i]=(char*)malloc(sizeof(char)*tempi);
  memcpy((void*)  (info->label[i]), (const void*) (buf+nbuf), tempi);
  nbuf+=tempi;  /*label*/
}

info->unit=(char**)malloc(sizeof(char*)*info->ndims);
for(i=0; i< info->ndims; i++){
	memcpy((void*) &tempi, (const void*) (buf+nbuf), sizeof(int));
	nbuf+=sizeof(int);                                            
  info->unit[i]=(char*)malloc(sizeof(char)*tempi);
  memcpy((void*)  (info->unit[i]), (const void*) (buf+nbuf), tempi);
  nbuf+=tempi;  /*unit*/
}



info->nwind=(int*)malloc(sizeof(int)*info->ndims);
memcpy((void*) info->nwind, (const void*) (buf+nbuf), sizeof(int)*info->ndims);
nbuf+=sizeof(int)*info->ndims;                          /*nwind*/

info->fwind=(int*)malloc(sizeof(int)*info->ndims);
memcpy((void*) info->fwind, (const void*) (buf+nbuf), sizeof(int)*info->ndims);
nbuf+=sizeof(int)*info->ndims;                          /*fwind*/

info->jwind=(int*)malloc(sizeof(int)*info->ndims);
memcpy((void*) info->jwind, (const void*) (buf+nbuf), sizeof(int)*info->ndims);
nbuf+=sizeof(int)*info->ndims;                          /*jwind*/
}

memcpy((void*) &info->nkeys, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                      /*nkeys*/
if(info->nkeys>0){

info->keyname=(char**)malloc(sizeof(char*)*info->nkeys);
for(i=0; i< info->nkeys; i++){
	memcpy((void*) &tempi, (const void*) (buf+nbuf), sizeof(int));
	nbuf+=sizeof(int);                                            
  info->keyname[i]=(char*)malloc(sizeof(char)*tempi);
  memcpy((void*)  (info->keyname[i]), (const void*)(buf+nbuf), tempi);
  nbuf+=tempi;  /*keyname*/
}

info->keytype=(char**)malloc(sizeof(char*)*info->nkeys);
for(i=0; i< info->nkeys; i++){
	memcpy((void*) &tempi, (const void*) (buf+nbuf), sizeof(int));
	nbuf+=sizeof(int);                                            
  info->keytype[i]=(char*)malloc(sizeof(char)*tempi);
  memcpy((void*)  (info->keytype[i]), (const void*) (buf+nbuf), tempi);
  nbuf+=tempi;  /*keytype*/
}

info->keyfmt=(char**)malloc(sizeof(char*)*info->nkeys);
for(i=0; i< info->nkeys; i++){
	memcpy((void*) &tempi, (const void*) (buf+nbuf), sizeof(int));
	nbuf+=sizeof(int);                                            
  info->keyfmt[i]=(char*)malloc(sizeof(char)*tempi);
  memcpy((void*)  (info->keyfmt[i]), (const void*) (buf+nbuf), tempi);
  nbuf+=tempi;  /*keyfmt*/
}


memcpy((void*) &info->nh, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                      /*nh*/
if(info->nh>0){
info->headers=(int*)malloc(sizeof(int)*info->nh*info->nkeys);
memcpy((void*) info->headers, (const void*) (buf+nbuf), sizeof(int)*info->nkeys*info->nh);
nbuf+=sizeof(int)*info->nkeys*info->nh;                          /*headers*/


info->drn=(int*)malloc(sizeof(int)*info->nh);
memcpy((void*) info->drn, (const void*) (buf+nbuf), sizeof(int)*info->nh);
nbuf+=sizeof(int)*info->nh;                          /*drn*/
}

}


memcpy((void*) &info->ncoord, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                      /*ncoord*/
if(info->ncoord>0){
info->coord=(long long*)malloc(sizeof(long long)*info->ncoord);
memcpy((void*) info->coord, (const void*) (buf+nbuf), sizeof(long long)*info->ncoord);
nbuf+=sizeof(long long)*info->ncoord;                          /*coord*/
}

memcpy((void*) &info->ntraces, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                      /*ntraces*/
memcpy((void*) &info->nextra_keys, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                      /*nextra_keys*/
memcpy((void*) &info->nkeys_in, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                      /*nkeys_in*/

if(info->nextra_keys>0){
info->exp=(char**)malloc(sizeof(char*)*info->nextra_keys);
for(i=0; i< info->nextra_keys; i++){
	memcpy((void*) &tempi, (const void*) (buf+nbuf), sizeof(int));
	nbuf+=sizeof(int);                                            
  info->exp[i]=(char*)malloc(sizeof(char)*tempi);
  memcpy((void*)  (info->exp[i]), (const void*) (buf+nbuf), tempi);
  nbuf+=tempi;  /*keyfmt*/
}
}

memcpy((void*) &info->wrote_data, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                      /*wrote_data*/
memcpy((void*) &info->wrote_headers, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                      /*wrote_headers*/
memcpy((void*) &info->su_input, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                      /*su_input*/
                                                                                
if(info->su_input==1){
#ifdef SU_SUPPORT

memcpy((void*) info->key_index, (const void*) (buf+nbuf), sizeof(int)*(SU_NKEYS+MAX_EXTRA_KEYS));
nbuf+=sizeof(int)*(SU_NKEYS+MAX_EXTRA_KEYS);                         /*key_index*/
memcpy((void*) info->float_fract, (const void*) (buf+nbuf), sizeof(float)*(SU_NKEYS+MAX_EXTRA_KEYS));
nbuf+=sizeof(float)*(SU_NKEYS+MAX_EXTRA_KEYS);                       /*float_fract*/

memcpy((void*) &info->nmem, (const void*) (buf+nbuf), sizeof(int));
nbuf+=sizeof(int);                                       /*nmem*/

if(info->nmem>0){
info->su_tr_block=(char*)malloc(sizeof(float)*(info->n[0]*info->nmem));
memcpy((void*) info->su_tr_block, (const void*) (buf+nbuf), sizeof(float)*(info->n[0]*info->nmem));
nbuf+=sizeof(float)*info->n[0]*info->nmem;                     /*su_tr_block*/
info->su_hdr_block=(char*)malloc(HDRBYTES*info->nmem);
memcpy((void*) info->su_hdr_block, (const void*) (buf+nbuf), (HDRBYTES*info->nmem));
nbuf+=HDRBYTES*info->nmem;                                     /*su_hdr_block*/
}

memcpy((void*) &info->extra_keys, (const void*) (buf+nbuf),sizeof(int));
nbuf+=sizeof(int);                                       /*extra_keys*/
memcpy((void*) &info->beg_trace, (const void*) (buf+nbuf),sizeof(int));
nbuf+=sizeof(int);                                       /*beg_trace*/
memcpy((void*) &info->last_trace, (const void*) (buf+nbuf),sizeof(int));
nbuf+=sizeof(int);                                       /*last_trace*/
memcpy((void*) &info->nextra_mapping, (const void*) (buf+nbuf),sizeof(int));
nbuf+=sizeof(int);                                       /*nextra_mapping*/

for(i=0; i < MAX_EXTRA_KEYS;i++){
   memcpy((void*)  (buf+nbuf), (const void*) (info->extra_name[i]), sizeof(char)*257);
    nbuf+=sizeof(char)*257;                               /*extra name*/
}
memcpy((void*)  (buf+nbuf), (const void*) (info->extra_offset), sizeof(int)*MAX_EXTRA_KEYS);
nbuf+=sizeof(int)*MAX_EXTRA_KEYS;                                    /*extra offset*/
for(i=0; i < MAX_EXTRA_KEYS;i++){
  memcpy((void*)  (buf+nbuf), (const void*) (info->extra_type[i]), sizeof(char)*5);
  nbuf+=sizeof(char)*5;                                 /*extra type*/
}

memcpy((void*)  (buf+nbuf), (const void*) &(info->su_ndims), sizeof(int));
nbuf+=sizeof(int);                                       /*su_ndims*/
if(info->su_ndims>0){
info->n_su=(int*)malloc(sizeof(int)*info->su_ndims);
info->o_su=(float*)malloc(sizeof(float)*info->su_ndims);
info->d_su=(float*)malloc(sizeof(float)*info->su_ndims);
memcpy((void*)  (info->n_su), (const void*) (buf+nbuf), sizeof(int)*info->su_ndims);
nbuf+=sizeof(int)*info->su_ndims;                         /*n_su*/
memcpy((void*)  (info->o_su), (const void*) (buf+nbuf), sizeof(float)*info->su_ndims);
nbuf+=sizeof(float)*info->su_ndims;                       /*o_su*/
memcpy((void*)  (info->d_su), (const void*) (buf+nbuf), sizeof(float)*info->su_ndims);
nbuf+=sizeof(float)*info->su_ndims;                       /*d_su*/
}
#endif
}

   *offset=nbuf ;

return(0);
}










_XFUNCPROTOBEGIN
int SEP3D_nconvert(sep_3d  *info)
_XFUNCPROTOEND
{
int nbuf,tempi,i;

nbuf=0;


nbuf+=sizeof(char)*(strlen(info->name)+1)+sizeof(int);  /*name*/
nbuf+=sizeof(int);                                     /*data_format*/
nbuf+=sizeof(int);                                      /*file_format*/
nbuf+=sizeof(int);                                      /*usage*/

nbuf+=sizeof(int);                                      /*ndims*/
nbuf+=sizeof(int)*info->ndims;                          /*n*/
nbuf+=sizeof(int)*info->ndims;                          /*o*/
nbuf+=sizeof(int)*info->ndims;                          /*d*/
for(i=0; i< info->ndims; i++)
  nbuf+=sizeof(char)*(strlen(info->label[i])+1)+sizeof(int);  /*label*/
for(i=0; i< info->ndims; i++)
  nbuf+=sizeof(char)*(strlen(info->unit[i])+1)+sizeof(int);  /*unit*/

nbuf+=sizeof(int)*info->ndims;                          /*nwind*/
nbuf+=sizeof(int)*info->ndims;                          /*fwind*/
nbuf+=sizeof(int)*info->ndims;                          /*jwind*/
nbuf+=sizeof(int);                                      /*nkeys*/

if(info->nkeys >0){
for(i=0; i< info->nkeys; i++)
  nbuf+=sizeof(char)*(strlen(info->keyname[i])+1)+sizeof(int);  /*keyname*/
for(i=0; i< info->nkeys; i++)
  nbuf+=sizeof(char)*(strlen(info->keytype[i])+1)+sizeof(int);  /*keytype*/
for(i=0; i< info->nkeys; i++)
  nbuf+=sizeof(char)*(strlen(info->keyfmt[i])+1)+sizeof(int);   /*keyfmt*/

nbuf+=sizeof(int);                                      /*nh*/
nbuf+=sizeof(int)*info->nkeys*info->nh;                 /*headers*/
nbuf+=sizeof(int)*info->nh;                             /*drn*/

}

nbuf+=sizeof(int);                                      /*ncoord*/
nbuf+=sizeof(long long)*info->ncoord;                      /*coord*/

nbuf+=sizeof(int);                                      /*ntraces*/
nbuf+=sizeof(int);                                      /*nextra_keys*/
nbuf+=sizeof(int);                                      /*nkeys_in*/


for(i=0; i< info->nextra_keys; i++)
  nbuf+=sizeof(char)*(strlen(info->exp[i])+1)+sizeof(int);  /*exp*/

nbuf+=sizeof(int);                                      /*wrote_data*/
nbuf+=sizeof(int);                                      /*wrote_headers*/
nbuf+=sizeof(int);                                      /*su_input*/
                                                                                
if(info->su_input==1){
#ifdef SU_SUPPORT

nbuf+=sizeof(int)*(SU_NKEYS+MAX_EXTRA_KEYS);                         /*key_index*/
nbuf+=sizeof(float)*(SU_NKEYS+MAX_EXTRA_KEYS);                       /*float_fract*/

nbuf+=sizeof(int);                                       /*nmem*/

nbuf+=sizeof(float)*info->n[0]*info->nmem;                     /*su_tr_block*/
nbuf+=HDRBYTES*info->nmem;                                     /*su_hdr_block*/

nbuf+=sizeof(int);                                       /*extra_keys*/
nbuf+=sizeof(int);                                       /*beg_trace*/
nbuf+=sizeof(int);                                       /*last_trace*/
nbuf+=sizeof(int);                                       /*nextra_mapping*/

nbuf+=sizeof(char)*MAX_EXTRA_KEYS*257;                               /*extra name*/
nbuf+=sizeof(int)*MAX_EXTRA_KEYS;                                    /*extra offset*/
nbuf+=sizeof(char)*MAX_EXTRA_KEYS*5;                                 /*extra type*/

nbuf+=sizeof(int);                                       /*su_ndims*/
nbuf+=sizeof(int)*info->su_ndims;                         /*n_su*/
nbuf+=sizeof(float)*info->su_ndims;                       /*o_su*/
nbuf+=sizeof(float)*info->su_ndims;                       /*d_su*/
#endif
}


return(nbuf);
}


