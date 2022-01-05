#define SET_SDOC 1
#include <seplib.h>
#include <superset.h>
#include "superset_internal.h"
#ifdef SU_SUPPORT
/*************************/
/*  SU SUPPORT          */
/*************************/
#include <segy.h>
/*$

=head1 NAME

tputtr -  put given trace next in output stream

=head1 SYNOPSIS

C<ierr= tputtr("in", trace)>

=head1 INPUT PARAMETERS

=over 4

=item sep3dname    -  char*  

      tag of sep_3d structure

=item trace        -  segy*  

      trace to read in

=back

=head1 RETURN VALUES

 0  = if it fails

 1  =  if it succeeds

=head1 DESCRIPTION

Read in the next trace

=head1 SEE ALSO

L<tgettr>

=head1 LIBRARY

B<sepsu>


=cut

>*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int tputtr(char *sep3dname, segy *trace)
_XFUNCPROTOEND
#else
int tputtr(sep3dname ,trace)
char *sep3dname;
segy *trace;
#endif
{
int i1,fwind[2],jwind[2],nwind[2],itr;
sep_3d *info;
Value val;
float *tempr;
int *tempi,cur_trace,sz;
int suoutput;
char *type;
float o,d;
char temp_ch[128],temp2_ch[128];



if(0==getch("suoutput","d",&suoutput)) suoutput=0;
if(suoutput==1){
    fputtr(stdout, trace);
    return(1);
}



if(0==strcmp("out",sep3dname)) info = tag_info_sep3d(sep3dname, OUTPUT);
else  info = tag_info_sep3d(sep3dname, SCRATCH);

if(info->beg_trace==BEG_TR_DEFAULT){/*this is our first call to put_trace */

  if(SUCCESS!=susep_output_init(sep3dname,trace))
    seperr("trouble initialize tag %s \n",sep3dname);


  info->beg_trace=0;
  info->end_trace=MIN(info->beg_trace+ info->nmem, info->ntraces)-1;

}

if(info->last_trace==info->ntraces-1){ /*we  have more traces going out than in*/
  info->ntraces=1234567890;
  info->end_trace=MIN(info->beg_trace+ info->nmem, info->ntraces)-1;
	info->n[1]=info->ntraces;
}
/* put the trace  and header into storage */
cur_trace=info->last_trace+1-info->beg_trace;


memcpy( (void*)(info->su_hdr_block+cur_trace*HDRBYTES), (const void*)trace,
  (size_t) HDRBYTES);
sz=info->n[0]*sizeof(float);
memcpy( (void*)(info->su_tr_block+cur_trace*sz),( void*)trace->data, sz);
if(cur_trace+info->beg_trace==info->end_trace){
   if(0!=susep_rite_trace_block(sep3dname))
     seperr("trouble writing trace block \n");
  info->beg_trace=info->last_trace+2;
/*  info->end_trace=MIN(info->beg_trace+ info->nmem, info->ntraces)-1;*/
  info->end_trace=info->beg_trace+ info->nmem-1;

}

  info->last_trace++;
  return(1);

}



/*< 
susep_init_rite

USAGE 
ierr=susep_init_rite(info)

INPUT PARAMETERS
info  -   sep3d*    sep_3d tag to rite initialize

DESCRIPTION
Initialize rite out

>*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int susep_init_rite(sep_3d *info)
_XFUNCPROTOEND
#else
int susep_init_rite(info )
sep_3d *info;
#endif
{
int i2,i1,nk,nkeys,logic,ierr,ierr2,nkey_temp,i,suoutput;
sep_3d *inp;
Value val,zeroval;
char temp_ch[128],temp2_ch[128],temp3_ch[128];
int tempi,nk_in,suinput,nk_new;





if(info==SEPNULL) return(sepwarn(INVALID_DATA,"Invalid sep_3d structure passed \n"));

sep3d_set_rite_status(info->name,YES,YES);
if(0==strcmp(info->name,"out")){  /*if stdin out first just copy the keys */
	suinput=0;
	getch("suinput","d",&suinput);
	if(suinput!=1){
  	inp=tag_info_sep3d("in",INQUIRE);
  	if(SEPNULL==inp) nk=0; /*input doesn't exis so no extra keys too copy */
 		else{
				nk_in=nk=inp->nkeys;
   		 for(i1=0; i1 < SU_NKEYS; i1++){
     		 if(inp->key_index[i1]!=0){
       		 info->key_index[i1]=inp->key_index[i1];
        		if(inp->float_fract[i1]!=0.)
          		info->float_fract[i1]=1./inp->float_fract[i1];
     		 }
 		   }
 		 }
	}
	else{
		nk=0;
  	inp=SEPNULL;
	}
}
else { /*we are not wrting out to stdout */
  nk=0;
  inp=SEPNULL;
}

nkeys=nk;

zeroval.l=0;
/*check for non-zero keys in current block */
for(i1=0; i1 < SU_NKEYS; i1++){
  if(i1==71 || i1==72|| i1==73 || i1==74 || i1==38 || i1==77|| i1==35 || i1==39)
  ;else{ /*if it isn't a special key (o1,d1,n1,etc) */
  if(info->key_index[i1]!=0) logic=YES;
  else{
    logic=NO;
    i2=0;
    /* loop over curent block */
    while(i2 < info->end_trace-info->beg_trace+1 && logic==NO){ 
      gethval((const segy*)(info->su_hdr_block+HDRBYTES*i2), i1, &val);
      ierr=getch(hdr[i1].key,"s",temp_ch);
      if( (!(ierr==1 && 0==strcmp("skip",temp_ch))) &&
          (ierr==1 || valcmp(hdr[i1].type,val,zeroval))){

         /* we want to map this key */
         if(SEPNULL==inp) logic=YES; /*we will need to create the key */
         else{
           for(i1=0; i1 < nk_in; i1++){
             if(0==strcmp(inp->keyname[i1],hdr[i1].key)){ /*we have a match*/

             info->key_index[i1]=i1+1;
               if(0==strcmp(inp->keytype[i1],"scalar_float")){
                 sprintf(temp_ch,"%s-fract",hdr[i1].key);
                 if(1==getch(temp_ch,"f",&(info->float_fract[i1])))
                   info->float_fract[i1]=1./info->float_fract[i1];
                 else info->float_fract[i1]=1.;
               }
               logic=YES;
             }
           }
         }
         if(logic==YES){ /*we have a new key */
           nkeys++;
           info->key_index[i1]=nkeys;
           sprintf(temp_ch,"%s-fract",hdr[i1].key);
            ierr=getch(temp_ch,"f",&(info->float_fract[i1]));
         }
       }
      i2++;
    }
  }
}
}
if(nkeys==0){
  nkeys++;
  info->key_index[7]=nkeys;
}

/*now deal with potential extra mappings*/
if(0==getch("nextra_keys","d",&info->nextra_mapping)) info->nextra_mapping=0;
if(info->nextra_mapping < 0 || info->nextra_mapping>MAX_EXTRA_KEYS)
	seperr("only allow 0 to MAX_EXTRA_KEYS additional header mapping nextra_keys=%d \n",
   info->nextra_mapping);
for(i=0; i < info->nextra_mapping;i++){
	sprintf(temp_ch,"extra_name%d",i+1);
	if(0==getch(temp_ch,"s",&info->extra_name[i]))
    seperr("%s not specified\n",temp_ch);
	sprintf(temp_ch,"extra_offset%d",i+1);
	if(0==getch(temp_ch,"d",&info->extra_offset[i]))
    seperr("%s not specified\n",temp_ch);
  if(info->extra_offset[i] <0 || info->extra_offset[i] > 238)
    seperr("invalid %d=%s, must be between 0 and 238 \n",temp_ch,
     info->extra_offset[i]);
	sprintf(temp_ch,"extra_type%d",i+1);
	if(0==getch(temp_ch,"s",&info->extra_type[i]))
    seperr("%s not specified\n",temp_ch);
  if(info->extra_type[i][0]!='h' && info->extra_type[i][0]!='u' &&
  info->extra_type[i][0]!='l' && info->extra_type[i][0]!='v' &&
  info->extra_type[i][0]!='i' && info->extra_type[i][0]!='p' &&
  info->extra_type[i][0]!='f' && info->extra_type[i][0]!='d' )
   seperr("unacceptable type %s=%d must be (h,u,l,v,i,p,f, or d)\n");
  nkeys++;
}


ierr=sep3d_set_nkeys(info->name,nkeys);
if(ierr!=SUCCESS) return(ierr);

for(i1=1; i1 <= nk; i1++){
  ierr=sep3d_grab_key("in",i1,temp_ch,temp2_ch,temp3_ch);
	if(ierr!=SUCCESS) return(ierr);
  ierr=sep3d_set_key(info->name,i1,temp_ch,temp2_ch,temp3_ch);
	if(ierr!=SUCCESS) return(ierr);
}
for(i1=0; i1 < SU_NKEYS; i1++){
  if(info->key_index[i1] >nk){
    /*a new key not copied from the input */
    strcpy(temp3_ch,hdr[i1].key);
    getch(hdr[i1].key,"s",temp3_ch);
    if(hdr[i1].type[0]=='f'||hdr[i1].type[0]=='d'|| info->float_fract[i1]!=0.){
     strcpy(temp_ch,"scalar_float");
     strcpy(temp2_ch,"xdr_float");
    }
    else{
     strcpy(temp_ch,"scalar_int");
     strcpy(temp2_ch,"xdr_int");
    }
    ierr=sep3d_set_key(info->name,info->key_index[i1],temp3_ch,
      temp_ch,temp2_ch);
		if(ierr!=SUCCESS) return(ierr);
  }
}
/* now deal with the extra keys*/
for(i1=0; i1 < info->nextra_mapping;i1++){
    if(info->extra_type[i1][0] == 'f' || info->extra_type[i1][0] == 'd' ){
     strcpy(temp_ch,"scalar_float");
     strcpy(temp2_ch,"xdr_float");
    }
    else{
     strcpy(temp_ch,"scalar_int");
     strcpy(temp2_ch,"xdr_int");
    }
    ierr=sep3d_set_key(info->name,info->nkeys-info->nextra_mapping+i1+1,
      info->extra_name[i1], temp_ch,temp2_ch);
		if(ierr!=SUCCESS) return(ierr);
}

if(0==getch("suoutput","d",&suoutput)) suoutput=0;
if(suoutput!=1) ierr=sep3d_rite_format(info->name,info->name);

if(ierr!=SUCCESS) return(ierr);


info->extra_keys=nk;
return(SUCCESS);
}

/*<
susep_rite_trace_bloack


INPUT PARAMETERS
sep3dname -  char*  tag of the structure to write out


DESCRIPTION

Write out a block of traces and headers 


>*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int susep_rite_trace_block(char *sep3dname)
_XFUNCPROTOEND
#else
int susep_rite_trace_block(char tsep3dname)
char *sep3dname;
#endif
{
int *tempi,sz,i1,fwind[2],jwind[2],nwind[2],itr;
float *tempr;
sep_3d *info;
double a;
Value val;
int ierr;
char *tp;
int i;



info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

if(info->beg_trace == BEG_TR_INIT) return(sepwarn(INVALID_DATA,
   "you have not put anything in this trace block (%s) \n",sep3dname));

 if(info->extra_keys==EXTRA_KEY_DEFAULT) ierr=susep_init_rite(info);
 else ierr=0;
if(ierr!=0) return(sepwarn(INVALID_STRUC,"Trouble writing out format \n"));
      

if(info->end_trace-info->beg_trace+1==0) return(0);
  ierr=sep3d_set_nh(sep3dname,info->end_trace-info->beg_trace+1);
  ierr=SEP3D_alloc_coord(info,info->end_trace-info->beg_trace+1);
  for(itr=0; itr<=info->end_trace-info->beg_trace; itr++)
     info->coord[itr]=(double)(info->beg_trace+itr);
    
	if(ierr!=SUCCESS) return(ierr);


  tempi=(int*) alloc(sizeof(int)*info->nmem);

  for(i1=0; i1 < info->extra_keys; i1++){
    ierr=sep3d_grab_header_vals_i("in",i1+1,tempi);
		if(ierr!=SUCCESS){ free(tempi); return(ierr);}
    ierr=sep3d_set_header_vals_i(sep3dname,i1+1,tempi);
		if(ierr!=SUCCESS){ free(tempi); return(ierr);}
  }

  tempr=(float*) alloc(sizeof(float)*info->nmem);


  /*write out the header */
  for(i1=0; i1 < SU_NKEYS; i1++){
    if(info->key_index[i1] !=0){
      if(hdr[i1].type[0]=='f'||hdr[i1].type[0]=='d'||info->float_fract[i1]!=0.){
        for(itr=0; itr<=info->end_trace-info->beg_trace; itr++){
          gethval((const segy*)(info->su_hdr_block+HDRBYTES*itr),i1,&val);
          if(info->float_fract[i1]!=100.)
            tempr[itr]=info->float_fract[i1]*vtof(hdr[i1].type,val);
          else tempr[itr]=vtof(hdr[i1].type,val);
          ierr=sep3d_set_header_vals_i(sep3dname,info->key_index[i1],
           (int*)tempr);
					if(ierr!=SUCCESS){ free(tempr); free(tempi); return(ierr);}
        }
      }
       else{
        for(itr=0; itr<=info->end_trace-info->beg_trace; itr++){
        gethval((const segy*)(info->su_hdr_block+HDRBYTES*itr),i1,&val);
          tempi[itr]=vtoi(hdr[i1].type,val);
        }
  				ierr=sep3d_set_header_vals_i(sep3dname,info->key_index[i1],
            (int*)tempi);
					if(ierr!=SUCCESS){ free(tempr); free(tempi); return(ierr);}
      }
    }
  }

  /*now deal with extra keys */
for(i1=0; i1 < info->nextra_mapping;i1++){
  for(itr=0; itr<=info->end_trace-info->beg_trace; itr++){
    tp=(char*)info->su_hdr_block+HDRBYTES*itr;
    switch(info->extra_type[i1][0]){
      case 'h': tempi[itr]= *((short*)  (tp + info->extra_offset[i1])); break;
      case 'u': tempi[itr]= *((unsigned short*)(tp +info->extra_offset[i1]));break;
      case 'i': tempi[itr]= *((int*)   (tp + info->extra_offset[i1])); break;
      case 'p': tempi[itr]= *((unsigned int*)(tp+info->extra_offset[i1])); break;
      case 'l': tempi[itr]= *((long*)   (tp + info->extra_offset[i1])); break;
      case 'v': tempi[itr]= *((unsigned long*)(tp+info->extra_offset[i1])); break;
      case 'f': tempr[itr]= *((float*)  (tp + info->extra_offset[i1])); break;
      case 'd': tempr[itr]= *((double*) (tp + info->extra_offset[i1])); break;
  	}
	}
  if(info->extra_type[i1][0] == 'f' || info->extra_type[i1][0] == 'd' )
   ierr=sep3d_set_header_vals_i(sep3dname,info->nkeys-info->nextra_mapping+i1+1,
    (int*)tempr);
  else ierr=sep3d_set_header_vals_i(sep3dname,info->nkeys-info->nextra_mapping+i1+1, tempi);
  if(ierr!=SUCCESS){ free(tempr); free(tempi); return(ierr);}
}
    
 fwind[0]=0; jwind[0]=1; nwind[0]=info->ntraces;
 fwind[1]=info->beg_trace; jwind[1]=1; nwind[1]=info->end_trace-info->beg_trace+1;
	if(0!=sep3d_set_inorder(sep3dname))
		seperr("trouble setting traces in order \n");

  ierr=sep3d_rite(sep3dname,sep3dname,nwind,fwind,jwind,
      info->su_tr_block, nwind[1],1,1,0);
	if(ierr!=SUCCESS) {free(tempr); free(tempi); return(ierr);}

  free(tempi);
  free(tempr);


  return(SUCCESS);
}
#endif
