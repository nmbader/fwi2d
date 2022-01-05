#define SET_SDOC 1
#define SEP3d_AUX_I 1
#include <seplib.h>
#include <superset.h>
#ifdef SU_SUPPORT
#include "superset_internal.h"
#ifndef YES
#define YES 1
#endif
#ifndef NO
#define NO 0
#endif

/*
susep_input_init

USAGE
C<ierr= susep_input_init(tag)>


INPUT PARAMETERS
tag        -  char*  tag of susep structure


RETURN VALUES
SUCCESS      =  if successful

DESCRIPTION
Initialize an input tag

>*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int susep_input_init(char *sep3dname)
_XFUNCPROTOEND
#else
int susep_input_init(sep3dname)
char *sep3dname;
#endif
{
int i1,ierr,i2,i,ierr2;
sep_3d *info;
char temp_ch[1024],temp2_ch[128],temp3_ch[128];
Value *val;
unsigned int *int_pt;
int axes_out[2];

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));


ierr=sep3d_tag_init(sep3dname,sep3dname,"INPUT");
if(ierr!=SUCCESS) return(ierr);

info->su_input=1;
info->su_ndims=info->ndims;
info->n_su=(int*)alloc(sizeof(int)*info->su_ndims);
info->o_su=(float*)alloc(sizeof(float)*info->su_ndims);
info->d_su=(float*)alloc(sizeof(float)*info->su_ndims);
for(i=0;i < info->ndims; i++){
  info->n_su[i]=info->n[i];
  info->o_su[i]=info->o[i];
  info->d_su[i]=info->d[i];
}
	

axes_out[0]=1; axes_out[1]=info->ndims;
ierr=sep3d_change_dims(sep3dname,2,axes_out);
if(ierr!=SUCCESS) return(ierr);



if(info->file_format!=HEADER){
	return(sepwarn(INVALID_DATA,
   "invalid file type, must supply type HEADER as input (%s)\n",sep3dname));
}



/* ALLOCATE THE TEMPORARY STORAGE BLOCKS */

/*see if user has specified the ammount */
sprintf(temp_ch,"%s.nmem",sep3dname);
if(0==getch(temp_ch,"d",&(info->nmem)))
  if(0==getch("nmem","d",&(info->nmem))) info->nmem=100;

if(info->su_hdr_block!=SEPNULL) free(info->su_hdr_block);
info->su_hdr_block=(char*) alloc(HDRBYTES*info->nmem);
if(info->su_tr_block!=SEPNULL) free(info->su_hdr_block);
info->su_tr_block=(char*) alloc(info->n[0]*sizeof(float)*info->nmem);


/* NOW WE NEED TO SET UP THE CONVERSIONS */
for(i1=0; i1 < SU_NKEYS; i1++){
  ierr=getch(hdr[i1].key,"s",temp_ch); /*see if this key is being aliased */
  if(ierr!=1){ /*if it isn't look for the key in the hff */
			ierr2=sep3d_grab_key_index(sep3dname,hdr[i1].key,&(info->key_index[i1]));
			if(ierr2<SUCCESS) return(sepwarn(FAIL_OTHER,
				"trouble grabbing key index for %s, tag %s\n",hdr[i1].key,sep3dname));
			else if (ierr2!=SUCCESS) info->key_index[i1]=0;
  }
  else if(0!=sep3d_grab_key_index(sep3dname,temp_ch,&(info->key_index[i1])))
		 return(sepwarn(FAIL_OTHER,
        "trouble grabbing key index for %s, tag %s\n",temp_ch,sep3dname));

  if(info->key_index[i1]!=0){ /*the key exists in the hff, see if there is 
                               a conversion factor*/
    sprintf(temp_ch,"%s-%s",hdr[i1].key,"fract");
    if(1!=getch(temp_ch,"f",&(info->float_fract[i1]))){
      if(SUCCESS!=sep3d_grab_key(sep3dname,info->key_index[i1],temp_ch,
        temp2_ch,temp3_ch)) return(sepwarn(FAIL_OTHER,
          "trouble getting key info for index %d, tag %s \n",
         info->key_index[i1], sep3dname));

      if(0==strcmp(temp2_ch,"scalar_float")) info->float_fract[i1]=1.;

    }
  }
}

/* INITIALIZE THE HEADER BLOCK TO 0 */
for(i2=0; i2 < info->nmem;i2++){
  for(i=0; i < HDRBYTES; i++){
   info->su_hdr_block[i+i2*HDRBYTES]=(unsigned char) 0;
  }
}

val=(Value*) alloc(sizeof(Value));

/* now set the constant keys */
val->f=info->o[0];
for (i1=0;i1 <info->nmem; i1++)puthval((segy*)(info->su_hdr_block+HDRBYTES*i1),72,val);
val->u=(int)(info->o[0]*1000);
for (i1=0;i1 <info->nmem; i1++)puthval((segy*)(info->su_hdr_block+HDRBYTES*i1),35,val);
      val->f=info->d[0];
for (i1=0;i1 <info->nmem; i1++)puthval((segy*)(info->su_hdr_block+HDRBYTES*i1),71,val);
val->u=(int)(info->d[0]*1000000);
for (i1=0;i1 <info->nmem; i1++)puthval((segy*)(info->su_hdr_block+HDRBYTES*i1),39,val);
val->u=info->n[0];
for (i1=0;i1 <info->nmem; i1++)puthval((segy*)(info->su_hdr_block+HDRBYTES*i1),38,val);
/*int_pt=(unsigned int*)(info->su_hdr_block+114);*/
val->f=info->d[1];
for (i1=0;i1 <info->nmem; i1++)puthval((segy*)(info->su_hdr_block+HDRBYTES*i1),73,val);
val->f=info->o[1];
for (i1=0;i1 <info->nmem; i1++)puthval((segy*)(info->su_hdr_block+HDRBYTES*i1),74,val);
val->i=info->ntraces;
for (i1=0;i1 <info->nmem; i1++)puthval((segy*)(info->su_hdr_block+HDRBYTES*i1),77,val);


free(val);
return(SUCCESS);
}

/*--------------------------------------------------------------------------*/
/*<
susep_output_init

USAGE
ierr= susep_output_init(tag)


INPUT PARAMETERS
tag        -  char*  tag of susep structure


RETURN VALUES
0      =  end of traces

DESCRIPTION
Initialize an output tag

>*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int susep_output_init(char *sep3dname, segy *trace)
_XFUNCPROTOEND
#else
int susep_output_init(sep3dname,trace)
char *sep3dname;
segy *traces;
#endif
{
int i1,ierr,i;
sep_3d *info;
char temp_ch[1024],temp2_ch[128],temp3_ch[128];
float o,d;
int n,i2,is_out;
Value val;



info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

if(0==strcmp(sep3dname,"out")) is_out=YES;

if(0!=sep3d_par_init(sep3dname,"OUTPUT"))
  seperr("trouble initialzing sep3d tag %s \n",sep3dname);

sep3d_set_file_type(sep3dname,"HEADER");
sep3d_set_data_type(sep3dname,"FLOAT");

if(0!= sep3d_set_inorder(sep3dname))
 seperr("trouble setting trace in order flag \n");

/* FIRST THE BASICS */
if(0!=sep3d_set_ndims(sep3dname,2))
  seperr("trouble setting number of dimensions \n");

/*gethval(trace,72,&val); info->o[0]=(float)val.u /1000000.; o=info->o[0];*/
/*gethval(trace,71,&val); info->d[0]=(float)val.u /1000000.; d=info->d[0];*/
/*gethval(trace,38,&val); info->n[0]=(int)val.u; n=info->n[0];*/
/**/
info->o[0]=(float)trace->delrt/1000.;
info->d[0]=(float)trace->dt/1000000.;
info->n[0]=(int)trace->ns;

info->o[1]=(float)trace->f2;
info->d[1]=(float)trace->d2;
info->n[1]=MAX(1,(int)trace->ntr);

sep3d_set_axis(sep3dname,1,info->n[0],info->o[0],info->d[0],"none","none");
sep3d_set_axis(sep3dname,2,info->n[1],info->o[1],info->d[1],"none","none");

/*gethval(trace,73,&val); info->d[1]=(float)val.f; o=info->o[1];*/
/*gethval(trace,74,&val); info->o[1]=(float)val.f; d=info->d[1];*/
/*gethval(trace,77,&val); info->n[1]=(int)val.i; n=info->n[1];*/
/**/


/*sep3d_set_axis(sep3dname,2,n,o,d,"none","none");*/
      sep3d_set_ntraces(sep3dname,info->n[1]);

/* ALLOCATE THE TEMPORARY STORAGE BLOCKS */
sprintf(temp_ch,"%s.nmem",sep3dname);

if(0==getch(temp_ch,"d",&(info->nmem)))
  if(0==getch("nmem","d",&(info->nmem))) info->nmem=100;

if(info->su_hdr_block!=SEPNULL) free(info->su_hdr_block);
info->su_hdr_block=(char*) alloc(HDRBYTES*info->nmem);
if(info->su_tr_block!=SEPNULL) free(info->su_hdr_block);
info->su_tr_block=(char*) alloc(info->n[0]*sizeof(float)*info->nmem);

for(i2=0; i2 < info->nmem; i2++){
  for(i=0; i < HDRBYTES; i++){
    info->su_hdr_block[i+i2*HDRBYTES]=(unsigned char) 0;
  }
}




return(SUCCESS);
}

/*$

=head1 NAME

finish_susep - finish i/o for a SU like program using sep3d data

=head1 SYNOPSIS

finish_susep()

=head1 DESCRIPTION

Correct the number of traces in the sep3d output if it has
changed from input.

=head1 LIBRARY

B<sepsu>

=cut


*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int finish_susep( void)
_XFUNCPROTOEND
#else
int finish_susep()
#endif
{
int n,ierr,suoutput;
float o,d;
char temp_ch[128],temp2_ch[128];
sep_3d *outp;
	suoutput=0;
	getch("suoutput","d",&suoutput);
	

	
	if(suoutput==1) outp=SEPNULL;
  else outp=tag_info_sep3d("out",INQUIRE);
  if(outp!=SEPNULL){ /* if we actually wrote somethin out */
    if(outp->extra_keys==EXTRA_KEY_DEFAULT){ /*we haven't written to the tag*/
      outp->end_trace=outp->last_trace;
      if(SUCCESS!=susep_rite_trace_block("out"))
				seperr("(finish_susep) trouble writing trace block \n");
    }
		else if(outp->end_trace != outp->last_trace){ /*we have block still in mem*/
				outp->end_trace=outp->last_trace;	
				if(SUCCESS!=susep_rite_trace_block(outp->name))
					seperr("trouble writing trace block \n");
		}
    if(outp->ntraces!=outp->last_trace+1){ /* we changed the number of traces */
      ierr=sep3d_grab_axis("out",2,&n,&o,&d,temp_ch,temp2_ch);
			if(ierr!=SUCCESS) seperr("(finish_susep) trouble with grab_axis \n");
      ierr=sep3d_set_axis("out",2,outp->last_trace+1,o,d,temp_ch,temp2_ch);
			if(ierr!=SUCCESS) seperr("(finish_susep) trouble with set_axis \n");
      ierr=sep3d_set_ntraces("out",outp->last_trace+1);
			if(ierr!=SUCCESS)seperr("(finish_susep)trouble with sep3d_set_ntraces\n");
      ierr=sep3d_rite_format("out","out");
			if(ierr!=SUCCESS) seperr("(finish_susep) trouble with sep3d_rite_format \n");
    }
  }

return(SUCCESS);
}

/*<
initargs

USAGE
void initargs(int argc char **argv)


INPUT PARAMETERS
argc -  int number of arguments
argv -  char** argument list



DESCRIPTION
Initialize parameters

>*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void both_initargs(int argc, char **argv)
_XFUNCPROTOEND
#else
void both_initargs(int argc, char **argv)
char *susepname;
#endif
{
char datapth[1024],temp_ch[1024];
size_t targc;    /* total number of args   */
char **targv;    /* pointer to arg strings */
int i,len;
int suinput,suoutput;


initpar(argc,argv);
getch_add_string("in.ignore_gff=1");
if(argc>1){
  targc = argc +1;
  targv = (char **) ealloc1(targc, sizeof(char*));
  for( i=0; i < argc; i++){
    targv[i]=(char *)alloc(sizeof(char)*(strlen(argv[i])+1));
    strcpy(targv[i],argv[i]);
  }
  datapath(datapth);
  sprintf(temp_ch,"tmpdir=%s",datapth);
  targv[targc-1]=(char *)alloc(sizeof(char)*(strlen(temp_ch)+1));
  strcpy(targv[targc-1],temp_ch);
  initargs(targc,targv);
}
else initargs(argc,argv);
if(0==getch("suinput","d",&suinput)) suinput=0;
if(0==getch("suoutput","d",&suoutput)) suoutput=0;
if(suinput!=1) init_3d();
else {
	getch_add_string("noheader=yes");
}
if(suoutput==1) getch_add_string("head=/dev/null");
return;
}
#endif
