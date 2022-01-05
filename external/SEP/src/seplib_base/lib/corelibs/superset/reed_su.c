#define SET_SDOC 1
/*************************************************************/
/*******************SU SUPPORT********************************/
/*************************************************************/

#define YES 1
#define NO 0

#include <seplib.h>
#include <superset.h>
#ifdef SU_SUPPORT
#include "superset_internal.h"
#include <segy.h>

/*local routine*/
int susep_init_reed(sep_3d *info);
void off_puthval(segy *tr,char type, int ioff, Value *valp);

/* internal function prototype*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void setval( cwp_String , Value* , double);
_XFUNCPROTOBEGIN
#else
void setval( );
#endif


/*$

=head1 NAME

tgettr - get the next trace from a SEP3d dataset

=head1 SYNOPSIS

C<ierr= tgettr("in", trace)>


=head1 INPUT PARAMETERS

=over 4

=item sep3dname    -  char*  

      tag of sep_3d structure

=item trace        -  segy*  

      trace to read in

=back

=head1 RETURN VALUES

=over 4

=item returns 0 at end of traces

=item nread  =  n

      number of bytes read

=back

=head1 DESCRIPTION

Read in the next trace

=head1 SEE ALSO

L<tgettra>, L<tputtr>

=head1 LIBRARY

B<sepsu>

=cut

>*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int tgettr(char *sep3dname, segy *trace)
_XFUNCPROTOEND
#else
int tgettr(sep3dname ,trace)
char *sep3dname;
segy *trace;
#endif
{
sep_3d *info;
int sz,off;
int b2,e2,i;
Value val;
int suinput;

if(0==getch("suinput","d",&suinput)) suinput=0;
if(suinput==1) return(fgettr(stdin, (trace)));


if(0==strcmp(sep3dname,"in")) info = tag_info_sep3d(sep3dname, INPUT);
else info = tag_info_sep3d(sep3dname, SCRATCH);  /* get info on this tag */

if(info->beg_trace == BEG_TR_DEFAULT) { /* we need to initialize this tag */

  if(susep_input_init(sep3dname)) /*we must error condition this because ;
																	  it is the upper most level we have 
                                     of */
    seperr("(tgettr) trouble initialize tag %s \n",sep3dname);


   /*we need to allow conversion into the extra part of the segy header*/
   if(suinput==0) {
     if(SUCCESS!=susep_init_reed(info)){
        seperr("trouble initializing susep_init_reed \n");
     }
   }


  info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
  info->beg_trace=BEG_TR_INIT;
}
info->last_trace++; /* update the last trace read */

if(info->last_trace==info->n[1])return(0);/*we have finished reading the data */

if(info->last_trace < info->beg_trace || info->last_trace  > info->end_trace){
  /* we don't have this trace in memory so read in a trace block */
  b2=info->last_trace;
  e2=MIN(b2+info->nmem,info->ntraces)-1;
	if(SUCCESS!= susep_read_trace_block(sep3dname,b2,e2))
		seperr("(tgettr) Unable to read trace block (tag=%s)\n",sep3dname);
	if(SUCCESS!= SEP3D_alloc_coord(info,e2-b2+1))
		seperr("(tgettr) Unable to read trace block (tag=%s)\n",sep3dname);
  for(i=0; i < info->ncoord; i++) info->coord[i]=(double)(b2+i);
}

/*copy the buffered trace in header to segy trace */
sz=info->n[0]*sizeof(float);
off=info->last_trace-info->beg_trace;
memcpy( (void *)trace, (const void *)(info->su_hdr_block+(HDRBYTES*off)),
HDRBYTES);
memcpy( (void *) trace->data, (const void *) (info->su_tr_block+(sz*off)),
sz);


/*Sucess, return the size of the trace */
return(sz+HDRBYTES);

}

/*$

=head1 NAME

tgettr - read a specified  trace from a SEP3d dataset

=head1 SYNOPSIS

C<ierr= tgettr("in", trace,itr)>


=head1 INPUT PARAMETERS

=over 4

=item sep3dname    -  char*  

      tag of sep_3d structure

=item trace        -  segy*  

      trace to read in

=item itr          -  int    

      trace number

=back

=head1 RETURN VALUES

=over 4

=item returns 0 at end of traces

=item nread  =  

      number of bytes read

=back

=head1 DESCRIPTION

Read in  a given trace specified by itr

=head1 SEE ALSO

L<tgettr>, L<tputtr>

=head1 LIBRARY

B<sepsu>

=cut

>*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int tgettra( char *sep3dname,  segy *trace, int itr)
_XFUNCPROTOEND
#else
cwp_String type;
int tgettra(sep3dname,trace,itr)
char *sep3dname;
segy *trace;
int itr;
#endif
{
sep_3d *info;

if(0==strcmp(sep3dname,"in")) info = tag_info_sep3d(sep3dname, INPUT);
else info = tag_info_sep3d(sep3dname, SCRATCH);  /* get info on this tag */

if(info->beg_trace == BEG_TR_DEFAULT) { /* we need to initialize this tag */
  if(0!=susep_input_init(sep3dname))
    seperr("trouble initialize tag %s \n");
  info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
	if(info==SEPNULL) seperr("(tgettra) could not initialize tag %s \n",sep3dname);
}

if(itr != info->last_trace +1){ /* we are not reading in order */
  sepwarn(0,"resetting the nmem to 1 \n");
  info->nmem=1;
  info->last_trace=itr-1;
}
info->beg_trace=BEG_TR_INIT;
return(tgettr(sep3dname,trace));
}
/*<
setval

Usage
setval( cwp_String type, Value *valp, double a)


INPUT PARAMETERS
type   - cwp_String   type of input
a      - double       value to set valp to

OUTPUT PARAMETERS
valp   -  Value*     Value to set


DESCRIPTION
Set the mixed trace variable to the value a 

>*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void setval( cwp_String type, Value *valp, double a)
_XFUNCPROTOEND
#else
void setval(type, valp,a)
cwp_String type;
Value *valp;
double a;
#endif
{
  switch (*type) {
  case 's':
    err("can't set char header word");
  break;
  case 'h':
    valp->h = a;
  break;
  case 'u':
    valp->u = a;
  break;
  case 'l':
    valp->l = (long) (a);
  break;
  case 'v':
    valp->v = (unsigned long) (a);
  break;
  case 'i':
    valp->i = a;
  break;
  case 'p':
    valp->p = a;
  break;
  case 'f':
  valp->f = a;
  break;
  case 'd':
    valp->d = a;
  default:
    err("unknown type %s", type);
  break;
  }
  return;
}
/*<
susep_red_trace_bloack


INPUT PARAMETERS
sep3dname   -   char    tag to read into
beg_trace   -   int     first trace to read
end_trace   -   int     last trace to read

DESCRIPTION
Read a trace block into memory

>*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int susep_read_trace_block(char *sep3dname,int beg_trace, int end_trace)
_XFUNCPROTOEND
#else
int susep_read_trace_block(char sep3dname,beg_trace,end_trace)
char *sep3dname;
int beg_trace,end_trace;
#endif
{
int *tempi,sz,i1,fwind[2],jwind[2],nwind[2],itr;
float *tempr;
sep_3d *info;
double a;
int  ierr,index;
Value val;

info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

info = tag_info_sep3d(sep3dname, INPUT);

if(info->beg_trace == BEG_TR_DEFAULT)
  return(sepwarn(INVALID_STRUC,"tag:%s susep_input_init has not been called \n",
     sep3dname));

  tempi=(int *) alloc(sizeof(int)*info->nmem); /*integer temp block */
  tempr=(float *) alloc(sizeof(float)*info->nmem); /*float temp block */

  /*we need to read a new block */
  info->beg_trace=beg_trace;
  info->end_trace=end_trace;
  sz=info->end_trace-info->beg_trace+1;

  /* read in the header header */
  fwind[0]=info->last_trace; jwind[0]=1; nwind[0]=sz;
	ierr=sep3d_read_headers(sep3dname,sep3dname,nwind,fwind, jwind,&i1);

	if(ierr!=SUCCESS){ 
    sepwarn(0,"trouble reading in headers  for tag %s\n",sep3dname);
    free(tempr); free(tempi); return(ierr);
  }


  /*set the header keys */
  for(i1=0; i1 < SU_NKEYS; i1++){ 
    if(info->key_index[i1] !=0){ /*we have this key */
      if(info->float_fract[i1] !=0.){ /*this is a float key */
				ierr=sep3d_grab_header_vals_i(sep3dname,info->key_index[i1],(int*)tempr);
				if(ierr!=SUCCESS){ free(tempr); free(tempi); return(ierr);}
          for(itr=0; itr< sz; itr++){
            a=(double)(tempr[itr]*info->float_fract[i1]);
            setval(hdr[i1].type,&val,a);
            puthval((segy*)(info->su_hdr_block+itr*HDRBYTES),i1,&val);
          }
      }
      else{
				ierr=sep3d_grab_header_vals_i(sep3dname,info->key_index[i1],(int*)tempi);
				if(ierr!=SUCCESS){ free(tempr); free(tempi); return(ierr);}
        for(itr=0; itr< sz; itr++){
          a=(double)(tempi[itr]);
          setval(hdr[i1].type,&val,a);
          puthval((segy*)(info->su_hdr_block+itr*HDRBYTES),i1,&val);
        }
      }
    }
  }

  /* now do extra mapping if needed*/
for(i1=0; i1 < info->nextra_mapping;i1++){
  /*get key index*/
    if(SUCCESS !=sep3d_grab_key_index(sep3dname,info->extra_name[i1],&index))
        seperr("trouble grabbing key index for %s, tag %s\n",info->extra_name[i1],sep3dname);
    if(0==strcmp(info->keytype[index],"scalar_float")){
				ierr=sep3d_grab_header_vals_i(sep3dname,index,(int*) tempr);
				if(ierr!=SUCCESS){ free(tempr); free(tempi); return(ierr);}
        for(itr=0; itr< sz; itr++){
          a=(double)(tempr[itr]);
          setval(info->extra_type[i1],&val,a);
          off_puthval((segy*)(info->su_hdr_block+itr*HDRBYTES),info->extra_type[i1][0],
           info->extra_offset[i1],&val);
        }
    }
    else if(0==strcmp(info->keytype[index],"scalar_int")){
				ierr=sep3d_grab_header_vals_i(sep3dname,index,tempi);
				if(ierr!=SUCCESS){ free(tempr); free(tempi); return(ierr);}
        for(itr=0; itr< sz; itr++){
          a=(double)(tempi[itr]);
          setval(info->extra_type[i1],&val,a);
          off_puthval((segy*)(info->su_hdr_block+itr*HDRBYTES),info->extra_type[i1][0],
           info->extra_offset[i1],&val);
        }
    }
    else {
      sepwarn(-1,
         "Only know how to convert scalar_float and scalar_int keys with extra mapping\n");
    }

}



/*grab the data record number for the trace*/
  ierr=sep3d_grab_drn(sep3dname,tempi);
	if(ierr!=SUCCESS) return(ierr);
  fwind[0]=0; nwind[0]=info->n[0];
  ierr= sep3d_read_list(sep3dname,nwind[0],nwind[0],0,1,4,sz,tempi,
  info->su_tr_block);
	if(ierr!=SUCCESS){ free(tempr); free(tempi); return(ierr);}
  return(SUCCESS);
}

/*<
susep_init_reed

USAGE
ierr=susep_init_reed(info)

INPUT PARAMETERS
info  -   sep3d*    sep_3d tag to rite initialize

DESCRIPTION
Initialize rite out

>*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int susep_init_reed(sep_3d *info)
_XFUNCPROTOEND
#else
int susep_init_reed(info )
sep_3d *info;
#endif
{
int i2,i1,nk,nkeys,logic,ierr,ierr2,nkey_temp,i,suoutput;
sep_3d *inp;
Value val,zeroval;
char temp_ch[128],temp2_ch[128],temp3_ch[128];
int tempi,nk_in,suinput,nk_new;


if(info==SEPNULL) return(sepwarn(INVALID_DATA,"Invalid sep_3d structure passed \n"));

/*now deal with potential extra mappings*/
if(0==getch("nextra_keys","d",&info->nextra_mapping)) info->nextra_mapping=0;
if(info->nextra_mapping < 0 || info->nextra_mapping>MAX_EXTRA_KEYS)
  seperr("only allow 0 to 15 additional header mapping nextra_keys=%d \n",
   info->nextra_mapping);
for(i=0; i < info->nextra_mapping;i++){
  sprintf(temp_ch,"extra_name%d",i+1);
  if(0==getch(temp_ch,"s",&info->extra_name[i]))
    seperr("%s not specified\n",temp_ch);
  /*check to make sure the key exists*/
   ierr2=sep3d_grab_key_index(info->name,info->extra_name[i],&ierr);
      if(ierr2!=SUCCESS) return(sepwarn(FAIL_OTHER,
        "trouble grabbing key index for %s, tag %s\n",info->name,info->extra_name[i]));

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
}


return(SUCCESS);

}

void off_puthval(segy *tr, char type, int ioff, Value *valp)
{
  char *tp = (char*) tr;


  switch(type) {
  case 's': (void) strcpy(tp + ioff, valp->s);  break;
  case 'h': *((short*)  (tp + ioff)) = valp->h; break;
  case 'u': *((unsigned short*) (tp + ioff)) = valp->u; break;
  case 'i': *((int*)   (tp + ioff)) = valp->i; break;
  case 'p': *((unsigned int*)   (tp + ioff)) = valp->p; break;
  case 'l': *((long*)   (tp + ioff)) = valp->l; break;
  case 'v': *((unsigned long*)  (tp + ioff)) = valp->v; break;
  case 'f': *((float*)  (tp + ioff)) = valp->f; break;
  case 'd': *((double*) (tp + ioff)) = valp->d; break;
  default: seperr("%s: %s: mysterious data type", __FILE__,__LINE__); break;
  }

  return;
}




#endif


/*  $Id: reed_su.c,v 1.2 2004/04/08 22:32:28 bob Exp $ */
