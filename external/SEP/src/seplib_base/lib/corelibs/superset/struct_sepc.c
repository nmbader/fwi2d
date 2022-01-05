#define SET_SDOC 1
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Sat Nov 16 15:43:54 PST 2002


Purpose: To handle axis information in the superset data

*/	 

int sep3d_tag_num=0;
#include<string.h>
#include <superset_internal.h>
#include <sep3dc.h> 
void my_axes_clean(sep3d *sep3dc);
void my_keys_clean(sep3d *sep3dc);
void my_char_dealloc(char **str_array, int nlen);
extern int sep_thread_num(void);

/*
=head1 NAME

sep3d_grab_sep3d - Grab the internal structure


=head1 SYNOPSIS

C<logic= sep3d_grab_sep3d(char *sep3dt,sep3d *sep3dc)>

=head1 INPUT PARAMETERS

=over 4

=item sep3dt  -  char*

     sep3d internal structure tag

=item struct  -  sep3d

      structure


=back

=head1 DESCRIPTION

Grab the internal copy of  a structure to the public sep3dc

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

#ifdef SEP_MPI
#include<mpi.h>
#endif

int sep3d_grab_sep3d(char *sep3dt, sep3d *sep3dc){
int i,i1;



if(1==valid_structure(sep3dc)  || 0!=strcmp(sep3dc->first,"1MClean"))
   sep3d_clean(sep3dc);

if(0!=sep3d_grab_ndims(sep3dt,&sep3dc->ndims))
 return(sepwarn(NOT_MET,
   "trouble obtaining number of dimensions from %s\n",sep3dt));


if(0!=sep3d_axes_allocate(sep3dc,sep3dc->ndims))
 return(sepwarn(NOT_MET, "trouble allocataing axes\n"));

for(i=0; i < sep3dc->ndims; i++){
  if(0!=sep3d_grab_axis(sep3dt,i+1, &(sep3dc->n[i]), &(sep3dc->o[i]),
    &(sep3dc->d[i]), sep3dc->label[i], sep3dc->unit[i]))
     return(sepwarn(NOT_MET,"trouble grabbing axis %d from tag %s \n",
       i+1,sep3dt));
}


if(sep3dc->ndims>0){
  if(0!= sep3d_grab_wind(sep3dt,sep3dc->nwind,sep3dc->fwind,sep3dc->jwind))
     seperr("trouble grabbing window parameters");
}



if(0!=sep3d_grab_nkeys(sep3dt,&(sep3dc->nkeys)))
 return(sepwarn(NOT_MET,
   "trouble obtaining number of keys from %s\n",sep3dt));

if(sep3dc->nkeys >0){
  if(0!=sep3d_key_allocate(sep3dc,sep3dc->nkeys))
   return(sepwarn(NOT_MET, "trouble allocating keys\n"));

  for(i=0; i < sep3dc->nkeys; i++){
    if(0!=sep3d_grab_key(sep3dt,i+1, sep3dc->keyname[i], sep3dc->keytype[i],
      sep3dc->keyfmt[i]))
       return(sepwarn(NOT_MET,"trouble grabbing key %d from tag %s \n",
         i+1,sep3dt));
  }
}


if(0!=sep3d_grab_usage(sep3dt,sep3dc->usage) ||
  0!=sep3d_grab_file_type(sep3dt,sep3dc->file_format) ||
  0!=sep3d_grab_data_type(sep3dt,sep3dc->data_format) )
   return(sepwarn(NOT_MET,"Trouble grabbing usage,file_format,data_format \n"));

if(0!= sep3d_grab_long_ntraces(sep3dt,&sep3dc->ntraces))
  return(sepwarn(NOT_MET,"trouble grabing into C the number of traces\n"));

if(strcmp(sep3dc->tag,sep3dt)!=0) strncpy(sep3dc->tag,sep3dt,sizeof(sep3dc->tag));




return(SUCCESS);
}

/*
=head1 NAME

sep3d_set_sep3d - Store to the internal structure


=head1 SYNOPSIS

C<logic= sep3d_set_sep3d(sep3d *sep3dc)>

=head1 INPUT PARAMETERS

=over 4

=item sep3dc  -  sep3d

      structure


=back

=head1 DESCRIPTION

Store the public sep3dc contents into the  private structure

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3d_set_sep3d(sep3d *sep3dc){
char tag[256];
int i1;
if(1==valid_structure(sep3dc))
 return(sepwarn(NOT_MET,"structure not valid can't sync it \n"));

strncpy(tag,sep3dc->tag,sizeof(tag));


if(0!= sep3d_par_init(tag,sep3dc->usage))
  return(sepwarn(NOT_MET,"trouble initializing tag\n"));


if(sep3dc->ndims>0){
  if(0!=sep3d_set_ndims(tag,sep3dc->ndims))
    return(sepwarn(NOT_MET,"trouble setting ndims tag=%s \n",tag));

  for(i1=0; i1 < sep3dc->ndims; i1++){
   if(0!=sep3d_set_axis(tag,i1+1,sep3dc->n[i1],sep3dc->o[i1],sep3dc->d[i1],
      sep3dc->label[i1],sep3dc->unit[i1]))
    return(sepwarn(NOT_MET,"trouble setting axis %d for tag %s \n",
     i1+1,tag));
  }
  if(0!= sep3d_set_wind(sep3dc->tag,sep3dc->nwind,sep3dc->fwind,sep3dc->jwind))
     seperr("trouble grabbing window parameters");
}
if(sep3dc->nkeys>0){
  if(0!= sep3d_set_nkeys(tag,sep3dc->nkeys))
   return(sepwarn(NOT_MET,"trouble allocating nkeys for tag=%s \n",
    sep3dc->nkeys,tag));

  for(i1=0; i1 < sep3dc->nkeys; i1++){
    if(0!=sep3dc->keyname[i1],"data_record_number"){
     if(0!= sep3d_set_key(tag,i1+1,sep3dc->keyname[i1],sep3dc->keytype[i1],
       sep3dc->keyfmt[i1]))
       return(sepwarn(NOT_MET,"trouble setting key %d for tag %s \n",
         i1+1,tag));
    }
  }
}
if(0!= sep3d_set_file_type(tag,sep3dc->file_format)
   ||sep3d_set_data_type(tag,sep3dc->data_format))
    return(sepwarn(NOT_MET,"Trouble setting file and data type for tag %s \n",
    tag));
if(0!= sep3d_set_ntraces(tag,sep3dc->ntraces))
  return(sepwarn(NOT_MET,"trouble putting into C the number of traces\n"));


return(SUCCESS);
}

/*
=head1 NAME

init_sep3d_tag - Initialize a structure from a tag


=head1 SYNOPSIS

C<logic= init_sep3d_tag(char *tag,sep3d *sep3dc,char *usage)>

=head1 INPUT PARAMETERS

=over 4

=item tag  -  sepfile

     seplib tag to initialize from

=item sep3dc  -  sep3d

      structure

=item usage  -  char*

      usage for tag ('INPUT','OUTPUT','SCRATCH')


=back

=head1 DESCRIPTION

Create a sep3d structure from a tag

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

int init_sep3d_tag(char *tag,sep3d *sep3dc,char *usage){
/*
char  temp_ch[1024];
int local_tag,i;


sprintf(temp_ch,"%s.local_tag",tag);
local_tag=0;
getch(temp_ch,"d",&local_tag);
*/

sep3d_initialize(sep3dc);


/*
 if(local_tag==1){
*/
   if(0!= sep3d_tag_init(tag,tag,usage)) seperr("trouble initializing tag\n");
/*
 }else{
    if(0!= sep3d_tag_init_thread(tag,tag,usage,0)) 
      seperr("trouble initializing tag\n");
 }
*/
if(0!=sep3d_grab_sep3d(tag,sep3dc))
  return(sepwarn(NOT_MET,"Trouble grabbing sep3d structure \n"));


return(SUCCESS);
}
/*
=head1 NAME

init_sep3d_tag - Initialize a structure from another structure


=head1 SYNOPSIS

C<logic= init_sep3d_struct(sep3d sep3din,sep3d *sep3dc,char *usage)>

=head1 INPUT PARAMETERS

=over 4

=item sep3din  -  sep3d

     structure to initialize from

=item sep3dc  -  sep3d

      structure

=item usage  -  char*

      usage for tag ('INPUT','OUTPUT','SCRATCH')


=back

=head1 DESCRIPTION

Create a sep3d structure from another structure

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/
int init_sep3d_struct(sep3d sep3din,sep3d *sep3dout,char *usage){
char tag[256];

sep3d_initialize(sep3dout);

sprintf(tag,"sep3dc%d",sep3d_tag_num); sep3d_tag_num++;
if(0!=sep3d_struct_init(sep3din.tag,tag,usage))
  return(sepwarn(NOT_MET,"Trouble with sep3d_struct_init  \n"));
 

if(0!=sep3d_grab_sep3d(tag,sep3dout))
  return(sepwarn(NOT_MET,"Trouble grabbing sep3d structure \n"));

return(SUCCESS);
}
/*
=head1 NAME

init_sep3d_par - Initialize a structure from parameters


=head1 SYNOPSIS

C<logic= init_sep3d_par(sep3d *sep3dc,char *usage,char *data_type,char *file_type,int ndim,int nkeys)>

=head1 INPUT PARAMETERS

=over 4

=item sep3c  -  sep3d

     structure to create

=item usage  -  char*

      usage for tag ('INPUT','OUTPUT','SCRATCH')

=item data_type - char*

     datatype ("FLOAT","COMPLEX","INTEGER","BYTE")

=item file_type - char*

     filetype ("REGULAR","HEADERS","GRID")

=item ndim - int

     number of dimensions

=item nkeys - int

     number of keys

=back

=head1 DESCRIPTION

Create a sep3d structure from another structure

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/
int init_sep3d_par(sep3d *sep3dc,char *usage,char *data_type,char *file_type,int ndim,int nkeys){
char  temp_ch[1024];
sep3d_initialize(sep3dc);
sprintf(sep3dc->tag,"sep3dc%d",sep3d_tag_num); sep3d_tag_num++;



if(0!=sep3d_par_init(sep3dc->tag,usage))
  return(sepwarn(NOT_MET,"Trouble with sep3d_par_init  \n"));



if(0!=sep3d_set_ndims(sep3dc->tag,ndim))
 return(sepwarn(NOT_MET,"Trouble setting ndims \n"));

if(0!=sep3d_set_nkeys(sep3dc->tag,nkeys))
 return(sepwarn(NOT_MET,"Trouble setting nkeys \n"));

if(0!=sep3d_set_file_type(sep3dc->tag,file_type))
  return(sepwarn(NOT_MET,"Trouble setting file type \n"));


if(0!=sep3d_set_data_type(sep3dc->tag,data_type))
  return(sepwarn(NOT_MET,"Trouble setting data type \n"));

strncpy(temp_ch,sep3dc->tag,sizeof(sep3dc->tag));
if(0!=sep3d_grab_sep3d(temp_ch,sep3dc))
  return(sepwarn(NOT_MET,"trouble grabbing structure \n"));


return(SUCCESS);
}

int valid_structure(sep3d *sep3dc){
if( sep3dc->n==SEPNULL2 ||
 sep3dc->o==SEPNULL2 ||
 sep3dc->d==SEPNULL2 ||
 sep3dc->label==SEPNULL2 ||
 sep3dc->unit==SEPNULL2 )
  return(1);


return(SUCCESS);
}

/*
=head1 NAME

sep3d_initialize - Initialize a structure


=head1 SYNOPSIS

C<logic= sep3d_initialize(sep3d *sep3dc)

=head1 INPUT PARAMETERS

=over 4

=item sep3c  -  sep3d

     structure 

=back

=head1 DESCRIPTION

Initialize a structure

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

void sep3d_initialize(sep3d *sep3dc){

strcpy(sep3dc->tag,"");
sep3dc->n=SEPNULL2;
sep3dc->o=SEPNULL2;
sep3dc->d=SEPNULL2;
sep3dc->fwind=SEPNULL2;
sep3dc->jwind=SEPNULL2;
sep3dc->nwind=SEPNULL2;
sep3dc->label=sep3dc->unit=SEPNULL2;
sep3dc->ndims=0;
sep3dc->keyname=sep3dc->keytype=sep3dc->keyfmt=SEPNULL2;
sep3dc->nkeys=0;
sep3dc->drn=-1;
sep3dc->ntraces=0;
strcpy(sep3dc->first,"1MClean");
}

/*
=head1 NAME

sep3d_key_allocate - Allocate keys for a structure


=head1 SYNOPSIS

C<logic= sep3d_key_allocate(sep3d *sep3dc,int nkeys)

=head1 INPUT PARAMETERS

=over 4

=item sep3c  -  sep3d

     structure 

=item nkeys  -  int

     number of keys 

=back

=head1 DESCRIPTION

Allocate keys in a structure

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

int sep3d_key_allocate(sep3d *sep3dc,int nkey){
int i;

my_keys_clean(sep3dc);
sep3dc->nkeys=nkey;
sep3dc->keyname=(char**) malloc(sizeof(char*)*nkey);
sep3dc->keyfmt=(char**) malloc(sizeof(char*)*nkey);
sep3dc->keytype=(char**) malloc(sizeof(char*)*nkey);
for(i=0; i < nkey; i++){
  sep3dc->keyname[i]=(char*) malloc(sizeof(char)*256);
  strcpy(sep3dc->keyname[i],"UNSPECIFIED");
  sep3dc->keytype[i]=(char*) malloc(sizeof(char)*256);
  strcpy(sep3dc->keytype[i],"UNSPECIFIED");
  sep3dc->keyfmt[i]=(char*) malloc(sizeof(char)*256);
  strcpy(sep3dc->keyfmt[i],"UNSPECIFIED");
}


return(SUCCESS);
}
/*
=head1 NAME

sep3d_axes_allocate - Allocate axes for a structure


=head1 SYNOPSIS

C<logic= sep3d_key_allocate(sep3d *sep3dc,int naxes)

=head1 INPUT PARAMETERS

=over 4

=item sep3c  -  sep3d

     structure 

=item naxes  -  int

     number of axes 

=back

=head1 DESCRIPTION

Allocate axes in a structure

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3d_axes_allocate(sep3d *sep3dc,int ndim){
int i;

my_axes_clean(sep3dc);
sep3dc->ndims=ndim;
sep3dc->nwind=(int*)malloc(sizeof(int)*ndim);
sep3dc->fwind=(int*)malloc(sizeof(int)*ndim);
sep3dc->jwind=(int*)malloc(sizeof(int)*ndim);
sep3dc->n=(int*)malloc(sizeof(int)*ndim);
sep3dc->o=(float*)malloc(sizeof(float)*ndim);
sep3dc->d=(float*)malloc(sizeof(float)*ndim);
sep3dc->label=(char**) malloc(sizeof(char*)*ndim);
sep3dc->unit=(char**) malloc(sizeof(char*)*ndim);
for(i=0; i < ndim; i++){
  sep3dc->n[i]=1;
  sep3dc->nwind[i]=1;
  sep3dc->fwind[i]=0;
  sep3dc->jwind[i]=1;
  sep3dc->o[i]=0.;
  sep3dc->d[i]=1.;
  sep3dc->label[i]=(char*) malloc(sizeof(char)*256);
  strcpy(sep3dc->label[i],"UNSPECIFIED");
  sep3dc->unit[i]=(char*) malloc(sizeof(char)*256);
  strcpy(sep3dc->unit[i],"UNSPECIFIED");
}



return(SUCCESS);
}


void my_axes_clean(sep3d *sep3dc){

if(SEPNULL2!=sep3dc->nwind) free(sep3dc->nwind);
if(SEPNULL2!=sep3dc->fwind) free(sep3dc->fwind);
if(SEPNULL2!=sep3dc->jwind) free(sep3dc->jwind);
if(SEPNULL2!=sep3dc->n) free(sep3dc->n);
if(SEPNULL2!=sep3dc->o){ free(sep3dc->o); }
if(SEPNULL2!=sep3dc->d) free(sep3dc->d);
if(SEPNULL2!=sep3dc->label){
  my_char_dealloc(sep3dc->label,sep3dc->ndims); free(sep3dc->label);
}
if(SEPNULL2!=sep3dc->unit){
   my_char_dealloc(sep3dc->unit,sep3dc->ndims); free(sep3dc->unit);
}
sep3dc->n=SEPNULL2;
sep3dc->o=SEPNULL2;
sep3dc->d=SEPNULL2;
sep3dc->label=SEPNULL2;
sep3dc->unit=SEPNULL2;
}

void my_keys_clean(sep3d *sep3dc){

if(SEPNULL2!=sep3dc->keyname){
   my_char_dealloc(sep3dc->keyname,sep3dc->nkeys); free(sep3dc->keyname);
}
if(SEPNULL2!=sep3dc->keytype){
   my_char_dealloc(sep3dc->keytype,sep3dc->nkeys); free(sep3dc->keytype);
}
if(SEPNULL2!=sep3dc->keyfmt){
   my_char_dealloc(sep3dc->keyfmt,sep3dc->nkeys); free(sep3dc->keyfmt);
}
sep3dc->keyname=SEPNULL2;
sep3dc->keytype=SEPNULL2;
sep3dc->keyfmt=SEPNULL2;
}
/*
=head1 NAME

sep3dc_delete - Delete a structure


=head1 SYNOPSIS

C<logic= sep3dc_delete(sep3d *sep3dc)

=head1 INPUT PARAMETERS

=over 4

=item sep3c  -  sep3d

     structure 

=back

=head1 DESCRIPTION

Delete  a structure

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3dc_delete(sep3d *sep3dc){

if(SUCCESS!=sep3d_delete(sep3dc->tag))
  return(sepwarn(NOT_MET,"trouble deleting structure \n"));

sep3d_clean(sep3dc);

return(SUCCESS);
}
/*
=head1 NAME

sep3dc_clean - Clean a structure


=head1 SYNOPSIS

C<logic= sep3dc_clean(sep3d *sep3dc)

=head1 INPUT PARAMETERS

=over 4

=item sep3c  -  sep3d

     structure 

=back

=head1 DESCRIPTION

Clean a structure

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/



void sep3d_clean(sep3d *sep3dc){
int i;
if(0==strcmp(sep3dc->first,"1MClean")){

my_axes_clean(sep3dc);
my_keys_clean(sep3dc);


}
sep3d_initialize(sep3dc);

}

void my_char_dealloc(char **str_array, int nlen){
int i,ierr;


for(i=0 ; i < nlen;i++) {
  if(SEPNULL2!=str_array[i]) free(str_array[i]);
}

}
/*
=head1 NAME

sep3dc_grab_headers -  Grab a window of the headers


=head1 SYNOPSIS

C<logic=sep3dc_grab_headers(char *tag, sep3d *sep3dc, int *nh, int *nwind, int *fwind, int *jwind)


=head1 INPUT PARAMETERS

=over 4

=item tag  -  sepfile

     tag to grab the headers from 

=item sep3c  -  sep3d

     structure 

=item nh  -  int*

     snumber of header

=item nwind  -  int*

     Number of elements along each axis (starting with axis 2) 

=item fwind  -  int*

     First elements along each axis (starting with axis 2) 

=item jwind  -  int*

    Sampling elements along each axis (starting with axis 2) 

=back

=head1 DESCRIPTION

Grab a window the headers 

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3dc_grab_headers(char *tag, sep3d *sep3dc, int *nh, int *nwind, int
*fwind, int *jwind){

if(1==valid_structure(sep3dc))
 return(sepwarn(NOT_MET,"structure not valid can't grab headers \n"));



   if(0!=sep3d_read_headers(tag,sep3dc->tag,nwind,fwind,jwind,nh))
     seperr("trouble reading headers");
                                                                                
return(0);
}




/*
=head1 NAME

sep3dc_copy -  Copy the contents of one structure to another


=head1 SYNOPSIS

C<logic=sep3dc_copy(sep3d *in,sep3d *out)


=head1 INPUT PARAMETERS

=over 4

=item  in  -  sep3d

     structure 

=item  out  -  sep3d

     structure 

=head1 DESCRIPTION

 Copy the contents of one structure to another

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/
int sep3dc_copy(sep3d *input, sep3d *output){
int ierr,nh,i;
int *trnum;
int esize,ntr;
char  temp_ch[1024];

if(1==valid_structure(input))
 return(sepwarn(NOT_MET,"structure not valid can't grab data \n"));

ierr=sep3d_copy_struct(input->tag,output->tag);

                                                                                
strncpy(temp_ch,output->tag,sizeof(output->tag));
 if(0!=sep3d_grab_sep3d(temp_ch,output))
   return(sepwarn(NOT_MET,"Trouble grabbing structure \n"));


return(SUCCESS);


}

/*
=head1 NAME

sep3dc_read_data -  Read the data assoicated with a given window



=head1 SYNOPSIS

C<logic= sep3dc_read_data(char *tag, sep3d *sep3dc, void *data, int nt, int ft, int jt)>



=head1 INPUT PARAMETERS

=over 4

=item  tag  -  sepfile

     file to read  the data from

=item  sep3dc  -  sep3d

     structure 

=item  data  -  void*

     data 

=item  nt  -  int

     number of samples along axis 1

=item  ft  -  int

     first sample along axis 1

=item  jt  -  int

     sampling along axis 1

=head1 DESCRIPTION


Read in  a portion of the dataset

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

int sep3dc_read_data(char *tag, sep3d *sep3dc, char *data,
  int nt, int ft, int jt){
int *trnum;
int esize,ntr,i;

if(1==valid_structure(sep3dc))
 return(sepwarn(NOT_MET,"structure not valid can't grab data \n"));

esize=sep3dc_get_esize(sep3dc);


     if(0!=sep3d_read_data(tag,sep3dc->tag,nt,ft,jt,data))
       seperr("trouble reading data\n");
                                                                                
return (0);

}
/*
=head1 NAME

sep3dc_write_description -  Write a description of the dataset



=head1 SYNOPSIS

C<logic= sep3dc_write_description(char *tag, sep3d *sep3dc)>



=head1 INPUT PARAMETERS

=over 4

=item  tag  -  sepfile

     file to read  the data from

=item  sep3dc  -  sep3d

     structure 

=back

=head1 DESCRIPTION


Write a description of a dataset to  a file

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

int  sep3dc_write_description(char *tag, sep3d *sep3dc){

if(0!= valid_structure(sep3dc))
  return(sepwarn(NOT_MET,
  "you are attempting to write out an invalid strucutre\n"));

  if(0!= sep3d_set_sep3d(sep3dc))
    return(sepwarn(NOT_MET, "trouble setting sep3dc (1)\n"));


if(0!=sep3d_set_ntraces(sep3dc->tag,sep3dc->ntraces))
  return(sepwarn(NOT_MET, "trouble writing out structure to disk1\n"));
if(0!=sep3d_rite_format(tag,sep3dc->tag))
  return(sepwarn(NOT_MET, "trouble writing out structure to disk2\n"));
if(0!=sep3d_rite_ntraces(tag,sep3dc->tag))
  return(sepwarn(NOT_MET, "trouble writing out structure to disk3\n"));

return(0);
}
/*
=head1 NAME

sep3dc_write_status - Set write status



=head1 SYNOPSIS

C<logic= sep3dc_write_status(char *sep3dc, int data, int headers)>



=head1 INPUT PARAMETERS

=over 4

=item  sep3dc  -  sep3d

     structure 

=item  data  -  int

    whether [1] or not [0] we have written out data

=item  headers  -  int

    whether [1] or not [0] we have written out headers

=back

=head1 DESCRIPTION


Set whether the data and/or headers of this dataset have been written (necessary

to know where to write the number of traces in a dataset)

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

int sep3dc_set_write_status(sep3d *sep3dc,int data,int header){

if(0!=valid_structure(sep3dc))
  return(sepwarn(NOT_MET,
  "you are attempting to write out an invalid strucutre\n"));

  if(0!= sep3d_set_sep3d(sep3dc))
    return(sepwarn(NOT_MET,
     "trouble setting sep3dc structure \n"));
 

  if(0!=sep3d_set_rite_status(sep3dc->tag,data,header))
    return(sepwarn(NOT_MET,
     "trouble writing out structure to disk \n"));

  return(SUCCESS);
}
/*
=head1 NAME

sep3dc_rite_num_traces - write the number of traces



=head1 SYNOPSIS

C<logic= sep3dc_rite_num_traces(char *tag, sep3dc *struct)>



=head1 INPUT PARAMETERS

=over 4

=item  tag  -  sepfile

     tag

=item  sep3dc  -  sep3d

     structure 

=back

=head1 DESCRIPTION


Write the number of traces to a file

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

int sep3d_rite_num_traces(char *tag,sep3d *sep3dc){

if(0!=sep3d_rite_ntraces(tag,sep3dc->tag)){
  return(sepwarn(NOT_MET,
   "trouble writing number of traces for %s\n",tag));
}

return(SUCCESS);
}
int sep3dc_write_data(char *tag,sep3d *sep3dc,char *data,int *nwind,int *fwind,
int *jwind,int nh,int write_headers,int write_grid){
int esize,write_data,i;

if(0!=valid_structure(sep3dc))
  return(sepwarn(NOT_MET,
  "you are attempting to write out an invalid strucutre\n"));

if(0==strcmp(sep3dc->data_format,"COMPLEX")) esize=8;
else esize=4;

if(SEPNULL2 != data) write_data=1;
else write_data=0;



      if(0!=sep3d_rite(tag,sep3dc->tag, nwind,fwind,jwind,data, nh,write_data,write_headers,write_grid))
        seperr("trouble writing out the data sep3dc_write_data\n");


return(SUCCESS);
}

int print_sep3dc(sep3d *sep3dc){

sep3d_print_info(sep3dc->tag);

return(SUCCESS);
}
/*
=head1 NAME

sep3dc_ndims - the last axis of length >1



=head1 SYNOPSIS

C<ndims= sep3dc_ndims(sep3dc *struct)>



=head1 INPUT PARAMETERS

=over 4

=item  sep3dc  -  sep3d

     structure 

=back

=head1 DESCRIPTION


The last axis of length >1

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

int sep3dc_ndims(sep3d *sep3dc){
int ndim;

if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

ndim=sep3dc->ndims;
while(ndim!=1&& sep3dc->n[ndim-1] ==1) ndim--;
return(ndim);
}
/*
=head1 NAME

sep3dc_ge_space - whether space1 contains space 2



=head1 SYNOPSIS

C<ndims= sep3dc_ge_space(sep3dc *space1, sep3dc *space2)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  space2  -  sep3d

     structure 

=back

=head1 DESCRIPTION


Whether space1 contains space2

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3dc_ge_space(sep3d *s1, sep3d *s2){
if(0!=valid_structure(s1))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
if(0!=valid_structure(s2))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
  if(0!= sep3d_set_sep3d(s1))
    return(sepwarn(NOT_MET,
     "trouble setting sep3dc structure \n"));

  if(0!= sep3d_set_sep3d(s2))
    return(sepwarn(NOT_MET,
     "trouble setting sep3dc structure \n"));

return(sep3d_ge_space(s1->tag,s2->tag));

}
/*
=head1 NAME

sep3dc_conform - whether space1 and  space2 conform



=head1 SYNOPSIS

C<ndims= sep3dc_conform(sep3dc *space1, sep3dc *space2)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  space2  -  sep3d

     structure 

=back

=head1 DESCRIPTION


Whether space1 and space2 have the same n, o, and d

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3dc_conform(sep3d *s1, sep3d *s2){
  if(0!= sep3d_set_sep3d(s1))
    return(sepwarn(NOT_MET,
     "trouble setting sep3dc structure \n"));

  if(0!= sep3d_set_sep3d(s2))
    return(sepwarn(NOT_MET,
     "trouble setting sep3dc structure \n"));

if(0!=valid_structure(s1))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
if(0!=valid_structure(s2))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
return(sep3d_conform(s1->tag,s2->tag));

}
/*
=head1 NAME

sep3dc_key_index - The index of a given key



=head1 SYNOPSIS

C<index= sep3dc_key_index(sep3dc *space1, char *key)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  key  -  char*

     keyname 

=back

=head1 DESCRIPTION


The internal key index associated with a given key (-1 if not found)

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3dc_key_index(sep3d *sep3dc, char *name){
int ierr,index;
if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
  if(0!= sep3d_set_sep3d(sep3dc))
    return(sepwarn(NOT_MET,
     "trouble setting sep3dc structure \n"));


  ierr=sep3d_grab_key_index(sep3dc->tag,name,&index);
 if(ierr!=0) return(-1);
  else return(index);
}

/*
=head1 NAME

sep3dc_key_index - The index of a given axis



=head1 SYNOPSIS

C<index= sep3dc_key_index(sep3dc *space1, char *axis_name)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  axis_name  -  char*

     The label of an axis

=back

=head1 DESCRIPTION


The internal axis index associated with a given key (-1 if not found)

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

int sep3dc_axis_index(sep3d *sep3dc, char *name){
int ierr,index;
if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

  if(0!= sep3d_set_sep3d(sep3dc))
    return(sepwarn(NOT_MET,
     "trouble setting sep3dc structure \n"));


  ierr=sep3d_grab_axis_index(sep3dc->tag,name,&index);
 if(ierr!=0) return(-1);
  else return(index);
}
/*
=head1 NAME

sep3dc_grab_key_vals - Grab header values based on keyname



=head1 SYNOPSIS

C<logic= sep3dc_grab_key_vals(sep3dc *space1, char *key, float *headers)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  key  -  char*

     The keyname

=item  headers  -  float*

     Header values

=back

=head1 DESCRIPTION


Grab header values based on the keyname

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

int sep3dc_grab_key_vals(sep3d *sep3dc, char *name, float *header){

if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

  return(sep3d_grab_header_vals_s(sep3dc->tag,name,(int*)header));
}
/*
=head1 NAME

sep3dc_set_key_vals - Set header values based on keyname



=head1 SYNOPSIS

C<logic= sep3dc_set_key_vals(sep3dc *space1, char *key, float *headers)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  key  -  char*

     The keyname

=item  headers  -  float*

     Header values

=back

=head1 DESCRIPTION


Set header values based on the keyname

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3dc_set_key_vals(sep3d *sep3dc, char *name, float *header){

if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

  return( sep3d_set_header_vals_s(sep3dc->tag,name,(int*)header));
}


/*
=head1 NAME

sep3dc_grab_key_vali - Grab header values based on index



=head1 SYNOPSIS

C<logic= sep3dc_grab_key_vali(sep3dc *space1, int key, float *headers)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  key  -  int

     The keyindex

=item  headers  -  float*

     Header values

=back

=head1 DESCRIPTION


Grab header values based on the keyindex

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/



int sep3dc_grab_key_vali(sep3d *sep3dc,int num, float *header){

if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

  return( sep3d_grab_header_vals_i(sep3dc->tag,num,(int*)header));
}

/*
=head1 NAME

sep3dc_set_key_vals - Set header values based on keyinex



=head1 SYNOPSIS

C<logic= sep3dc_set_index(sep3dc *space1, int key, float *headers)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  key  -  int

     The keyindex

=item  headers  -  float*

     Header values

=back

=head1 DESCRIPTION


Set header values based on the keyindex

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3dc_set_key_vali(sep3d *sep3dc, int num, float *header){

if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

  return( sep3d_set_header_vals_i(sep3dc->tag,num,(int*)header));
}
/*
=head1 NAME

sep3dc_grab_grid_values - Grab grid values



=head1 SYNOPSIS

C<logic= sep3dc_grab_grid_values(sep3dc *space1, int *grid)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  grid  -  int*

     Grid values

=back

=head1 DESCRIPTION


Grab grid values of the window currently in memory

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

int sep3dc_grab_grid_values(sep3d *sep3dc, int *vals){

if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
return(sep3d_grab_grid_block(sep3dc->tag,vals));
}
/*
=head1 NAME

sep3dc_set_grid_values - Set grid values



=head1 SYNOPSIS

C<logic= sep3dc_grab_grid_values(sep3dc *space1, int *grid)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  grid  -  int*

     Grid values

=back

=head1 DESCRIPTION


Set grid values of the window currently in memory

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3dc_set_grid_values(sep3d *sep3dc, int *vals){

if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
return(sep3d_set_grid_vals(sep3dc->tag,vals));
}
/*
=head1 NAME

sep3dc_coord_copy - Copy the coordinates from one structure to another



=head1 SYNOPSIS

C<logic= sep3dc_coord_copy(sep3dc *space1, sep3dc *space2)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  space2  -  sep3d

     structure 

=back

=head1 DESCRIPTION


Copy the coordinate values from struct1 to struct2

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3dc_coord_copy(sep3d *s1, sep3d *s2){
if(0!=valid_structure(s2))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
if(0!=valid_structure(s1))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
    return(sep3d_copy_coords(s1->tag,s2->tag));
}

int sep3dc_set_ncoord(sep3d *s, int ncoord){

 return(sep3d_alloc_coord(s->tag,ncoord));
}

int sep3dc_grab_coord_vals(sep3d *s1, int *coords)
{
return(sep3d_grab_coord_vals(s1->tag,coords));
}

int sep3dc_grab_coordh(sep3d *s1, long long *coords)
{
return(sep3d_grab_coordh(s1->tag,coords));
}

int sep3dc_set_coord_vals(sep3d *s1, int *coords)
{
return(sep3d_set_coord_vals(s1->tag,coords));
}

int sep3dc_set_coordh(sep3d *s1, long long *coords)
{
return(sep3d_set_coordh(s1->tag,coords));
}

/*
=head1 NAME

sep3dc_header_copy - Copy the headers from one structure to another



=head1 SYNOPSIS

C<logic= sep3dc_header_copy(sep3dc *space1, sep3dc *space2)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  space2  -  sep3d

     structure 

=back

=head1 DESCRIPTION


Copy the headers values from struct1 to struct2

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/

int sep3dc_header_copy(sep3d *s1, sep3d *s2){
if(0!=valid_structure(s2))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
if(0!=valid_structure(s1))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
    return(sep3d_copy_headers(s1->tag,s2->tag));
}

int sep3dc_set_drn(sep3d *s1, int *drn){
if(0!=valid_structure(s1))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
    return(sep3d_set_drn(s1->tag,drn));
}
int sep3dc_grab_drn(sep3d *s1, int *drn){
if(0!=valid_structure(s1))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
    return(sep3d_grab_drn(s1->tag,drn));
}

int sep3dc_drn_copy(sep3d *s1, sep3d *s2){
if(0!=valid_structure(s2))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
if(0!=valid_structure(s1))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
    return(sep3d_drn_copy(s1->tag,s2->tag));
}



/*
=head1 NAME

sep3dc_grid_copy - Copy the grid from one structure to another



=head1 SYNOPSIS

C<logic= sep3dc_grid_copy(sep3dc *space1, sep3dc *space2)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  space2  -  sep3d

     structure 

=back

=head1 DESCRIPTION


Copy the grid values from struct1 to struct2

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3dc_grid_copy(sep3d *s1, sep3d *s2){
if(0!=valid_structure(s2))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
if(0!=valid_structure(s1))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
    return(sep3d_copy_grid(s1->tag,s2->tag));
}


/*
=head1 NAME

sep3dc_set_number_headers - Set the number of headers in memory



=head1 SYNOPSIS

C<logic= sep3dc_set_number_headers(sep3dc *space1, int num)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  num  -  int

     number of headers

=back

=head1 DESCRIPTION

Set the number of headers

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/


int sep3dc_set_number_headers(sep3d *sep3dc, int num){

if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

  return(sep3d_set_nh(sep3dc->tag,num));
}

/*
=head1 NAME

sep3dc_update_ntraces - Update number of traces



=head1 SYNOPSIS

C<logic= sep3dc_update_ntraces(sep3dc *space1)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=back

=head1 DESCRIPTION

Copy the number of traces written to the number of traces in the dataset

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/



int sep3dc_update_ntraces(sep3d *sep3dc){

if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

  return(sep3d_count_ntraces(sep3dc->tag));
}

/*
=head1 NAME

sep3dc_inorder -  Signify that the data and headers are synched



=head1 SYNOPSIS

C<logic= sep3dc_inorder(sp3dc *space1)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=back

=head1 DESCRIPTION

Signify that the headers and data are synched

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/



int sep3dc_grab_inorder(sep3d *sep3dc,int *inorder){

if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));


  if(0!=sep3d_grab_inorder(sep3dc->tag,inorder))
    return(-1);
  return(0);
}


int sep3dc_inorder(sep3d *sep3dc){

if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

  return(sep3d_set_inorder(sep3dc->tag));
}


/*
=head1 NAME

sep3dc_reshape - Reshape a dataset



=head1 SYNOPSIS

C<logic= sep3dc_reshape(sep3dc *space1, int *n)>



=head1 INPUT PARAMETERS

=over 4

=item  space1  -  sep3d

     structure 

=item  n  -  int*

     Reshape a dataset

=back

=head1 DESCRIPTION

Reshape a dataset. See the sep3d_reshape manpage for examples

=head1 SEE ALSO


=head1 LIBRARY

B<superset>

=cut
*/



int  sep3dc_reshape(sep3d *sep3dc,int ndim, int *n){
char  temp_ch[1024];

if(0!=valid_structure(sep3dc))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

  if(0!= sep3d_set_sep3d(sep3dc))
    return(sepwarn(NOT_MET,
     "trouble setting sep3dc structure \n"));

if(0!=sep3d_change_dims(sep3dc->tag,ndim,n)){
  return(sepwarn(NOT_MET,
  "trouble changing dimensions for tag  %s \n",sep3dc->tag));
}
strncpy(temp_ch,sep3dc->tag,sizeof(sep3dc->tag));
if(0!=sep3d_grab_sep3d(temp_ch,sep3dc))
   return(sepwarn(NOT_MET,"Trouble grabbing structure \n"));
return(0);
}

                                                                                

int  sep3dc_close_tag(char *tag,sep3d *sep3din){
if(0!=valid_structure(sep3din))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

 if(0!=sep3d_close_tag(tag,sep3din->tag))
  return(sepwarn(-1, "trouble closeing tag \n"));


  return(SUCCESS);
}

int  sep3dc_get_esize(sep3d *sep3din){
int esize;
if(0!=valid_structure(sep3din))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

if(0==strcmp(sep3din->data_format,"FLOAT")) esize=4;
else if(0==strcmp(sep3din->data_format,"COMPLEX")) esize=8;
else if(0==strcmp(sep3din->data_format,"BYTE")) esize=1;
else if(0==strcmp(sep3din->data_format,"INTEGER")) esize=4;

return(esize);
}

                                                                                


int  sep3dc_grab_header_block(sep3d *sep3din,void *block){
if(0!=valid_structure(sep3din))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

 if(0!=sep3d_grab_header_block(sep3din->tag,block))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
  return(SUCCESS);
}

int  sep3dc_set_header_block(sep3d *sep3din,void *block){
if(0!=valid_structure(sep3din))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));

 if(0!=sep3d_set_header_block(sep3din->tag,block))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
  return(SUCCESS);
}

                                                                                
int  sep3dc_rite_file_stat(sep3d *sep3din,int idat,int ihead){
if(0!=valid_structure(sep3din))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));
 if(0!=sep3d_rite_file_stat(sep3din->tag,idat,ihead))
  return(sepwarn(-1,
  "you are attempting to get info on an invalid structure \n"));


  return(SUCCESS);
}
int sep3dc_broadcast_data(sep3d *sep3dc, int ntr, int from, char *data){

if(1==valid_structure(sep3dc))
 return(sepwarn(NOT_MET,"structure not valid can't broadcast headers \n"));

if(0!=sep3d_broadcast_data(from, sep3dc_get_esize(sep3dc),
   sep3dc->n[0]*ntr,data))
      seperr("trouble distributing  headers\n");
                                                                                
return(0);
}

int sep3dc_pass_headers(sep3d *sep3dc, int ifrom,int ito){
int nlen;
char string[512];


if(ifrom==sep_thread_num()){
  if(1==valid_structure(sep3dc))
   return(sepwarn(NOT_MET,"structure not valid can't broadcast headers \n"));
}


if(0!= sep3d_pass_string(sep3dc->tag,ifrom,ito))
  return(sepwarn(NOT_MET,"trouble broadcasting string"));

if(0!=sep3d_pass_headers(sep_thread_num(),sep3dc->tag,ifrom,ito))
      seperr("trouble distributing  headers\n");

strncpy(string,sep3dc->tag,sizeof(sep3dc->tag));
if(0!=sep3d_grab_sep3d(string,sep3dc))
  return(sepwarn(NOT_MET,"trouble grabbing sep3d \n"));
                                                                                
return(0);
}
int sep3dc_broadcast_headers(sep3d *sep3dc, int ifrom){
char string[256];

if(ifrom==sep_thread_num()){
  if(1==valid_structure(sep3dc))
   return(sepwarn(NOT_MET,"structure not valid can't broadcast headers \n"));
}
if(0!= sep3d_broadcast_string(sep3dc->tag,ifrom))
  return(sepwarn(NOT_MET,"trouble broadcasting string"));
if(0!=sep3d_broadcast_headers(sep3dc->tag,0))
      seperr("trouble distributing  headers\n");
strncpy(string,sep3dc->tag,sizeof(sep3dc->tag));
sep3d_grab_sep3d(string,sep3dc);

                                                                                

return(0);
}


int sep3dc_grab_nh(sep3d *sep3dc, int *nh){

if(1==valid_structure(sep3dc))
 return(sepwarn(NOT_MET,"structure not valid can't broadcast headers \n"));

if(0!=sep3d_grab_nh(sep3dc->tag,nh))
      seperr("trouble grabbing nh\n");
                                                                                
return(0);
}

int sep3dc_clear(sep3d *sep3dc){
  int i;
  for(i=0; i< sep3dc->ndims; i++){
    sep3dc->nwind[i]=-1;
    sep3dc->fwind[i]=-1;
    sep3dc->jwind[i]=-1;
  }
  if(0!=sep3d_set_nh(sep3dc->tag,0))
    return(sepwarn(NOT_MET,"trouble setting nh=0"));
  if(0!=sep3d_alloc_coord(sep3dc->tag,0))
    return(sepwarn(NOT_MET,"trouble setting ncoord=0"));
  return(0);
}

