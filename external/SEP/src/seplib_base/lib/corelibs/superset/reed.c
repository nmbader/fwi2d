#define SET_SDOC 1
#define NO_DRN -100
/*#define NBUF 5000000*/
#define YES 1
#define NO 0

#include<string.h>
#include "superset_internal.h"
#include <superset.h> 
int wind_to_helix(int ndim,int *ngrid,int *nwind, int *fwind, int *jwind, long long *index);

extern int sep_thread_num(void);


/*
sep3d_read_headers
<
USAGE
ierr=sep3d_read_headers(char *tag, char *sep3dname, int *nwind, int *fwind, int *jwind,int *nh)


INPUT PARAMETERS
tag          -   char*    tag to read
sep3dname    -   char*    pointer to the sep_3d structured
nwind        -   int*     number of elements to copy
fwind        -   int*     first elements to copy    
jwind        -   int*     skip between elements to copy

OUTPUT PARAMETERS
nh           -   int*      number of headers in this section

RETURN VALUES
0    =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Reads headers into sep_3d structure from a tag

>
*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_read_headers(char *tag, char *sep3dname, int *nwind, int *fwind, int *jwind,int *nh)
_XFUNCPROTOEND
#else
int sep3d_read_headers(tag,sep3dname,nwind,fwind,jwind, nh)
char *tag,*sep3dname; 
int *nwind,*fwind,*jwind,*nh;
#endif 
{ 
int *header_list;
int i1,ntr,esize,*idim,n2,i2,temp2i,nkeys,*grid;
int *n_wind,*f_wind,*j_wind;
sep_3d *info;
int *temp_headers,drn,c,d,count,i,same,wind_elem;
int ierr,sz,last,nkeys_read,ntot;
char temp_ch[1024],temp2_ch[8192];
long long *coords;



info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));


for(i=0,ntot=1; i < info->ndims-1; i++){
  if(info->n[i+1]>1){
    info->nwind[i+1]=nwind[i]; 
    info->fwind[i+1]=fwind[i];
    info->jwind[i+1]=jwind[i];
  }
  else{
    info->nwind[i+1]=1;
    info->fwind[i+1]=0;
    info->jwind[i+1]=1;
  }
  ntot=ntot*info->nwind[i+1];
}
if(info->file_format==REGULAR){
   n_wind=(int*)malloc(sizeof(int)*info->ndims);
   f_wind=(int*)malloc(sizeof(int)*info->ndims);
   j_wind=(int*)malloc(sizeof(int)*info->ndims);
   f_wind[0]=0; j_wind[0]=n_wind[0]=1;
   for(i=1; i < info->ndims; i++){
     n_wind[i]=info->nwind[i];
     f_wind[i]=info->fwind[i];
     j_wind[i]=info->jwind[i];
   }
    
   if(0!= sep3d_set_wind(sep3dname,n_wind,f_wind,j_wind)){
     return(sepwarn(NOT_MET,"trouble setting sep window\n"));
   }
   if(0!=SEP3D_wind_coords(info,nh)){
     return(sepwarn(NOT_MET,"trouble setting up  coordinates for regular dataset \n"));
   }
   free(n_wind); free(f_wind); free(j_wind);
   return(SUCCESS); 
}



/*first check sanity */
if(info->file_format==UNSPECIFIED) 
	return(sepwarn(INVALID_DATA,"can not read undefined dataype (%s) \n",
   	sep3dname));
if(info->file_format==REGULAR)return(sepwarn(INVALID_DATA,
   "can not get headers from a regular dataset  (%s) \n",sep3dname));
if(info->nkeys==0) return(sepwarn(NOT_MET,
  "(reed_headers_from_grid) nkeys not set before attempting read (%s) \n",
 sep3dname));
      


if(info->file_format==GRID){ /*if we dealing with the grid file */
	sz=1;
	for(i=0; i< info->ndims-1;i++) sz=sz*info->nwind[i+1];
  grid=(int*)malloc(sizeof(int)*sz);
	
	ierr=sep_get_grid_window(tag,&(info->ndims),&(info->n[1]),
		&info->nwind[1],&info->fwind[1],&info->jwind[1],grid);
   if(ierr!=SUCCESS){
      sprintf(temp2_ch,"ERROR=%d\n",ierr);
      for(i=1; i < info->ndims; i++){
        sprintf(temp2_ch,"%sidim=%d ndim=%d nwind=%d fwind=%d jwind=%d \n",
            temp2_ch,i,info->n[i],info->nwind[i],info->fwind[i],info->jwind[i]);
      }
		return(sepwarn(FAIL_OTHER, 
			"trouble reading grid window from tag %s:\n%s \n",tag,temp2_ch));
  }

   if(0!=sep3d_set_grid_vals(sep3dname,grid)){
    free(grid);
		return(sepwarn(FAIL_OTHER, 
			"trouble setting grid window from tag %s \n",sep3dname));
   }
   ntr=info->ncoord;
   if(ntr==0) {
     free(grid); *nh=0; 
     if(0!=sep3d_clear_headers(info->name))
       seperr("trouble deleting headers \n");
     return(SUCCESS);
	}

 
		header_list=(int*) alloc(sizeof(int)*ntr);
		count=0;
		for(i1=0; i1 <  sz; i1++){;
    	if(grid[i1]>0){
       	header_list[count]=grid[i1];
       	count++;
    	}
		}
    free(grid);

} else{  /*we  don't have a grid so calculate header_list from window params*/



	idim=(int *) alloc(sizeof(int) * (info->ndims-1));
	wind_elem=1;
	for(i1=0; i1 < info->ndims-1; i1++){ 
          wind_elem=wind_elem*info->nwind[i1+1];
}
	header_list=(int *) malloc(sizeof(int) * wind_elem);
	ntr=wind_elem;
  if(0!=SEP3D_alloc_coord(info,ntr)){
      return(sepwarn(FAIL_OTHER,"trouble allocating coords \n"));
  }
	for(i2=0; i2< wind_elem; i2++){
                /* convert to helical coordinate sys*/
		h2c(i2,&info->nwind[1],info->ndims-1,idim); 
		for(i1=0; i1< info->ndims-1; i1++)
                  idim[i1]=info->fwind[i1+1]+info->jwind[i1+1]*idim[i1];
		c2h(&temp2i,&(info->n[1]),info->ndims-1,idim); /*convert back to reg*/
		header_list[i2]=temp2i+1;
    info->coord[i2]=temp2i;
	}
	free(idim);
}


ierr=sep3d_set_nh(sep3dname,ntr);
if(ierr!=SUCCESS){ free(header_list);return(ierr);}

if(info->nkeys_in==0) { /* regular dataset with headers created on the fly*/
	temp_headers=(int *) alloc (sizeof(int));
	drn=NO_DRN;
  nkeys_read=0;
}
else{

  /*get the tag for the header file and read in the headers */
  if(0!=fget_header_format_tag(tag,temp_ch)){
	  free(header_list);
	  return(sepwarn(FAIL_OTHER,
	   "read_headers:trouble reading header_format_tag from %s \n",tag));
  }
  /*check to see if a data record number exists */
  nkeys_read=info->nkeys_in;
  if(0==sep_get_key_index(tag,"data_record_number",&drn)){ /*if drn exists */
     nkeys_read++;
     drn--;
   }
	 else drn=NO_DRN;
   temp_headers=(int *) alloc (sizeof(int)* ntr * (nkeys_read));
	


}

if(info->nkeys_in>0){ /*not a faked regular dataset */
  ierr=sep3d_read_list(temp_ch,nkeys_read,nkeys_read,0,1,4,ntr,header_list,(char*)temp_headers);
  if(ierr!=0){ 
     free(header_list); return(ierr);
  }
}

/*PUT HEADERS INTO INTERNAL STRUCTURE */

/* if same_record_number =1 then we should ignore the drn, so we will just fix it */
if(0==auxpar("same_record_number","d",&same,tag)) same=1;

/*What is stored on disk is not what we want to store in the structure.
  Two possibilities: DRN in data and/or additional on the fly keys */
c=-1; /*pointer for the read in array*/
d=-1; /*pointer for the disk stored array */
for(i2=0;i2<info->nh;i2++){
	for(i1=0; i1 < nkeys_read; i1++){
			c++;
			if(i1 != drn){ 
        d++; info->headers[d] = temp_headers[c];
      }
			else if(same==0){
       info->drn[i2]=temp_headers[c]; 
      }
			else{ 
info->drn[i2]=header_list[i2]; }
	}
  if(info->nkeys!=info->nkeys_in) d+=info->nkeys-info->nkeys_in;
} 
free(temp_headers);


/*handle the case where we didn't have headers */
if(drn==NO_DRN) {
     for(i2=0; i2<info->nh; i2++){
         info->drn[i2]=header_list[i2];
     }
}

if(0!=calc_additional_headers(info,&info->nwind[1],&info->fwind[1],
   &info->jwind[1]))
    return(sepwarn(NOT_MET,"trouble calculating additional headers \n"));

*nh=info->nh;
free(header_list);
return (SUCCESS) ;
} 



/*
sep3d_read_list
<
USAGE
ierr=sep3d_read_list(char *tag, int ngrid, int nwind, int fwind, int jwind, int esize, int ntr, int *list, char *array)


INPUT PARAMETERS
tag            -   char*    tag to read from
ngrid          -   int      size of data axis 1
nwind          -   int      number of elements to grab
fwind          -   int      first element ot grab
jwind          -   int      sampling between elements
esize          -   int      esize of dataset
list           -   int*     list of traces to read
ntr            -   int      number of traces toread

OUTPUT PARAMETERS
array          -   char*     data to read


RETURN VALUES
0    =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Read portion of axis from a trace list


>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_read_data(char *tag, char *sep3dname ,int nwind,int fwind,int jwind,
 char *array)
_XFUNCPROTOEND
#endif
{
sep_3d *info;
int esize,ntr,*list,ierr,i,idim,icart[12];
char temp_ch[2048],temp2_ch[2048];;



info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */

if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));


info->nwind[0]=nwind; info->fwind[0]=fwind; info->jwind[0]=jwind;
esize=SEP3D_get_esize(info);
if(info->file_format==REGULAR){
  for(i=1,ntr=1; i <info->ndims; i++){
    if(info->nwind[i] == -1) ntr=ntr*info->n[i];
    else ntr=info->nwind[i]*ntr;
  }
  if(0!=SEP3D_alloc_coord(info,ntr)){
   return(sepwarn(FAIL_OTHER,"trouble allocating coords \n"));
  }
  list=(int*) malloc(sizeof(int)*ntr);
  if(info->nwind[0]==-1) { /*read in the begining of the data*/
    for(i=0; i < ntr; i++) list[i]=i+1;
  } 
  else{
     if(0!=SEP3D_wind_coords(info,&ntr))
       return(sepwarn(NOT_MET,"trouble initializing coordinates \n"));
     for(i=0; i < ntr; i++){
        list[i]= (int)(info->coord[i])+1;
     }
  }
}
else{
  if(info->nh<1) {
    return(sepwarn(NOT_MET,
     "irregular dataset (%s) and headers haven't been read\n", sep3dname));
  }
  list=(int*) malloc(sizeof(int)*info->nh);
  memcpy((void*)list,info->drn,info->nh*sizeof(int));
  ntr=info->nh;
}  


if(0!=sep3d_read_list(tag,info->n[0],nwind,fwind,jwind, esize,ntr,list,array)){
  free(list);
  sprintf(temp_ch,"Failed window read parameters thread=%d \n",sep_thread_num());
 for(i=0; i < info->ndims; i++){
   sprintf(temp2_ch,"axis=%d ng=%d nw=%d fw=%d jw=%d\n",i,info->n[i],
     info->nwind[i],info->fwind[i],info->jwind[i]);
   strcat(temp_ch,temp2_ch);
 }
 return(sepwarn(NOT_MET,"%s",temp_ch));

}
free(list);
return(SUCCESS);

}




_XFUNCPROTOBEGIN
int sep3d_read_list(char *tag, int ngrid,int nwind,int fwind,int jwind,int esize,
 int ntr, int *list, char *array)
_XFUNCPROTOEND
{
int tempi,i1,n,block,n2;
char temp_ch[10024];
int read,ierr,i2,new,fw[2],nw[2],ng[2],jw[2];
int buf_read,ib,ia,loc1,loc2,i0,loca,locb,two;


if(ntr >0 && nwind >0){
 if(0==auxpar("in","s",temp_ch,tag)) 
   return(sepwarn(NOT_MET,"attempt to read from tag=%s which has no in specified\n",tag));
  if(0==strcmp(temp_ch,"-1")) return(sepwarn(NOT_MET,
   "attempt to read from tag=%s which has no in specified\n",tag));

}

fw[0]=fwind; jw[0]=jwind; nw[0]=nwind; ng[0]=ngrid;
jw[1]=1;



/*do a buffered read  */
if(jwind!=1 || fwind !=0 || nwind != ngrid) buf_read=YES;
else buf_read=NO;

/*if(buf_read==YES) n2=NBUF/ngrid/esize-1;*/
n2=ntr;

read=0;two=2;
 
for(i2=0;i2<ntr;i2++){
   for (i1=i2; i1<(ntr-1); i1++){
       if((i1-i2)==n2||(list[i1]!=(list[i1+1]-1))) break;
   }
  block=i1-i2+1;
  new=(list[i2]-1);
  fw[1]=new; nw[1]=block; ng[1]=nw[1]+fw[1];


   
  if(0!=sreed_window_new(tag,two,ng,nw,fw,jw,esize,(array+i2*nw[0]*esize))){
    sprintf(temp_ch,"ng[0]=%d nw[0]=%d fw[0]=%d jw[0]=%d ng[1]=%d nw[1]=%d fw[1]=%d jw[1]=%d \n",
    ng[0],nw[0],fw[0],jw[0],ng[1],nw[1],fw[1],jw[1]);
    return(sepwarn(NOT_MET,"trouble reading block from tag %s\nread:%s",tag,temp_ch));
  }
   
  i2=i1;
}

return(SUCCESS);
}

_XFUNCPROTOBEGIN
int wind_to_helix(int ndim,int *ngrid,int *nwind, int *fwind, int *jwind, long long *index)
_XFUNCPROTOEND
{
int i,n123,idim;
int *block_wind,*block_dat;
int icoord,j,ih;

block_wind=(int*)malloc(sizeof(int)*ndim);
block_dat=(int*)malloc(sizeof(int)*ndim);


for(i=0,n123=1; i < ndim; i++) n123=n123*nwind[i];
for(i=1,block_wind[0]=1;i < ndim; i++) block_wind[i]=block_wind[i-1]*nwind[i-1];
for(i=1,block_dat[0]=1; i < ndim; i++) block_dat[i]=block_dat[i-1]*ngrid[i-1];

for(ih=0; ih< n123; ih++){
  /*convert from helix to wind space*/
  index[ih]=0;
  for(idim=0; idim < ndim; idim++){
    icoord=(ih/block_wind[idim])%nwind[idim]*jwind[idim]+fwind[idim];
    index[ih]+=(long long)icoord*block_dat[idim];
  } 
}

free(block_wind); free(block_dat);
return(SUCCESS);
}
