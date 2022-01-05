#define SET_SDOC 1
/*<

Date Created:Sun Aug 16 21:14:51 PDT 1998

*/	 

#include <seplib.h> 
#include <superset.h>
#include "superset_internal.h"
#define NBUF 100000
#ifndef YES
#define YES 1
#endif
#ifndef NO
#define NO 0
#endif
static int drn_warning=0, grid_warning=0,header_warning=0,order_warn=0;
int rite_list_l_l(char *tag, int ngrid,int nwind,int fwind,int jwind,int esize,
 int ntr, long long *list, char *array);


/*
<
sep3d_rite

USAGE
ierr= sep3d_rite(char *tag, char *sep3dname, int *nwind, int *fwind, int *jwind, char *data,int ntraces, int write_data,int write_headers,int write_grid) 

INPUT PARAMETERS
tag        -    char*   tag to write out too
sep3dname  -    char*   pointer to sep_3d structure
nwind      -    int*    number of traces/headers to write
fwind      -    int*   	first traces/headers to write
jwind      -    int*    skip of traces/headers to write
data       -    char*   Data to write
ntraces    -    int     Number of traces to write
write_data -    int     Whether(1) or not to write the data
write_headers - int     Whether(1) or not to write the headers
write_data -    int     Whether(1) or not to write the grid

RETURN VALUES
SUCCESS  =    if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Writes out headers and possibly data.  Use a C or F90 overlay 

>
*/



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_rite(char *tag, char *sep3dname, int *nwind, int *fwind, int *jwind, void *data, int ntraces, int write_data,int write_headers, int write_grid) 
_XFUNCPROTOEND
#else
int sep3d_rite(tag,sep3dname,nwind,fwind,jwind,data,ntraces, write_data,write_headers,write_grid) 
char *tag,*sep3dname; 
int *nwind,*fwind,*jwind, write_data;
void *data;
int write_headers,write_grid;
#endif 
{ 
int i1,esize,*idim,n2,i2,temp2i,*header_list,nelem,loc,*temp_headers;
int loc2;
sep_3d *info;
char temp_ch[1024];
int a,b,ierr,nelem_window,nelem_header,remainder,i,ia;
int *cloc,nkey_out,*grid;



info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));

for(i=0; i < info->ndims; i++){
   if(info->n[i] >1){
     info->nwind[i]=nwind[i]; 
     info->fwind[i]=fwind[i]; 
     info->jwind[i]=jwind[i];
   }
   else{
     info->nwind[i]=1;
     info->fwind[i]=0;
     info->jwind[i]=1;
   }
}


return(SEP3D_rite(tag,info,data,write_data,write_headers,write_grid));
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int SEP3D_rite(char *tag, sep_3d *info, void *data, int write_data,int write_headers, int write_grid) 
_XFUNCPROTOEND
#endif
{
int i1,esize,*idim,n2,i2,temp2i,nelem,loc,*temp_headers;
int loc2;
char temp_ch[1024];
int a,b,ierr,nelem_window,nelem_header,remainder,i,ia;
int *cloc,nkey_out,*grid;
int ntraces,nh;
long long  *header_list;


if(info->file_format!=GRID){
    if(0!=SEP3D_wind_coords(info,&nh))
      return(sepwarn(NOT_MET,"trouble creating coordinate map tag=%s \n",
        info->name));

    if(nh==0) return(SUCCESS);
}
ntraces=info->ncoord;



 esize=SEP3D_get_esize(info);
nelem_window=1; 
for(i=0; i<info->ndims;i++){
   if(i>0) nelem_window=nelem_window*info->nwind[i];
}

  if(info->file_format==REGULAR){
    if(write_data==0) return(SUCCESS);
		cloc=(int*)alloc(sizeof(int)*(info->ndims-1));	
		header_list=(long long*)alloc(sizeof(long long)*nelem_window);
		for(i=0; i < nelem_window; i++){
			header_list[i]=(long long)info->coord[i]+1;
		}
     if(0!=rite_list_l_l(tag,info->n[0],info->nwind[0],info->fwind[0],info->jwind[0],
       esize,nelem_window,header_list,data)){
		    free(cloc);
                    free(header_list);
        return(sepwarn(NOT_MET,"trouble  writing list tag=%s \n",
         tag));
     }
     free(cloc);
     free(header_list);
     return(SUCCESS);
    }


if(write_data==1 ) info->wrote_data=YES;
if(write_headers==1 ) info->wrote_headers=YES;

/*first lets check the users sanity */
if(info->headers!=SEPNULL && header_warning==0 && write_headers==0){
	header_warning=1; sepwarn(0,
		"WARNING: You  have allocated  headers in [%s] but asked not to write it \n",
     tag);
}
if(write_data==0 && write_headers ==0 && write_grid==0)
	return(sepwarn(INVALID_DATA,
   "Can not call sep3d_rite and not write out grid, data, or headers \n"));

if(write_headers==0 && write_grid==1)
	return(sepwarn(INVALID_DATA,
   "Can not write out the grid without also writing out the headers \n"));




if(write_grid==1){
if(order_warn==1 )
		return(sepwarn(INVALID_DATA,
		 "you can't write the grid after doing a header write where headers were skipped (%s)\n",tag));
}
else if(write_data==1 && info->file_format==HEADER){
	if(nelem_window!=info->nh)
	 return(sepwarn(INVALID_DATA,
  "The number of header elements(%d) stored is not equal to the window size(%d) tag=(%s)\n",
      nelem_window,info->nh,tag));
}



/*If we are writing out the grid or headers we need a starting position
  for the headers.  We can get the current position in the file by
  seeking from current position 0 blocks  */
  

if(write_headers==1){
	if(SUCCESS!=fget_header_format_tag(tag,temp_ch))
	 return(sepwarn(FAIL_OTHER,"trouble getting header format tag for %s\n",tag));
	/*get current position  in a tricky manner (posible error if 
	we have written out a portion of a key or wrong number of keys*/

        a=info->nkeys+1;
	if(1==info->in_order) a=info->nkeys;
	ierr=file_position(temp_ch,a*4,&loc,&remainder);
	if(ierr!=0) return(ierr);
	if(remainder!=0)return(sepwarn(INVALID_DATA,
		"inconsistent write of headers. Nkeys must be constant through write(%s)\n",tag
	));
	if(0!=sep3d_drn_set(info->name)){
          if(info->drn!=SEPNULL){
            for( i=0; i < info->nh; i++) info->drn[i]=loc+i+1;

          }
        }
}

if(write_grid==1){ /*TIME TO WORK ON THE GRID*/

	/*the grid and header must by synched so we can ignore the actual number
   in the grid, just check whether it is positive */
	nelem_header=loc;
	header_list=(long long*)alloc(sizeof(long long)*nelem_window);
	ia=0;
  grid=(int*) malloc(sizeof(int)*nelem_window);
   if(0!=sep3d_grab_grid_block(info->name,grid)){
    free(header_list);free(grid);
    return(sepwarn(NOT_MET,"trouble grabbing grid block for tag=%s \n",
     info->name));
    }

	for(i=0; i < nelem_window;i++){
		if(grid[i] >0) {
			nelem_header++;
			grid[i]=nelem_header;
			header_list[ia]=nelem_header;ia++;
		}
	}
	nelem_header-=loc;
	if(nelem_header!=info->nh){
		free(header_list);free(grid);
		return(sepwarn(INVALID_DATA,"The grid specifies (%d) a different number of headers than what is set (%d) for tag (%s) \n",nelem_header,info->nh,tag));

	}

	ierr=sep_put_grid_window(tag,&info->ndims,&(info->n[1]),&info->nwind[1],
   &info->fwind[1], &info->jwind[1],grid);
	if(ierr!=SUCCESS) {
     free(grid);free(header_list);
     return(ierr);
   }

  free(grid);
}
else nelem_header=nelem_window;



/*now we will deal with the headers */
if(write_headers==1){
	if(write_grid==0){
		if(order_warn==0){
			/*first check if we are doing an asynchronous write of the headers */
			c2h(&loc2,&(info->n[1]),info->ndims-1,&info->fwind[1]);
			if(loc!=loc2) order_warn=1;
		}
		cloc=(int*)alloc(sizeof(int)*(info->ndims-1));	
		header_list=(long long*)alloc(sizeof(long long)*nelem_window);
		for(i=0; i < nelem_window; i++){
			h2c(i,&(info->n[1]),info->ndims-1,cloc);
			for(ia=0;ia<info->ndims-1;ia++) 
        cloc[ia]=info->fwind[ia+1]+cloc[ia]*info->jwind[ia+1];
			c2h(&loc2,&(info->n[1]),info->ndims-1,cloc);
			header_list[i]=loc2+1;
		}
		free(cloc);
	}
	/* if we have a drn we nead to create a temporary buffer */
	if(1!=info->in_order){ /* IF WE DON'T HAVE A DRN */
		nkey_out=info->nkeys+1;
		temp_headers=(int*)alloc(nelem_header*(info->nkeys+1)*sizeof(int));
		a=-1; b=-1;
		for(i2=0; i2 < nelem_header; i2++){
			for(i1=0; i1< info->nkeys;i1++){
				a++; b++;
				temp_headers[a]=info->headers[b];
			}
			a++;
			temp_headers[a]=info->drn[i2];
		}
	}
	else{ nkey_out=info->nkeys; temp_headers=info->headers;}

	if(order_warn==0 && write_grid==0){
		/*check for asynchronous write that will break a grid write */
		loc2=loc;
		for(i=0; i < nelem_header; i++){
			if(header_list[i]!=loc2+1) { order_warn=1; break;}
			else loc=header_list[i];
		}
	}

  info->ntraces_wrote+=nelem_header;
	ierr=rite_list_l_l(temp_ch,nkey_out,nkey_out,0,1,4,nelem_header,
         header_list,(char*)temp_headers);
	if(0==sep3d_drn_set(info->name)) free(temp_headers);
/*AAAA*/
/*	free(temp_headers);*/
	if(ierr!=SUCCESS){ free(header_list);return(ierr);}
}
if(write_data==1){
	if(write_headers==1){
		if(ntraces!=info->nh) {
			free(header_list);
			return(sepwarn(INVALID_DATA,
			"sep3d_rite:number of coordinates not equal to the number of headers nh=%d ntraces=%d \n",info->nh,ntraces));

		}
	}
	if(write_headers==1 && 0==sep3d_drn_set(info->name)){
   	/*then we are going to write using the supplied drn */
		ierr=sep3d_rite_list(tag,info->n[0],info->n[0],info->fwind[0],info->jwind[0],esize,ntraces,
			info->drn, (char*)data);
	}
	else if(write_headers==1){
   	/*then we are going to write using the same list as the headers */
		ierr=rite_list_l_l(tag,info->n[0],info->n[0],info->fwind[0],info->jwind[0],esize,ntraces,
			header_list, (char*)data);
	}
	else{
   	/*we are just writing out data*/
		cloc=(int*)alloc(sizeof(int)*(info->ndims-1));	
		header_list=(long long*)alloc(sizeof(int)*nelem_window);
		for(i=0; i < nelem_window; i++){
			h2c(i,&(info->n[1]),info->ndims-1,cloc);
			for(ia=0;ia<info->ndims-1;ia++) cloc[ia]=info->fwind[ia+1]+cloc[ia]*info->jwind[ia+1];
			c2h(&loc2,&(info->n[1]),info->ndims-1,cloc);
			header_list[i]=loc2;
		}
		free(cloc);
	}
	free(header_list);
	if(ierr!=0){ return(ierr);}
}
else free(header_list);
return (SUCCESS);
} 




/*
sep3d_rite_list
<
USAGE
ierr=sep3d_rite_list(char *tag, int ngrid, int nwind, int fwind, int jwind, int esize, int ntr, int *list, char *array)


INPUT PARAMETERS
tag            -   char*    tag to read from
ngrid          -   int      size of data axis 1
nwind          -   int      number of elements to grab
fwind          -   int      first element ot grab
jwind          -   int      sampling between elements
esize          -   int      esize of dataset
list           -   int*     list of traces to read
ntr            -   int      number of traces toread
array          -   char*     data to read


RETURN VALUES
0    =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Writes out  portion of axis from a trace list

>
*/




_XFUNCPROTOBEGIN
int rite_list_l_l(char *tag, int ngrid,int nwind,int fwind,int jwind,int esize,
 int ntr, long long *list, char *array)
_XFUNCPROTOEND
{
int tempi,i1,n,read,block,n2;
int ierr,i2;
int buf_rite,ib,ia,loc1,loc2,is;
long long new;

n2=ntr;
if(n2==0) return(SUCCESS);
if(n2<1) return(sepwarn(INVALID_DATA,
  "Can not even hold a single trace in memory  (%s) \n",tag));

read=0;

for(i2=0;i2<ntr;i2++){
   for (i1=i2; i1<(ntr-1) ; i1++){
			 if((i1-i2)==n2||(list[i1]!=(list[i1+1]-1))){
          break;
       }
   }
  block=i1-i2+1;
  new=(list[i2]-1);
	if(0> sseek_l_l(tag,new*(long long)esize*(long long)ngrid,0))
		return(sepwarn(FAIL_OTHER,"sseek failed (block=%d, to=%lld) (tag=%s) (iloc=%d)\n",
     esize*ngrid,new,tag,list[i2]));

 	ierr=srite(tag,&array[i2*ngrid*esize],ngrid*esize*block);
	
  if(block*ngrid*esize!=ierr)
	   return(sepwarn(FAIL_OTHER,
   "(rite_window) trouble writing in data (err=%d) starting_trace=%d n_bytes=%d (%s)\n",
   ierr,new,esize*block*ngrid,tag));
  i2=i1;
}

return(SUCCESS);
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_rite_list(char *tag, int ngrid,int nwind,int fwind,int jwind,int esize,
 int ntr, int *list, char *array)
_XFUNCPROTOEND
#else
int sep3d_rite_list(tag,ngrid,nwind,fwind,jwind,esize,ntr,list,array)
char *tag; 
int ngrid, nwind, fwind,jwind, ntr,*list,esize;
char *array;
#endif 
{
int tempi,i1,n,read,block,n2;
int ierr,i2,new;
int buf_rite,ib,ia,loc1,loc2,is;

n2=ntr;
if(n2==0) return(SUCCESS);
if(n2<1) return(sepwarn(INVALID_DATA,
  "Can not even hold a single trace in memory  (%s) \n",tag));

read=0;

for(i2=0;i2<ntr;i2++){
   for (i1=i2; i1<(ntr-1) ; i1++){
			 if((i1-i2)==n2||(list[i1]!=(list[i1+1]-1))){
          break;
       }
   }
  block=i1-i2+1;
  new=(list[i2]-1);
	if(0> sseek_block(tag,new,esize*ngrid,0))
		return(sepwarn(FAIL_OTHER,"sseek failed (block=%d, to=%d) (tag=%s) (iloc=%d)\n",
     esize*ngrid,new,tag,list[i2]));

 	ierr=srite(tag,&array[i2*ngrid*esize],ngrid*esize*block);
	
  if(block*ngrid*esize!=ierr)
	   return(sepwarn(FAIL_OTHER,
   "(rite_window) trouble writing in data (err=%d) starting_trace=%d n_bytes=%d (%s)\n",
   ierr,new,esize*block*ngrid,tag));
  i2=i1;
}

return(SUCCESS);
}



