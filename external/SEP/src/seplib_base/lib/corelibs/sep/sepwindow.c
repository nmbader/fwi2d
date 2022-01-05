/*$

=head1 NAME

sreed_window - read a window of a seplib dataset

=head1 SYNOPSIS

C<int sreed_window(tag_history,n_dim_cube,n_cube,n_wind,f_wind,j_wind,esize,values)>

=head1 INPUT PARAMETERS

=over 4

=item   char* - tag_history     

        tag of History File

=item    int*  - n_dim_cube	     

         Number of Dimensions in the cube

=item    int*  - n_cube		       

         vector of length n_dim_cube length of the axes in the cube

=item    int*  - n_wind          

         vector of length n_dim_cube axes length after windowing

=item    int*  - f_wind          

         vector of length n_dim_cube index of first elements C<(0<= o < n)>

=item    int*  - j_wind          

         vector of length n_dim_cube sampling rate along axes C<(1<= j < n)>

=item    int   - esize	   

         number of bytes per element

=back

=head1 OUTPUT PARAMETERS

=over 4

=item   void* - values          

        array of values to be read

=back


=head1 RETURN VALUE

 0 = if successful
 
 -1 = if fails for other reasons

 -2 = if the values in n_wind, f_wind, j_wind are incorrect

=head1 DESCRIPTION

It reads a subset (window) of a Seplib cube

=head1  SEE ALSO

L<srite_window>


=head1 LIBRARY

B<sep>

=cut


>*/

/*$

=head1 NAME

srite_window - write a seplib window


=head1 SYNOPSIS

C<int srite_window(tag_history,n_dim_cube,n_cube,n_wind,f_wind,j_wind,esize,values)>

=head1 INPUT PARAMETERS

=over 4


=item    char* - tag_history     

         tag of History File

=item    int*  - n_dim_cube	     

          Number of Dimensions in the cube

=item    int*  - n_cube		       

         Vector of length n_dim_cube length of the axes in the cube

=item    int*  - n_wind          

         Vctor of length n_dim_cube axes length after windowing

=item    int*  - f_wind          

         Vector of length n_dim_cube index of first elements C<(0<= o < n)>

=item    int*  - j_wind          

         Vector of length n_dim_cube sampling rate along axes C<(1<= j < n)>

=item    int   - esize	   

         number of bytes per element

=item    void* - values          

         array of values to be written.


=back

=head1 RETURN VALUE

 0 = if successful

 -1 = if fails for other reasons

 -2 = if the values in n_wind, f_wind, j_wind are incorrect

=head1 DESCRIPTION

It writes a subset (window) of a Seplib cube

=head1 SEE ALSO

L<srite_window>

=head1 LIBRARY

B<sep>

=cut


>*/

/*
KEYWORDS
   read, write, window

SEE ALSO
   sreed srite 

AUTHOR
   Biondo Biond1 , September 1995
	 Modified: Robert Clapp : Added Prototypes
	 Modified: Robert Clapp : Fixed to +2GB  (Aug 98)
	 Modified: Robert Clapp : Fixed to +8GB  (june 99)
	 Modified: Robert Clapp : Buffered reads, improve readibility
         Modified: Adam Halpert : Changed "Illegal window axis" calculation (Feb 10)

*/
#include <sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <stdio.h>
#include "sepstream.h"
#include "streamlist.h"
#include "sep_main_internal.h"
#include <sep_main_external.h>
#include <math.h>
#include <assert.h>

enum{
READ=-1,/*Call core sepwindow with read option */
BUFFER_READ=-2,/*Call core sepwindow with read option */
WRITE=1,/*Call core sepwindow with write option */
CREATE_GRID =0,/*We have a regular grid so just construct the header numbers*/
BUF_SIZE=5000000,/*Internal buffer for buffered io */
YES = 1,/*logical int definiton */
NO = 0 /*logical int definiton */
};


#if defined(SGI)
#include <sys/types.h>
#include <sys/socket.h>
#include <rpc/rpc.h>
#else
#include <rpc/types.h>
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int sep_reed_rite_window_rec(int,const char*,const int,const int,const int,const int,const int,const int, 
 const int*,const int*,const int*,const int*,const int*,const int*,const int,void*);
void sep_comp_index_grid(int,int, int*);
_XFUNCPROTOEND 
#else
int sep_reed_rite_window_rect();
void sep_comp_index_grid();
#endif


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int sreed_window( const char *tag_history, const int *n_dim_cube, const int *n_cube, 
		const int *n_wind,const int *f_wind,const int *j_wind, const int esize, void *values)
_XFUNCPROTOEND 
#else 
int sreed_window(tag_history, n_dim_cube, n_cube, n_wind, f_wind, j_wind, esize, values)
char *tag_history;
int  *n_dim_cube;
int  *n_cube;
int  *n_wind, *f_wind, *j_wind;
int  esize;
char *values;
#endif

{
  //fprintf(stderr,"in tag history %s \n",tag_history);
    return (sep_window(READ,tag_history, n_dim_cube, n_cube, n_wind, f_wind, j_wind,esize, values));
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int srite_window(const char *tag_history, const int *n_dim_cube, const int *n_cube, 
		const int *n_wind,const int *f_wind,const int *j_wind, const int esize, void *values)
_XFUNCPROTOEND 
#else 
int srite_window(tag_history, n_dim_cube, n_cube, n_wind, f_wind, j_wind, esize, values)
char *tag_history;
int  *n_dim_cube;
int  *n_cube;
int  *n_wind, *f_wind, *j_wind;
int  esize;
char *values;
#endif

{
    return (sep_window(WRITE,tag_history, n_dim_cube, n_cube, n_wind, f_wind, j_wind, esize, values));
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int sep_window( int mode,const char *tag_history,
                const int *n_dim_cube, const int *n_cube, 
		const int *n_wind,const int *f_wind,const int *j_wind, const int esize, void *values)
_XFUNCPROTOEND 
#else 
int sep_window(mode,tag_history, n_dim_cube, n_cube, n_wind, f_wind, j_wind, esize, values)
int mode;
char *tag_history;
int  *n_dim_cube;
int  *n_cube;
int  *n_wind, *f_wind, *j_wind;
int  esize;
char *values;
#endif
{
int i,ndim,iax;
int max_dim_found=NO;
int *nwind,*fwind,*jwind;
int *ngrid, *narray, *nbuffer,*fbuffer,*jbuffer;
int last_buf_ax,first_buf_ax;
long long block;
int  seek_to,cur_pos;
int  block_int,seek_int,array_pt,nelem,n_axis,n_buf_axis;
char buffer[BUF_SIZE];
int nloop,tempi,nloop_buf;
int seek_from, copy_size;
int pos,pos2,*ipos1,*ipos2,ia,ib,ierr2;
double derr;
char *pval,*pval2;
int *int_pt,seek_block_elem;
int last_rel; /* last axis with n>1 */
sep_off_t seek_to_off_t,current_pos,big_block,seek_distance,block_distance;
int blocks,remainder;
int beg_seek;  /*the ammount we need to move forward on the biggest cont axies*/
int seek_block_size; /*the block size of how much we need to move forward */
long long ierr;
ndim=*n_dim_cube;block=esize;

ngrid=malloc(sizeof(int)*ndim);
nwind=malloc(sizeof(int)*ndim);
fwind=malloc(sizeof(int)*ndim);
jwind=malloc(sizeof(int)*ndim);
last_rel=0;
for(i=0; i < ndim; i++){
	nwind[i]=1; fwind[i]=0; jwind[i]=1; ngrid[i]=1;
	if(n_wind[i] > 1) last_rel=i;
}

n_axis=0; 
first_buf_ax=-1;
last_buf_ax=0;
seek_block_size=-1;
seek_block_elem=1;
beg_seek=0;
for(i=0;i<ndim;i++){
	if(f_wind[i] < 0||j_wind[i]<=0||f_wind[i]+j_wind[i]*(n_wind[i]-1)>=n_cube[i]){
		free(ngrid); free(nwind);free(fwind);free(jwind);
		return(sepwarn(-2,
       "Invalid window parameters (%s) iax=%d ng=%d nw=%d fw=%d jw=%d \n",
       tag_history,i+1,n_cube[i],n_wind[i],f_wind[i],j_wind[i]));

		}
	if(mode==CREATE_GRID){
		nwind[n_axis]=n_wind[i];
		fwind[n_axis]=f_wind[i];
		jwind[n_axis]=j_wind[i];
		ngrid[n_axis]=n_cube[i];
		n_axis++;
	}
	else{
		if(max_dim_found==NO){
			if(block * (long long) n_wind[i] > MAX_INT_SIZE || i>last_rel){
				 max_dim_found=YES; last_buf_ax=i;
			}
			else if(n_wind[i]==n_cube[i]){ /*the windows spans this dimension */
				if(mode!=BUFFER_READ){
					block=(long long)block* (long long)n_wind[i];
				}
				else if(block*n_wind[i] <= BUF_SIZE){ /*we are buffer ing and it fits*/
					block=(sep_file_size_t)block* (sep_file_size_t)n_wind[i];
				}
				else{  
					max_dim_found=YES; ;last_buf_ax=i;
				}
			   }
			else if(j_wind[i]==1){ /*we have a consecutive chunk of this window*/
				if(mode!=BUFFER_READ){
					block=(sep_file_size_t)block* (sep_file_size_t)n_wind[i];
				}
				else if(block*(sep_file_size_t)n_cube[i] <= BUF_SIZE){
					block=(sep_file_size_t)block* (sep_file_size_t)n_cube[i];
				}
				else {
					max_dim_found=YES; last_buf_ax=i;}
			}
			else if((mode==READ || mode==BUFFER_READ) && i<=last_rel&&
         (long long)BUF_SIZE >= (long long)n_cube[i]*(long long)block){
				if(mode==READ){
					mode=BUFFER_READ;
					first_buf_ax=i;
				}
				block=block*(sep_file_size_t)n_cube[i];
			}
			else{
        last_buf_ax=i; max_dim_found=YES;}/*we definitely have to loop */
		}
		if(max_dim_found==YES){
			nwind[n_axis]=n_wind[i]; fwind[n_axis]=f_wind[i]; jwind[n_axis]=j_wind[i];
			ngrid[n_axis]=n_cube[i]; n_axis++;
		}
	  else if(n_wind[i] != n_cube[i] && j_wind[i]==1){
			if((mode==READ)){ /*we might be able to buffer this axis*/
        if(BUF_SIZE>= (sep_file_size_t) n_cube[i]*block && i <= last_rel){
						mode=BUFFER_READ;
						first_buf_ax=i;
						block=block*((long long)(n_cube[i]/n_wind[i]));
				}
				else{  /*we are reading a continuous section of an axis */ 
         max_dim_found=YES; beg_seek=f_wind[i]; 
				 seek_block_size=block/n_wind[i];				
				 seek_block_elem=n_cube[i];
				}
			}
			else if(mode!=BUFFER_READ) {
				max_dim_found=YES; beg_seek=f_wind[i];
			  seek_block_size=block/n_wind[i];				
			  seek_block_elem=n_cube[i];
			}
		}

	}
}

if(seek_block_size==-1) seek_block_size=block;





if(last_buf_ax==0){
	if(mode==BUFFER_READ)
		last_buf_ax=ndim; /*we can hold the whole data in memory*/
}
	
	n_buf_axis=last_buf_ax-first_buf_ax;

nloop=1; tempi=1;
for(i=0; i < n_axis; i++){
	 nloop=nloop*nwind[i]; /*total number of loops */
	 if((long long)tempi*(long long)ngrid[i] > MAX_INT_SIZE){
			free(ngrid); free(nwind);free(fwind);free(jwind);
			return(sepwarn(-2,"THE NUMBER OF ELEMENTS IS TO BIG TO BE HANDLED \n"));
		}
		else tempi=tempi*ngrid[i]; 
}



/* if we are doing a buffered read we to do a similar calculation*/
fbuffer=(int*) malloc(sizeof(int)*n_buf_axis);
jbuffer=(int*) malloc(sizeof(int)*n_buf_axis);
nbuffer=(int*) malloc(sizeof(int)*n_buf_axis);
narray=(int*) malloc(sizeof(int)*n_buf_axis);
if(mode==BUFFER_READ){
	copy_size=esize;
	for(i=0; i < first_buf_ax;i++) copy_size=copy_size*n_wind[i];
	nloop_buf=1;
	n_buf_axis=0;
	for(i=first_buf_ax; i < last_buf_ax; i++){
		fbuffer[n_buf_axis]=f_wind[i]; jbuffer[n_buf_axis]=j_wind[i];
		narray[n_buf_axis]=n_wind[i]; nbuffer[n_buf_axis]=n_cube[i];
		nloop_buf=nloop_buf*narray[n_buf_axis]; 
/*   ,n_buf_axis,fbuffer[n_buf_axis], jbuffer[n_buf_axis],  */
/*   narray[n_buf_axis],nbuffer[n_buf_axis],  nloop_buf,copy_size);*/
		n_buf_axis++;
	}
}


/*NOW IT IS TIME TO DO THE REAL WORK */
ipos1=(int *) malloc(sizeof(int)*MAX(n_axis,n_buf_axis));
ipos2=(int *) malloc(sizeof(int)*MAX(n_axis,n_buf_axis));

/*for(i=0; i<n_axis;i++){*/
/*	fprintf(stderr,"window pars fw=%d %d %d %d ng=%d \n",i,fwind[i],jwind[i],nwind[i],ngrid[i]);*/
/*}*/


if(mode==BUFFER_READ) beg_seek=0;

pval=(char*)values;
for(i=0; i <nloop; i++){
	h2c(i,nwind,n_axis,ipos1) ; /*calculate our current position in the window*/
	/*scale by the window parameters */
	for(ia=0; ia < n_axis; ia++){
		ipos2[ia]=fwind[ia]+jwind[ia]*ipos1[ia];
	}
	c2h(&pos,ngrid,n_axis,ipos2); /* convert to grid parameters */
	if(mode==CREATE_GRID){
		int_pt=(int*)pval;
		*int_pt=pos+1;
		pval=(char*)pval+block;
	}
	else{

		seek_to_off_t=(sep_off_t)pos*(sep_off_t)seek_block_elem+(sep_off_t)beg_seek;
		if(fabs(seek_to_off_t) >MAX_INT_SIZE){
			/*we need to be tricky.  We will figure out our relative position
       and then seek to where we want to go.  If the seek distance is greater
       than MAX_INT_SIZE we will do it in 2+ steps */
			big_block=(sep_off_t)seek_block_elem * (sep_off_t) seek_block_size;
			if(big_block >MAX_INT_SIZE)
				return(sepwarn(-1,"can not handle seek, blocks are to large \n"));
			if(0!=file_position(tag_history,seek_block_elem*seek_block_size,
        &blocks,&remainder));
       current_pos=big_block*(sep_off_t)blocks+(sep_off_t) remainder;
			seek_distance=seek_to_off_t*seek_block_size-current_pos;
			block_distance=seek_distance/(sep_off_t)seek_block_size; /*units to relative 
      seek in seek_block_size quantities*/
/*      fprintf(stderr,"now time to check %g %g %g \n",(long long)current_pos,(long long)seek_distance,block_distance);*/
			if(fabs(block_distance)>MAX_INT_SIZE){ /*if larger than max int */
				seek_to=(int)(seek_distance/big_block);
				derr=sseek_block_d(tag_history,seek_to,(int)big_block,SEEK_SET);
				if(derr<0)
					return(sepwarn(-1,"trouble seeking %d blocks of size %d in tag %d\n",
					seek_to,(int)big_block,tag_history));
				seek_distance-=(sep_off_t)seek_to * big_block;
				block_distance=seek_distance/(sep_off_t)seek_block_size; 
				seek_to=(int)block_distance;
				seek_from=SEEK_CUR;
			}
			else{
				 seek_to=(int)block_distance;
				seek_from=SEEK_CUR;
			}
		}
		else{ 
                seek_to=(int)seek_to_off_t; seek_from=SEEK_SET;}
	  derr=sseek_block_d(tag_history,seek_to,seek_block_size,seek_from);
/* fprintf(stderr,"seek stuff %g %d %d %d \n",derr,seek_to,seek_block_size,seek_from);*/
	  if(derr<0){
			  free(ngrid); free(nwind); free(fwind); free(jwind); free(ipos2);
        free(narray);free(nbuffer); free(fbuffer); free(jbuffer);free(ipos1);

				fprintf(stderr,"sepwindow, SEEK ERROR: big_block=%d cur=%17.12g \n seek_to_off=%17.12g  block_distance=%17.12g \n",(int)big_block,(double)current_pos,(double)seek_to_off_t,(double)block_distance);
			return(sepwarn(-1,"trouble seeking(2) %d blocks of size %d in tag %s\n",
				seek_to,seek_block_size,tag_history));
		}
		if(mode==WRITE){
			if(block!=srite_big(tag_history,pval,block)){
				free(ngrid); free(nwind); free(fwind); free(jwind); free(ipos2);
      	free(narray);free(nbuffer); free(fbuffer); free(jbuffer);free(ipos1);
				return(sepwarn(-1,"srite error to tag %s \n",tag_history));
			}
			pval=(char*)pval+block;
		}
		else if(mode==READ){
		  ierr=sreed_big(tag_history,pval,block);
			if(block!=ierr){
				free(ngrid); free(nwind); free(fwind); free(jwind); free(ipos2);
      	free(narray);free(nbuffer); free(fbuffer); free(jbuffer);free(ipos1);
				return(sepwarn(-1,"sreed error from tag %s \n",tag_history));
			}
			pval=(char*)pval+block;
		}
		else{
/*			seek_to=pos*seek_block_elem+beg_seek;*/
/*			fprintf(stderr,"asking to read  %d  loop=%d \n",block,nloop_buf);*/
/*			ierr2=sseek_block(tag_history,seek_to,seek_block_size,SEEK_SET);*/
			ierr=sreed(tag_history,buffer,block);
/*  fprintf(stderr,"check the seek read  %d %d \n",ierr,ierr2);*/
			if(block!=ierr){
				fprintf(stderr,"mode=%d i %d, nloop=%d  nloop_buf=%d \n",mode,i,nloop,nloop_buf);
				fprintf(stderr,"pos %d of elem %d  beg=%d\n",pos,seek_block_elem,beg_seek);
				fprintf(stderr,"seeking %d of size %d \n",seek_to,seek_block_size);
				fprintf(stderr,"reading %lld got=%lld \n",block,ierr);



				free(ngrid); free(nwind); free(fwind); free(jwind); free(ipos2);
      	free(narray);free(nbuffer); free(fbuffer); free(jbuffer);free(ipos1);
				return(sepwarn(-1,"sreed  buffer error for tag %s \n",tag_history));
			}
			for(ia=0; ia < nloop_buf; ia++){
				/*calculate our current offset in the array */
				h2c(ia,narray,n_buf_axis,ipos1) ; 
				/*scale by the window parameters */
				for(ib=0; ib < n_buf_axis; ib++) {
					ipos2[ib]=fbuffer[ib]+jbuffer[ib]*ipos1[ib];
				}
				c2h(&pos2,nbuffer,n_buf_axis,ipos2); /* convert to grid parameters */
/*				fprintf(stderr,"buff looop from=%d to=%d \n",ia,pos2);*/
				memcpy((void *)pval,(const void*)(buffer+copy_size*pos2),copy_size);
/*				int_pt=(int*)pval; fprintf(stderr,"the first %d, %d \n",*int_pt,copy_size);*/
				pval=(char*)pval+copy_size;
			}
		}
	}
}

free(ngrid); free(nwind); free(fwind); free(jwind); free(ipos2);
free(narray);free(nbuffer); free(fbuffer); free(jbuffer);free(ipos1);
return(0);
}




struct _wind{
  long long nc[9];
  long long nw[9];
  long long jw[9],fw[9];
  long long block;
  int jit;
  int fit;
  int nit;
  int bit;
  long bs[9];
};
typedef struct _wind wind;
#define BUFSIZE 50000000
char wind_buf[BUFSIZE];

void compress_it(const int ndim, const int *ncube, const int *nwind, const int *fwind, const int *jwind, const int esize,
   wind *comp, int type);
int window_it(const char *tag, wind *pars, const int write, char *buffer);

int sreed_window_new(const char *tag_history, const int ndim, const int *n_cube, 
		const int *n_wind,const int *f_wind,const int *j_wind, const int esize, void *values){
		
 wind compress;
  compress_it(ndim,n_cube,n_wind,f_wind,j_wind,esize,&compress,0);
   return(window_it(tag_history,&compress,0,values));
}
int srite_window_new(const char *tag_history, const int ndim, const int *n_cube, 
		const int *n_wind,const int *f_wind,const int *j_wind, const int esize, void *values){
		
 wind compress;
  compress_it(ndim,n_cube,n_wind,f_wind,j_wind,esize,&compress,1);
   return(window_it(tag_history,&compress,1,values));
}

int window_it(const char *tag, wind *pars, int write,char *buffer){
  int ierr;
  long long i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i0;
  long long seek_me[9],seek_to[9],seek_it;
  long long pos=0;

  for(i8=0; i8  < pars->nw[8]; i8++){
    seek_me[8]=(pars->fw[8]+pars->jw[8]*i8)*pars->bs[8];
    for(i7=0; i7 < pars->nw[7];i7++){
      seek_me[7]=(pars->fw[7]+pars->jw[7]*i7)*pars->bs[7];
      seek_to[7]=seek_me[8]+seek_me[7];
      for(i6=0; i6 < pars->nw[6];i6++){
        seek_me[6]=(pars->fw[6]+pars->jw[6]*i6)*pars->bs[6];
        seek_to[6]=seek_me[7]+seek_me[6];
        for( i5=0; i5 < pars->nw[5]; i5++){
          seek_me[5]=(pars->fw[5]+pars->jw[5]*i5)*pars->bs[5];
          seek_to[5]=seek_me[6]+seek_me[5];
          
          for(i4=0; i4 < pars->nw[4]; i4++){
          
            seek_me[4]=(pars->fw[4]+pars->jw[4]*i4)*pars->bs[4];
            seek_to[4]=seek_me[4]+seek_me[5];
            for( i3=0; i3 < pars->nw[3]; i3++){
              seek_me[3]=(pars->fw[3]+pars->jw[3]*i3)*pars->bs[3];
              seek_to[3]=seek_me[3]+seek_me[4];
              for(i2=0; i2 < pars->nw[2]; i2++){
                seek_me[2]=(pars->fw[2]+pars->jw[2]*i2)*pars->bs[2];
                seek_to[2]=seek_me[2]+seek_me[3];
                for(i1=0; i1 < pars->nw[1];i1++){
                  seek_me[1]=(pars->fw[1]+pars->jw[1]*i1)*pars->bs[1];
                  seek_to[1]=seek_me[2]+seek_me[1];
                  for(i0=0; i0 < pars->nw[0]; i0++){
                    seek_it=seek_to[1]+(pars->fw[0]+pars->jw[0]*i0)*pars->bs[0];
                    if(seek_it!=sseek_l_l(tag,seek_it,0)){
                       fprintf(stderr,"Trouble seeking tag=%s to=%g \n",
                          tag,(double)seek_it);
                       return 1;
                    }
                    if(write==1){
                      if((int)pars->block!=srite(tag,&buffer[pos],(int)pars->block)){
                         fprintf(stderr,"Trouble writing tag=%s file_pos=%g buffer=%g\n",
                           tag,(double)seek_it,(double)pos);
                          return 1;
                      }
                      pos+=pars->block;
                      
                    }
                    else{
                         
                      if(pars->nit==-1){
                        if((int)pars->block!=sreed(tag,&buffer[pos],(int)pars->block)){
                         fprintf(stderr,"Trouble reading tag=%s file_pos=%g buffer=%g\n",
                           tag,(double)seek_it,(double)pos);
                         return 1;
                        }
                        pos+=pars->block;
                      }
                      else{
                                    assert(BUFSIZE>=pars->block);
                        if((int)pars->block!=sreed(tag,wind_buf,(int)pars->block)){
                         fprintf(stderr,"Trouble reading tag=%s file_pos=%g buffer=%g\n",
                           tag,(double)seek_it,(double)pos);
                         return 1;
                        }
    
                       for(i=0; i < pars->nit; i++){
               
                         /*fprintf(stderr,"do a copy to from size %d %d %d\n",
           pos+pars->bit*i,(pars->fit+pars->jit*i)*pars->bit ,pars->bit);*/
                          memcpy((void*)(buffer+pos+pars->bit*i),
                            (const void*)(wind_buf+(pars->fit+pars->jit*i)*pars->bit),
                              pars->bit);
                        }
                         pos+=pars->bit*pars->nit;
                         
                     }
                    }
                  }
                } 
              }   
            }      
          }
         }
       }
     }
   }
   return 0;
  }




void compress_it(const int ndim,const int *ncube,const int *nwind, const int *fwind,const int *jwind,
const int esize, wind *comp, const int type){
  
  wind big;
  long long big_int=(long long)1024*(long long)1024*(long long)1024*2-1;
  int idim=0;
  int found=0;
  int first=0;
  long long block=esize;
  int i;
  int nd;
  int all=0;
  
  /*check axes get rid of axes length 1*/
  long long myalloc=esize;
  for(i=0; i < ndim; i++){
    myalloc=myalloc*nwind[i];
    if(fwind[i] < 0  || jwind[i] <1 || (fwind[i]+1)+jwind[i]*(nwind[i]-1) >ncube[i])
       seperr("Illegal window axis=%d ncube=%d nw=%d fw=%d jw=%d \n",
         i,ncube[i],nwind[i],fwind[i],jwind[i]);
    big.nc[idim]=ncube[i]; big.nw[idim]=nwind[i];
    
    big.jw[idim]=jwind[i]; big.fw[idim]=fwind[i];
    if(ncube[i] >1) idim++;
  }
  nd=idim;
  idim=0;
  
  i=0;
  /*Figure out biggest block we can read*/
  while(found==0 && i < nd){
      if(big.nc[i]==big.nw[i]){
         block=block*((long long)big.nw[i]);
         if(block > big_int) {
            found=1; comp->block=block/((long long)big.nw[i]);
         }
      }
     else{ found=1; comp->block=block;}
     i++;
  }

  if(i>1){
   for(idim=0; idim < nd -i+1; idim++){
     big.nc[idim]=big.nc[idim+i-1];
     big.nw[idim]=big.nw[idim+i-1];
     big.fw[idim]=big.fw[idim+i-1];
     big.jw[idim]=big.jw[idim+i-1];
   }
  }

  nd=nd-i+1;
  comp->nit=-1;
  /*A special case when reading a portion of the first non-1 axis */
  if(type==0 && block < 100 && (long long) block*(long long)big.nc[0] <myalloc && big.nw[0] > 5){
     comp->jit=big.jw[0]; comp->nit=big.nw[0]; comp->fit=big.fw[0];
     comp->bit=block; block=block*big.nc[0];
     nd--;
    for(idim=0; idim < nd ; idim++){
     big.nc[idim]=big.nc[idim+1];
     big.nw[idim]=big.nw[idim+1];
     big.fw[idim]=big.fw[idim+1];
     big.jw[idim]=big.jw[idim+1];
    
   }
 }

 comp->block=block;
 
 for(i=0; i < nd; i++){
   comp->bs[i]=block;
   comp->nc[i]=big.nc[i];
   comp->nw[i]=big.nw[i];
   comp->jw[i]=big.jw[i];
   comp->fw[i]=big.fw[i];
   block=block*((long long)big.nc[i]);
  
 }
 if(found==0) nd=0;
 for( i=nd; i < 9;i++){
   comp->nc[i]=comp->nw[i]=comp->jw[i]=1;
   comp->bs[i]=block;
   comp->fw[i]=0;
 }
   
}    
     
  
  
  
  
  
 
