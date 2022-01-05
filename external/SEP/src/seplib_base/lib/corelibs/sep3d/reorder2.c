#include <sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include <sep3d.h>
int partition(int  **A, int first, int last) ;
void quicksort(int  **A, int first, int last) ;
int read_block(char *tagin,char *buffer,int iold,int tsize,int nread,int new);

/*$

=head1 NAME

sep_reorder_data_fast - reorder SEPlib dataset

=head1 SYNOPSIS

C<int sep_reorder_data_fast(char *tagin , char *tagout, int n2h, int tsize, int *order,int megabytes)>

=head1 DESCRIPTION
reorder traces (allocates something n2h*tsize  so be careful)

=head1 INPUT PARAMETERS

=over 4

=item tagin - char*  

      tag of input

=item tagout- char*  

      tag of output

=item n2h   - int    

      number of traces

=item tsize  - int    

      length of traces

=item order - int*   

      order of traces

=item megabyte - int

      maximum number of megabytes to use

=back


=head1 LIBRARY

B<sep3d>

=cut


*/

/*

Author: Robert G. Clapp
Modified: Jul 1999 Buffer output 
Modified: Sep 2003 Converted reorder_data
>*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_reorder_data_fast(char *tagin, char *tagout, int n2h, int tsize, int *order,int megabyte)
_XFUNCPROTOEND
#else
int sep_reorder_data_fast(tagin,tagout,n2h,tsize,order,megabyte)
int  n2h,tsize,*order,megabyte;
char *tagin,*tagout;
#endif
{
int **list;
int i,new,iloc,nread,ierr,iold;
char *buf_in,*buf_out;
int nblock,ndone,npart;

nblock=megabyte*1000*500/tsize;
if(nblock > n2h) nblock=n2h;



buf_in=(char*) malloc(sizeof(char)*nblock*tsize);
buf_out=(char*) malloc(sizeof(char)*nblock*tsize);
list=(int**) malloc(sizeof(int*)*nblock);

/* make a list[nh2][2] where 0 is the trace number, 1 is the place in array*/
for(i=0; i < nblock; i++) list[i]=(int*)malloc(sizeof(int)*2);

ndone=0;
while(ndone< n2h){
  npart=nblock; if(n2h - ndone < npart) npart=n2h-ndone;

  for(i=0; i < npart; i++){ 
   list[i][0]=order[i+ndone]; list[i][1]=i; 
   }

  /*sort the list according to trace number*/
  quicksort(list,0,npart-1);

  new=list[0][0]; /*initialize the first trace number we want to read*/
  i=1; nread=1; 
  iold=0; 

  /*loop over the traces*/
  while(i < npart){ 
    if(list[i][0] != new +nread){ /* if this trace isn't in sequence*/
      /* setup for a new read   */
      if(0!=read_block(tagin,buf_in,iold,tsize,nread,new))
       seperr("trouble reading in block \n");
      iold+=nread; nread=1; new=list[i][0];
    }
    /*traces in order*/
    else nread++;
    i++;
 }
  /*do the last block*/
 if(0!=read_block(tagin,buf_in,iold,tsize,nread,new))
    seperr("trouble reading in block \n");

  for(i=0; i < npart; i++){
/*    memcpy((void*)(buf_out+i*tsize),(const void*)(buf_in+tsize*list[i][1]),tsize);*/
    memcpy((void*)(buf_out+list[i][1]*tsize),(const void*)(buf_in+i*tsize),tsize);
  }
  if(npart*tsize!=srite(tagout,buf_out,npart*tsize))
    seperr("trouble writing out %d bytes to %s \n",npart*tsize,tagout);
 ndone+=npart; 
}

for(i=0; i< nblock; i++) free(list[i]);
free(buf_in); free(buf_out); free(list);
return(0);
}



int read_block(char *tagin,char *buffer,int iold,int tsize,int nread,int new){
int ierr;


    /*seek to the block begining*/
	  ierr=sseek_block(tagin,new-1,tsize,0);
	  if(new-1!=ierr)  /*give an error if we couldn't seek */
      return(sepwarn(-1,
        "trouble seeking to trace %d for trace read %d \n",new));

    /*read in the trace block*/
    ierr=sreed(tagin,(buffer+iold*tsize),tsize*nread);
    if(tsize*nread!=ierr)
      return(sepwarn(-1,"trouble reading %d bytes from %s read=%d \n",
        tsize*nread,tagin,ierr));
    
   return(0);
}

/*from 

http://ironbark.bendigo.latrobe.edu.au/courses/bcomp/c103/sem296/lectures/Lecture2.html
*/

void quicksort(int  **A, int first, int last) {
    int pivindx; /* index of the element separating the two sub-arrays*/
    if (last > first) {  /* More than one element to be sorted?*/
        pivindx = partition(A, first, last);
        quicksort(A, first, pivindx - 1);
        quicksort(A, pivindx + 1, last);
    }
}

int partition(int  **A, int first, int last) {
    int pivindx, top, i; 
    int pivot[2];

    /* Choose a pivot: select a random index between first and last.*/
    i = rand() % (last - first + 1) + first;

    /* Put the pivot first, remember pivot, initialise ready for loop.*/
    pivot[0] = A[i][0]; pivot[1] = A[i][1];                   
    A[i][0] = A[first][0]; A[i][1] = A[first][1];
    A[first][0] = pivot[0];A[first][1] = pivot[1]; /* pivot now first */

    pivindx = first;
    top = last;                       
    while (top > pivindx) {         /* Still unknown elements */
        /* top indicates the highest unknown element */
        if (A[top][0] >= pivot[0]) {
            top--;              /* where it belongs, count as >=*/
        } else {
            A[pivindx][0] = A[top][0];      /* shift down */
            A[pivindx][1] = A[top][1];      /* shift down */
            A[top][0] = A[pivindx + 1][0];/*shift displaced element up*/
            A[top][1] = A[pivindx + 1][1];/*shift displaced element up*/
            A[pivindx + 1][0] = pivot[0];   /* Put pivot back */
            A[pivindx + 1][1] = pivot[1];   /* Put pivot back */
            pivindx++;            /* Alter record of pivot location*/
        }
    }
    return pivindx;
}


