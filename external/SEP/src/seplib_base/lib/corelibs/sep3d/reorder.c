#include <sep3d.h>
#define MAX_TRACES 10000

/*$

=head1 NAME

sep_reorder_data - reorder SEPlib dataset

=head1 SYNOPSIS

C<int sep_reorder_data(char *tagin , char *tagout, int n2h, int tsize, int *order)>

=head1 DESCRIPTION
reorder traces

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

=back


=head1 LIBRARY

B<sep3d>

=cut


*/

/*

Author: Robert G. Clapp
Modified: Jul 1999 Buffer output 
>*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_reorder_data(char *tagin, char *tagout, int n2h, int tsize, int *order)
_XFUNCPROTOEND
#else
int sep_reorder_data(tagin,tagout,n2h,tsize,order)
int  n2h,tsize,*order;
char *tagin,*tagout;
#endif
{
int block_count,beg_block,size,new,block_size;
int position,type,offset,ierr,out_count;
char *buffer,*out_pt;

block_size= tsize; 
buffer=(char*) malloc( (long long) block_size * (long long) MAX_TRACES * sizeof(char ) ) ;
out_count=0;  out_pt=buffer;
size=1; 
for(beg_block=0;beg_block<n2h;beg_block++){
  for (block_count=beg_block; block_count<(n2h-2); block_count++){
		if(order[block_count]!=(order[block_count+1]-1)/*traces not consecutive*/
		|| out_count+block_count	-beg_block+1 >= MAX_TRACES)/*limit of buffer*/
			break;
  }
	size=block_count-beg_block+1; /*the number of traces we are processing*/
  new=(order[beg_block]-1);     /*c location of the first trace*/
	ierr=sseek_block(tagin,new,block_size,0);
	if(new!=ierr) /*seek to the block begining*/
		seperr("trouble seeking to the begining of block, tag %s, ntr=%d, size=%d err=%d \n", tagin,new,block_size,ierr);
  ierr=sreed(tagin,out_pt,size*block_size);
  if(block_size*size!=ierr)
		seperr("trouble reading in data (err=%d) starting_trace=%d n_bytes=%d \n",
      ierr,new,size*block_size);


	out_count+=size;
	if(out_count==MAX_TRACES || block_count==n2h-1){ /*We need to write out */
  	if(out_count*block_size!=srite(tagout,buffer,block_size*out_count))
			seperr("trouble writing out data to tag %s \n",tagout);
		out_count=0;
		out_pt=buffer;
	} 
	else  out_pt=(char*)(out_pt+block_size*size); /*we need to update out_pt */

  beg_block=block_count; /*reset our location */
}

free(buffer);
return(0);
}
