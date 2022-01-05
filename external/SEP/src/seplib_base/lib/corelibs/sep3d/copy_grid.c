/*$

=head1 NAME

sep_copy_grid - copy a grid  from one tag to another 

=head1 SYNOPSIS

C<ierr=sep_copy_grid( intag, outtag)>

=head1 INPUT PARAMETERS

=over 4

=item intag - char* 

      input tag

=item outtag - char* 

      output tag

=back

=head1 RETURN VALUES

 0= if sucessfully copies the grid

 1= if the grid doesn't exist in the intag

 -1= fails for some other reason

=head1 DESCRIPTION

When the grid exists, copies contents from intag to outtag

=head1 LIBRARY

B<sep3d>

=cut



>*/ 
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Thu Dec 18 11:09:05 PST 1997

Purpose: 

*/	 

#include <sep3d.h> 

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_copy_grid(char *intag, char *outtag)
_XFUNCPROTOEND
#else
int sep_copy_grid(intag,outtag) 
char *intag;
char *outtag;
#endif 
{ 
int read,n,ndim,ntot,i,readin;
float o,d;
char label[128];

char *in_grid_tag[1],*out_grid_tag[1],*buffer;

if(0!=sep_get_grid_format_tag(intag,in_grid_tag)) return(1);
if(0!=sep_get_grid_format_tag(outtag,out_grid_tag)) return(-1);

if(0!=sep_get_number_grid_axes(intag,&ndim)) return(-1);
	
ntot=1;
for (i=1; i<= ndim; i++){
	if(0!=sep_get_grid_axis_par(intag,&i,&n,&o,&d,label)) return(-1);
	ntot=ntot*n;
}

copy_history(*in_grid_tag,*out_grid_tag);
ntot=ntot*4;

buffer=(char*) malloc(sizeof(char) * 65536);

read=0;
while(read < ntot ){
	readin=MIN(ntot-read,65536);
	if(readin!=sreed(*in_grid_tag,buffer,readin)) 
		seperr("sep_copy_grid(): trouble reading in grid values \n");
	if(readin!=srite(*out_grid_tag,buffer,readin))
		seperr("sep_copy_grid(): trouble writing out grid values \n");
	read=read+readin;
}
return (0);
} 


