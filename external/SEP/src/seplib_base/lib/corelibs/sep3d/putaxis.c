/*$

=head1 NAME

sep_put_header_axis_par  - put n,o,d, label to output tag

=head1 SYNOPSIS

C<int sep_put_header_axis_par(tag_history,i_axis,n,o,d,label)>

=head1 INPUT PARAMETER

=over 4

=item   char* -  tag_history   

        tag of History File

=item   int* -    i_axis       

        Axis number of parameters to read

=item   int* -    n            

        value of nh-i (i = i_axis)

=item   float* -  o            

        value of oh-i (i = i_axis)

=item   float* -  d            

        value of dh-i (i = i_axis)

=item   char*l -  abel         

        value of labelh-i (i = i_axis)

=back

=head1 RETURN VALUE

 -1= if fails for other reasons

 0= if successful

 +1= if tag_history is an Sep77 History File

=head1 DESCRIPTION 

   Put the values that describe the i axis of the Header
   Coordinate System into the Header Format File pointed 
   by the History File corresponding to tag history.


=head1 SEE ALSO

L<sep_put_data_axis_par>, L<sep_put_grid_axis_par>, L<sep_get_header_axis_par>


=head1 LIBRARY

B<sep3d>

=cut


>*/
/*

KEYWORDS
   header 

SEE ALSO
   sep3d

AUTHOR
   SEP , July ... 1995
*/
#include "sep3d.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_put_header_axis_par (const char *tag_history, int *i_axis, int *n, float *o,                           float *d, const char *label)
_XFUNCPROTOEND
#else 
int sep_put_header_axis_par (tag_history, i_axis, n, o, d, label)
char *tag_history;
int *i_axis;
int *n;
float *o;
float *d;
const char *label;
#endif
{
/*    streaminf *info;*/
    char *tag_header[1];
    char ni[32], oi[32], di[32], labeli[32];
    int ierr=0;



    if(*i_axis == 1) {
	seperr("sep_put_header_axis_par i_axis == 1 is not legal --- Use sep_put_number_keys");
	}

    /* Check for returning errors */
    ierr=sep_get_header_format_tag(tag_history, tag_header);
    if(ierr!=0) {
        return ierr;
    }

		if(0!=auxputhead(*tag_header,"\t \tn%d=%d o%d=%g d%d=%g label%d=\"%s\"\n",
      *i_axis,*n,*i_axis,*o,*i_axis,*d,*i_axis,label))
				seperr("trouble writing data axis for tag %s \n",tag_history);

    free(*tag_header);
    return 0;	
}
/*$

=head1 NAME

sep_put_grid_axis_par  - put n,o,d, label to output tag

=head1 SYNOPSIS

C<int sep_put_grid_axis_par(tag_history,i_axis,n,o,d,label)>

=head1 INPUT PARAMETER

=over 4

=item   char* -  tag_history   

        tag of History File

=item   int* -    i_axis       

        Axis number of parameters to read

=item   int* -    n            

        value of nh-i (i = i_axis)

=item   float* -  o            

        value of oh-i (i = i_axis)

=item   float* -  d            

        value of dh-i (i = i_axis)

=item   char*l -  abel         

        value of labelh-i (i = i_axis)

=back

=head1 RETURN VALUE

 -1= if fails for other reasons

 0= if successful

 +1= if tag_history is an Sep77 History File

=head1 DESCRIPTION 

   Put the values that describe the i axis of the Grid
   Coordinate System into the Header Format File pointed 
   by the History File corresponding to tag history.


=head1 SEE ALSO

L<sep_put_data_axis_par>, L<sep_put_header_axis_par>, L<sep_get_grid_axis_par>


=head1 LIBRARY

B<sep3d>

=cut


>*/
/*
KEYWORDS
   Binning Coordinate System

SEE ALSO
   sep3d

AUTHOR
   SEP , July ... 1995
*/
#include "sep3d.h"

int sep_put_grid_axis_par (const char *tag_history, int *i_axis, int *n, float *o,                           float *d, const char *label)
_XFUNCPROTOEND
{
    char *tag_grid[1];
    char ni[32], oi[32], di[32], labeli[32];
    int ierr=0;

    if(*i_axis == 1) {
	seperr("sep_put_grid_axis_par i_axis == 1 is not legal");
	}

    /* Check for returning errors */
    ierr=sep_get_grid_format_tag(tag_history, tag_grid);
    if(ierr!=0) {
        return ierr;
    }
		if(0!=auxputhead(*tag_grid,"\t \tn%d_grid=%d o%d_grid=%g d%d_grid=%g label%d_grid=\"%s\"\n",
      *i_axis,*n,*i_axis,*o,*i_axis,*d,*i_axis,label))
				seperr("trouble writing data axis for tag %s \n",tag_history);
    

    free(*tag_grid);

    
    return 0;	
}

