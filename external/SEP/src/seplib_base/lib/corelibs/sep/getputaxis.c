/*$

=head1 NAME

sep_get_data_axis_par - grab data's n,o,d and label


=head1 SYNOPSIS

C<ierr=sep_get_data_axis_par(tag_history,i_axis,n,o,d,label)>

=head1 INPUT  PARAMETERS

=over 4

=item char* - tag_history  

      Tag of History File

=item int*  -  i_axis      

       Axis number of the parameters to read     

=back

=head1 OUTPUT PARAMETERS

=over 4

=item int* - n      

      value of n-i (i = i_axis)

=item float* - o    
       
      value of o-i (i = i_axis)

=item float* - d    

      value of d-i (i = i_axis)

=item char* - label 

      value of label-i (i = i_axis)

=back

=head1 RETURN VALUES

 -1 = if fails for other reasons

 0 = if successful

 +1 = if tag_history is not a Sep History File

=head1 DESCRIPTION

        Get the values that describe the Data Axes from the History
        File. 

=head1 COMMENTS

        Defaults the others (o=0 d=1 label="label-i axis").

=head1 SEE ALSO

L<sep_get_header_axis_par>, L<sep_get_grid_axis_par>, L<sep_put_data_axis_par>


=head1 LIBRARY

B<sep3d>

=cut


>*/

#include <sepConfig.h>
#include <stdio.h>
#include <sep_main_external.h>
#include "sep3d.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_get_data_axis_par(const char *tag_history, int *i_axis, 
                             int *n, float *o, float *d, char *label)
_XFUNCPROTOEND
#else
int sep_get_data_axis_par(tag_history,i_axis,n,o,d,label)
char *tag_history;
int *i_axis, *n;
float *o,  *d;
char *label;
#endif
{
int ierr;
char temp_ch[128];

sprintf(temp_ch,"n%d",*i_axis);
ierr=auxpar(temp_ch,"d",n,tag_history);
if(ierr==0) *n=1;
else if (ierr<0) return ierr;

sprintf(temp_ch,"o%d",*i_axis);
ierr=auxpar(temp_ch,"f",o,tag_history);
if(ierr==0) *o=0;
else if (ierr<0) return ierr;

sprintf(temp_ch,"d%d",*i_axis);
ierr=auxpar(temp_ch,"f",d,tag_history);
if(ierr==0) *d=1;
else if (ierr<0) return ierr;

sprintf(temp_ch,"label%d",*i_axis);
ierr=auxpar(temp_ch,"s",label,tag_history);
if(ierr==0) label[0]='\0';
else if (ierr<0) return ierr;


return 0;
}

/*$

=head1 NAME

sep_put_data_axis_par  - put n,o,d, label to output tag

=head1 SYNOPSIS

C<int sep_put_data_axis_par(tag_history,i_axis,n,o,d,label)>

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

   Put the values that describe the i axis of the Data
   Coordinate System into the Header Format File pointed 
   by the History File corresponding to tag history.


=head1 SEE ALSO

L<sep_put_header_axis_par>, L<sep_put_grid_axis_par>, L<sep_get_data_axis_par>


=head1 LIBRARY

B<sep3d>

=cut


>*/

/*
KEYWORDS
   data 

SEE ALSO
   sep3d

AUTHOR
   SEP , July ... 1995
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_put_data_axis_par(const char *tag_history, int *i_axis,
			   int *n, float *o, float *d, const char *label)
_XFUNCPROTOEND
#else 
int sep_put_data_axis_par(tag_history, i_axis, n, o, d, label)
const char *tag_history;
int *i_axis;
int *n;
float *o;
float *d;
const char *label;
#endif
{
    char ni[32], oi[32], di[32], labeli[32];


		if(0!=auxputhead((char *) tag_history,"\t \tn%d=%d o%d=%g d%d=%g label%d=\"%s\" \n",
      *i_axis,*n,*i_axis,*o,*i_axis,*d,*i_axis,label))
				seperr("trouble writing data axis for tag %s \n",tag_history);
    
    return 0;	
}
