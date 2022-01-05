

/*$

=head1 NAME

sep_get_header_axis_par - grab headers's n,o,d and label


=head1 SYNOPSIS

C<ierr=sep_get_header_axis_par(tag_history,i_axis,n,o,d,label)>

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

Get the values that describe the Header Axes from the History
File. 

=head1 COMMENTS

Defaults the others (o=0 d=1 label="label-i axis").

=head1 SEE ALSO

sep_get_data_axis_par, sep_get_grid_axis_par, sep_put_header_axis_par


=head1 LIBRARY

B<sep3d>

=cut
>*/

#include <stdio.h>
#include <stdlib.h>
#include <sep_pars_external.h>
#include <sep3d.h>

int sep_get_header_axis_par(const char *tag_history, int *i_axis, 
                             int *n, float *o, float *d, char *label)
{
int ierr;
char temp_ch[128];
char *tag_header[1];


    ierr=sep_get_header_format_tag(tag_history, tag_header);

    if(ierr!=0) {
        return ierr;
    }


sprintf(temp_ch,"n%d",*i_axis);
ierr=auxpar(temp_ch,"d",n,*tag_header);
if(ierr==0) *n=1;
else if (ierr<0) return ierr;

sprintf(temp_ch,"o%d",*i_axis);
ierr=auxpar(temp_ch,"f",o,*tag_header);
if(ierr==0) *o=0;
else if (ierr<0) return ierr;

sprintf(temp_ch,"d%d",*i_axis);
ierr=auxpar(temp_ch,"f",d,*tag_header);
if(ierr==0) *d=1;
else if (ierr<0) return ierr;

sprintf(temp_ch,"label%d",*i_axis);
ierr=auxpar(temp_ch,"s",label,*tag_header);
if(ierr==0) label[0]='\0';
else if (ierr<0) return ierr;


free(*tag_header);
return 0;


}
/*$
=head1 NAME

sep_get_grid_axis_par - grab data's n,o,d and label


=head1 SYNOPSIS

C<ierr=sep_get_grid_axis_par(tag_history,i_axis,n,o,d,label)>

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

Get the values that describe the GRid Axes from the History
File. 

=head1 COMMENTS
        Defaults the others (o=0 d=1 label="label-i axis").

=head1 SEE ALSO

sep_get_data_axis_par, sep_get_header_axis_par, sep_put_grid_axis_par


=head1 LIBRARY

B<sep3d>

=cut

>*/
 

int sep_get_grid_axis_par(const char *tag_history, int *i_axis, 
                             int *n, float *o, float *d, char *label)
{
int ierr;
char temp_ch[128];
char *tag_grid[1];


    ierr=sep_get_grid_format_tag(tag_history, tag_grid);

    if(ierr!=0) {
        return ierr;
    }


sprintf(temp_ch,"n%d_grid",*i_axis);
ierr=auxpar(temp_ch,"d",n,*tag_grid);
if(ierr==0) *n=1;
else if (ierr<0) return ierr;

sprintf(temp_ch,"o%d_grid",*i_axis);
ierr=auxpar(temp_ch,"f",o,*tag_grid);
if(ierr==0) *o=0;
else if (ierr<0) return ierr;

sprintf(temp_ch,"d%d_grid",*i_axis);
ierr=auxpar(temp_ch,"f",d,*tag_grid);
if(ierr==0) *d=1;
else if (ierr<0) return ierr;


sprintf(temp_ch,"label%d_grid",*i_axis);
ierr=auxpar(temp_ch,"s",label,*tag_grid);
if(ierr==0) label[0]='\0';
else if (ierr<0) return ierr;


free(*tag_grid);
return 0;


}

