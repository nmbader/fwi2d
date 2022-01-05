/*$
=head1 NAME

sep_get_number_data_axes - get number of data axes
 
=head1 SYNOPSIS

C<ierr= sep_get_number_data_axes (tag_history,n_axis)>

=head1 INPUT PARAMETER

=over 4

=item    char*-tag_history     

         tag of History File

=back

=head1 OUTPUT PARAMETER

=over 4

=item   int*-n_axis           

         Number of axes 

=back

=head1 RETURN VALUE

 -1= if fails for other reasons

 0= if successful

 +2= if tag_history is an Sep77 History File

=head1 DESCRIPTION 

	Get number of axes in history file

=head1 SEE ALSO

L<sep_get_number_header_axes>, L<sep_get_number_grid_axes>

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
int sep_get_number_data_axes (const char *tag_history, int *n_axis)
_XFUNCPROTOEND
#else 
int sep_get_number_data_axes (tag_history, n_axis)
char *tag_history;
int *n_axis;
#endif
{
    char ni[8];
    int  exist=1, i_axis=0, j_axis=0;
    int n;


    sprintf(ni, "n%d", i_axis+1);
    if(1!=auxpar(ni, "d", &n, tag_history))
			seperr("invalid seplib3d dataset,%s, n1 must exist and by a single integer \n",tag_history);


    /* Auxpar i_axis from History File  */
    while (exist) {
        if(n !=1) j_axis=i_axis;
	i_axis++;
        sprintf(ni, "n%d", i_axis);
        exist=auxpar(ni, "d", &n, tag_history);
    }
    if(j_axis==0) j_axis=1;
    *n_axis=j_axis;



    return 0;	
}
