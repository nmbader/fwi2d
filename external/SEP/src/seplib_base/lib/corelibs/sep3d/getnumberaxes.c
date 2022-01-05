#include <stdio.h>
#include <stdlib.h>
#include <sep_pars_external.h>
#include <sep3d.h>

int sep_get_number_header_axes (const char *tag_history, int *n_axis)
{
    char *tag_header[1];
    char ni[8];
    int ierr=0, exist=1, i_axis=1, j_axis=1;
    int n;
    /* Check for returning errors */
    ierr=sep_get_header_format_tag(tag_history, tag_header);
    if(ierr!=0) {
        return ierr;
    }

    sprintf(ni, "n%d", i_axis+1);
    ierr=auxpar(ni, "d", &n, *tag_header);


    /* Auxpar i_axis from header format file */
    while (exist) {
        if(n !=1) j_axis=i_axis;
	i_axis++;
        sprintf(ni, "n%d", i_axis);
        exist=auxpar(ni, "d", &n, *tag_header);
    }
    *n_axis=j_axis;

    free(*tag_header);

    return 0;	
}
/*$

=head1 NAME

sep_get_number_header_axes -get number of header axes
 
=head1 SYNOPSIS

int sep_get_number_header_axes (tag_history,n_axis)

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

L<sep_get_number_data_axes>, L<sep_get_number_grid_axes>

=head1 LIBRARY

B<sep3d>

=cut


>*/

/*
SEE ALSO
   sep3d

AUTHOR
   SEP , July ... 1995
*/

/*$
=head1 NAME

sep_get_number_grid_axes -get number of grid axes
 
=head1 SYNOPSIS

int sep_get_number_grid_axes (tag_history,n_axis)

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

Get number of axes in grid file

=head1 SEE ALSO

L<sep_get_number_header_axes>, L<sep_get_number_data_axes>

=head1 LIBRARY

B<sep3d>

=cut


>*/

/*
SEE ALSO
   sep3d

AUTHOR
   SEP , July ... 1995
*/
#include "sep3d.h"

int sep_get_number_grid_axes (const char *tag_history, int *n_axis)
{
    char *tag_grid[1];
    char ni[8];
    int ierr=0, exist=1, i_axis=1, j_axis=1;
    int n;


    /* Check for returning errors */
    ierr=sep_get_grid_format_tag(tag_history, tag_grid);
    if(ierr!=0) {
        return ierr;
    }

    sprintf(ni, "n%d_grid", i_axis+1);
    ierr=auxpar(ni, "d", &n, *tag_grid);

    /* Auxpar i_axis from grid format file */
    while (exist) {
        if(n !=1) j_axis=i_axis;
	i_axis++;
        sprintf(ni, "n%d_grid", i_axis);
        exist=auxpar(ni, "d", &n, *tag_grid);
    }
    *n_axis=j_axis;

    free(*tag_grid);

    return 0;	
}

