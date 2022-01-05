/*$

=head1 NAME
   sep_put_grid_window  - write a window of the grid

=head1 SYNOPSIS 
C<int sep_put_grid_window(tag_history,n_dim_grid,n_grid,n_wind,f_wind,j_wind,header_numbers)>

=head1 INPUT PARAMETER

=over 4


=item   char* - tag_history      

        tag of History File

=item   int* - n_dim_grid	  

        Number of Dimensions in the Grid Cube

=item   int* - n_grid		  

        vector of length n_dim_grid-1 length of the axes in the Grid Cube

=item   int* - n_wind            

        vector of length n_dim_grid axes length after windowing

=item   int* - f_wind            

        vector of length n_dim_grid-1 index of first elements (0<= o < n)

=item   int* - j_wind            

        vector of length n_dim_grid-1 sampling rate along axes (1<= j < n)

=item   int* - header_numbers  

        array containing Header Record Numbes of dimensions
        n_wind[0],n_wind[1],..n_wind[n_dim_grid-2]

=back

=head1 RETURN VALUE

 -2=if the values in n_wind, f_wind, j_wind are incorrect

 -1=if fails for other reasons

 0=if successful

 +1=if tag_history is a Sep77 History File but has no GFF file

=head1 DESCRIPTION

   Writes gridding information.

=head1 COMMENTS

   Thes routines are based on the sep_window routine that it is the kernel
   of the new sreed_window and srite_window routines.
   ATTENTION: Since the first dimension of the grid is the first axis
              of the data record only the parameter for dimensions >=2
              are specified in n_wind,f_wind,d_wind.


=head1 SEE ALSO

L<sep_get_grid_window>

=head1 LIBRARY

B<sep3d>

=cut
>*/




/*$

=head1 NAME

   sep_get_grid_window  - read a window of the grid

=head1 SYNOPSIS 

C<int sep_get_grid_window(tag_history,n_dim_grid,n_grid,n_wind,f_wind,j_wind,header_numbers)>

=head1 INPUT PARAMETER

=over 4


=item   char* - tag_history      

                tag of History File

=item   int* - n_dim_grid	  

        Number of Dimensions in the Grid Cube

=item   int* - n_grid		  

        vector of length n_dim_grid-1 length of the axes in the Grid Cube

=item   int* - n_wind            

        vector of length n_dim_grid axes length after windowing

=item   int* - f_wind            

        vector of length n_dim_grid-1 index of first elements (0<= o < n)

=item   int* - j_wind            

        vector of length n_dim_grid-1 sampling rate along axes (1<= j < n)

=back

=head1 OUTPUT PARAMETERS

=over 4

=item   int* - header_numbers  

        array containing Header Record Numbers of dimensions
        n_wind[0],n_wind[1],..n_wind[n_dim_grid-2]


=back

=head1 RETURN VALUE

 -2= if the values in n_wind, f_wind, j_wind are incorrect

 -1= if fails for other reasons

 0= if successful

 +1= if tag_history is a Sep77 History File but has no GFF file

=head1 DESCRIPTION

   Read gridding information.

=head1 COMMENTS

   Thes routines are based on the sep_window routine that it is the kernel
   of the new sreed_window and srite_window routines.
   ATTENTION: Since the first dimension of the grid is the first axis
              of the data record only the parameter for dimensions >=2
              are specified in n_wind,f_wind,d_wind.


=head1 SEE ALSO

L<sep_put_grid_window>

=head1 LIBRARY

B<sep3d>

=cut
>*/









/*


KEYWORDS
   header keys

SEE ALSO
   sep3d

AUTHOR
    Biondo Biondi , September 1995

*/
#include <string.h>
#include "sep3d.h"
#include "streamlist.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_get_grid_window(char *tag_history, int *n_dim_grid, int *n_grid, 
		int *n_wind,int *f_wind,int *j_wind, int *header_numbers)
_XFUNCPROTOEND
#else 
int sep_get_grid_window(tag_history, n_dim_grid, n_grid, n_wind, f_wind, j_wind, header_numbers)
char *tag_history;
int  *n_dim_grid;
int  *n_grid;
int  *n_wind, *f_wind, *j_wind;
int *header_numbers;
#endif

{
streaminf *info;
char *tag_grid[1];
int ierr,mode;
int new_n_dim_grid;

/* adapt parameters to Grid Format File */
new_n_dim_grid=*n_dim_grid-1;
/* Get tag_header and check for returning errors */
ierr=sep_get_grid_format_tag(tag_history, tag_grid);

if(ierr==0) {

/* Check if grid file exist*/
   info = tag_info(*tag_grid, TAG_IN);  /* get info on this tag */
   if((info->dataname) == 0 || !strncmp("-1",(info->dataname),2)) {

      mode=0;
      ierr=sep_window(mode,*tag_grid, &new_n_dim_grid, n_grid, n_wind, f_wind, j_wind, 4, header_numbers);


	}
   	else {

    mode=-1;
    ierr=sep_window(mode,*tag_grid, &new_n_dim_grid, n_grid, n_wind, f_wind, j_wind, 4, header_numbers);
    

	     }
    free(*tag_grid);
	}
return (ierr);
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_put_grid_window(char *tag_history, int *n_dim_grid, int *n_grid, 
		int *n_wind,int *f_wind,int *j_wind, int *header_numbers)
_XFUNCPROTOEND
#else 
int sep_put_grid_window(tag_history, n_dim_grid, n_grid, n_wind, f_wind, j_wind, header_numbers)
char *tag_history;
int  *n_dim_grid;
int  *n_grid;
int  *n_wind, *f_wind, *j_wind;
int *header_numbers;
#endif

{
    char *tag_grid[1];
    int ierr,mode;
    int new_n_dim_grid;


/* adapt parameters to Grid Format File */
    new_n_dim_grid=*n_dim_grid-1;


    /* Get tag_header and check for returning errors */
    ierr=sep_get_grid_format_tag(tag_history, tag_grid);
    if(ierr!=0) {
        return ierr;
	}


    mode=1;
    ierr=sep_window(mode,*tag_grid, &new_n_dim_grid, n_grid, n_wind, f_wind, j_wind, 4, header_numbers);


    free(*tag_grid);

    return (ierr);
    
}
