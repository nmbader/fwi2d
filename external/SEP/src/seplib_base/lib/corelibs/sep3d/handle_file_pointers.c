/*$

=head1 NAME

   copy_data_pointer -copy data pointer from one file to another

=head1 SYNOPSIS
 
   int sep_copy_data_pointer(tag_in,tag_out)

=head1 INPUT PARAMETER

=over 4

=item   char* - tag_in     

        tag of History File from which to copy data pointer

=item   char* - tag_out    

        tag of History File to which to copy data pointer

=back

=head1 RETURN VALUE

 0= if successful

 -1= if either tag_in or tag_out are not a SEP History File

=head1 DESCRIPTION
   Copy the Data Pointer from tag_in to tag_out.
   Useful to avoid duplications of Data File when the Header Files 
   are modified, but the Data is not.


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
   Biondo Biondi , November 1995

*/

#include "sep3d.h"
#include "streamlist.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_copy_data_pointer(char* tag_in, char* tag_out)
_XFUNCPROTOEND
#else
int sep_copy_data_pointer(tag_in, tag_out)
char *tag_in; 
char *tag_out;
#endif
{
    streaminf *info_in,*info_out;
 

    info_in = tag_info(tag_in, TAG_INQUIRE);  /* get info on this tag */
    info_out = tag_info(tag_out, TAG_INQUIRE);  /* get info on this tag */
    if((info_in == 0)||(info_out ==0)) return -1;  /* Not SEP History Files */

    sepstrput( info_out, "sets next: in", "s", info_in->dataname );

    return 0;
}

/*$

=head1 NAME

sep_set_no_headers - set no headers to an output tag
   
=head1 SYNOPSIS

#include <sep3d.h>

int sep_set_no_headers(char *tag_history)

=head1 INPUT PARAMETER

=over 4

=item   char *tag_history   :  

         tag of History File for which to null the gff pointer

=back


=head1 RETURN VALUE

 0 if successful

 -1 if either tag_history is not a SEP History File

=head1 DESCRIPTION

   Set the gff pointer to -1 in the History File tag_history

=head1 LIBRARY

B<sep3d>

=cut


*/
/*
AUTHOR
   Biondo Biondi; December 1995

*/

#include "sep3d.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_set_no_headers(char* tag_history)
_XFUNCPROTOEND
#else
int sep_set_no_headers(tag_history)
char *tag_history; 
#endif
{
    if(0!=auxputch("hff","s", "-1", tag_history))
			seperr("trouble writing hff = -1 to tag %s \n",tag_history);

    return 0;
}

/*$
=head1 NAME

sep_set_no_grid - set that the output tag does not have a grid

=head1 SYNOPSIS

#include <sep3d.h>

int sep_set_no_grid(char *tag_history)

=head1 INPUT PARAMETER

=over 4

=item  char *tag_history   
  
     tag of History File for which to null the gff pointer

=back


=head1 RETURN VALUE

 0 if successful

 -1 if either tag_history is not a SEP History File

=head1 DESCRIPTION

   Set the gff pointer to -1 in the History File tag_history

=head1 SEE ALSO

L<sep_set_no_headers>

=head1 LIBRARY

B<sep3d>

=cut

*/

#include "sep3d.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_set_no_grid(char* tag_history)
_XFUNCPROTOEND
#else
int sep_set_no_grid(tag_history)
char *tag_history; 
#endif
{
		if(0!= auxputch("gff","s", "-1", tag_history))
			seperr("trouble writing gff = -1 to tag %s \n",tag_history);
    return 0;
}

/*
NAME:
   sep_set_regular_grid

   
SYNOPSIS
   #include <sep3d.h>

   int sep_set_regular_grid(char *tag_history)

PARAMETER
   INPUT
   char *tag_history   :  tag of History File for which to null the gff pointer


RETURN VALUE
    0 if successful
   -1 if it fails for other reasons
   +1 if tag_history is neither a Sep77 nor a Sep3d History File
   +2 if tag_history is a Sep77 History File but not Sep3d
   +3 if tag_history is a Sep3d Hisory File but 
                   has no Grid Format File associated.

DESCRIPTION
   Set the in pointer to -1 in the Grid Format File associated with tag_history
KEYWORDS
   header 

SEE ALSO
   sep3d

AUTHOR
   Biondo Biondi; December 1995

*/

#include "sep3d.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_set_regular_grid(char* tag_history)
_XFUNCPROTOEND
#else
int sep_set_regular_grid(tag_history)
char *tag_history; 
#endif
{
    char *tag_grid[1];
    int ierr;

    ierr=sep_get_grid_format_tag(tag_history, tag_grid);
    if (ierr != 0) {
      return (ierr);
    }

    if(0!=auxputch("sets next: in","s", "-1", *tag_grid))
			seperr("trouble writing regular grid flag to tag %s\n",tag_history);

    return 0;
}


/*$


=head1 NAME

   sep_copy_gff - copy grid format file  from one SEP3d file to another

=head1 SYNOPSIS

   #include <sep3d.h>
   int sep_copy_gff(char *tag_in, char *tag_out)

=head1 PARAMETER

=over 4

=item    char *tag_in : 

         tag of input History File


=item    char *tag_out : 

         tag of output History File

=back

=head1 RETURN VALUE

 0 if successful

 +1 if tag_history is an Sep77 History File

 -1 if fails for other reasons

=head1 DESCRIPTION

Copy gff

=head1 SEE ALSO

L<sep_copy_hff>

=head1 LIBRARY

B<sep3d>

=cut

*/

/*

AUTHOR
   Bob , December .. 1995

*/
#include "sep3d.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_copy_gff(char *tag_in, char *tag_out)
_XFUNCPROTOEND
#else 
int sep_copy_gff(tag_in,tag_out)
char *tag_in,*tag_out;
#endif
{
    int i1,num_axis, naxis,ierr;
    char label[128];
    float oaxis,daxis;
    char *tag_grid_in[1];
    char *tag_grid_out[1];


	ierr=sep_get_number_grid_axes(tag_in,&num_axis);
	  if(ierr!=0) return ierr;

    /* Copy grid axis */
	for(i1=2;i1<=num_axis;i1++){
	ierr=sep_get_grid_axis_par(tag_in,&i1,&naxis,&oaxis,&daxis,label);
	  if(ierr!=0) return ierr;
	ierr=sep_put_grid_axis_par(tag_out,&i1,&naxis,&oaxis,&daxis,label);
	  if(ierr!=0) return ierr;
	}


	ierr = sep_get_grid_format_tag(tag_in, tag_grid_in);
	  if(ierr!=0) return ierr;
	ierr = sep_get_grid_format_tag(tag_out, tag_grid_out);
	  if(ierr!=0) return ierr;
	ierr=sep_copy_data_pointer(*tag_grid_in, *tag_grid_out);
	  if(ierr!=0) return ierr;



	 free(*tag_grid_in);
	 free(*tag_grid_out);
 
    return 0;
}




/*$

=head1 NAME

   sep_copy_hff - copy header format file  from one SEP3d file to another

=head1 SYNOPSIS

   #include <sep3d.h>

   int sep_copy_hff(char *tag_in, char *tag_out)

=head1 PARAMETER

=over 4

=item   char *tag_in 

        tag of input History File

=item   char *tag_out 

       tag of output History File

=back


=head1 RETURN VALUE

   +1 if tag_history is an Sep77 History File

    0 if successful

   -1 if fails for other reasons


=head1 DESCRIPTION

    copies hff


=head1 SEE ALSO

L<sep_copy_hff>

=head1 LIBRARY

B<sep3d>

=cut

AUTHOR
   Bob , December ... 1995

*/
#include "sep3d.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep_copy_hff(char *tag_in, char *tag_out)
_XFUNCPROTOEND
#else 
int sep_copy_hff(tag_in,tag_out)
char *tag_in,*tag_out;
#endif
{
    int i1,num_axis,nkeys, naxis,ierr;
    char label[128];
    float oaxis,daxis;
    char *tag_header_in[1];
    char *tag_header_out[1];


    /* Copy number of keys */
	ierr=sep_get_number_keys(tag_in,&nkeys);
	  if(ierr!=0) return ierr;
	ierr=sep_put_number_keys(tag_out,&nkeys);
	  if(ierr!=0) return ierr;

/*	Copy in*/

	ierr = sep_get_header_format_tag(tag_in, tag_header_in);
	  if(ierr!=0) return ierr;
	ierr = sep_get_header_format_tag(tag_out, tag_header_out);
	  if(ierr!=0) return ierr;
	ierr=sep_copy_data_pointer(*tag_header_in, *tag_header_out);
	  if(ierr!=0) return ierr;

	ierr=sep_get_number_header_axes(tag_in,&num_axis);
	  if(ierr!=0) return ierr;
    /* Copy number of headers */
	for(i1=2;i1<=num_axis;i1++){
/*	i1=2;*/
	ierr=sep_get_header_axis_par(tag_in,&i1,&naxis,&oaxis,&daxis,label);
	  if(ierr!=0) return ierr;
	ierr=sep_put_header_axis_par(tag_out,&i1,&naxis,&oaxis,&daxis,label);
	  if(ierr!=0) return ierr;
	}

    /* Copy keys */
   	ierr=sep_copy_header_keys(tag_in, tag_out);
	  if(ierr!=0) return ierr;

	 free(*tag_header_in);
	 free(*tag_header_out);
    return 0;
}

