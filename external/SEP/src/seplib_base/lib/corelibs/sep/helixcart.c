
/*$

=head1 NAME

h2c - convert from helix to cartersian coordinates

=head1 SYNOPSIS

C<h2c(int hindex, int *n, int ndim, int *cindex)>


=head1 INPUT PARAMETERS

=over 4

=item hindex -   int      

      helix location index

=item n -   int*     

      dimension of dataset

=item ndim -   int*     

      skip between elements to copy

=back

=head1 OUTPUT PARAMETERS

=over 4

=item cindex -   int*     

       location in multi-d mesh

=back



=head1 DESCRIPTION

Translate from helix coordinate to multi-dimension system

=head1 SEE ALSO

L<c2h>

=head1 LIBRARY

B<sep>

=cut

>
*/

#include<sep_main_external.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void h2c(int hindex, int *n, int ndim,int *cindex)
_XFUNCPROTOEND
#else
void h2c(hindex, n,ndim,cindex)
int hindex,  *n,  ndim, *cindex;
#endif
{
int block,i1;
float tempr;

block=1;
for(i1=0; i1 < ndim; i1++){
  cindex[i1]=(hindex/block)%n[i1];
  block=block*n[i1];
}

return;
}
/*$

=head1 NAME

c2h - convert from cartesian 

=head1 SYNOPSIS

C<c2h(int *hindex, int *n, int ndim, int *cindex)>


=head1 INPUT PARAMETERS

=over 4

=item n -   int*     

      dimension of dataset

=item ndim -   int*     

      skip between elements to copy

=item cindex-   int*     

      location in multi-d mesh

=back

=head1 OUTPUT PARAMETERS

=over 4

=item hindex -   int*     

      helix location index

=back

=head1 DESCRIPTION

Translate from multi-dimension system to  helix coordinate

=head1 SEE ALSO

L<h2c>

=head1 LIBRARY

B<sep>

=cut

>*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void c2h(int *hindex, int *n, int ndim,int *cindex)
_XFUNCPROTOEND
#else
void c2h(hindex, n,ndim,cindex)
int *hindex,  *n,  ndim, *cindex;
#endif
{
int block,i1;
float tempr;

block=1;
*hindex=0;
for(i1=0; i1 < ndim; i1++){
  *hindex+=cindex[i1]*block;
  block=block*n[i1];
}
return;
}

