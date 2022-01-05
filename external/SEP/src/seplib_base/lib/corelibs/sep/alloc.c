/*$

=head1 NAME

alloc  - allocate C array with error checking

=head1 SYNOPSIS

C<char *alloc(nbytes)>

=head1 INPUT PARAMETERS

=over 4

=item nbytes - integer 

      Number of bytes to allocate

=back

=head1 DESCRIPTION

	Alloc provides dynamic core allocation (via malloc (3))
	with error checking.
	Alloc returns a pointer to a block of at least nbytes
	characters suitably aligned (after possible pointer coercion)
	for storage of any type of object.

=head1 COMMENTS

 	In order to allocate an array of floating point numbers, use
 	the following command in the calling routine:

		float *x;
		x = (float *) alloc(nx*sizeof(float));

 	nx is the number of elements needed in the array. 

=head1 SEE ALSO

	malloc (3)

=head1 DIAGNOSTICS

	Alloc terminates program execution with an appropriate 
	error message if core could not be allocated.

=head1 KEYWORDS    

alloc malloc memory allocation

=head1  LIBRARY

B<sep>

=cut

*/

/*
 *	allocation with error detection
 *
 * modified 1/26/83  S. Levin : added unsigned coersion to agree with malloc(3)
 *                              changed to NULL instead of 0
 * modified 3/19/83  S. Levin : keep internal 1024 byte block for small request
 *			include source for malloc inline to block missing
 *				cfree() linkage references.
 * modified 2/7/84   S. Levin : deleted calloc() call to reduce paging -- no
 *				longer is memory necessarily zeroed out of alloc
 * modified      94  D. Nichols: Add posix definition of malloc
 * modified   7/20/97 R. Clapp: Moved malloc deifintions ifdef to sepcube.h
 *                              malloc is used everywhere, seems to be a waste;
 */
#include <sep_main_external.h>
#include <stdlib.h>
#include <stdio.h>


#define SMALLBLOCK 0
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
char *alloc (size_t size)
_XFUNCPROTOEND
#else
char *alloc (size)
size_t size;
#endif
{
	static char *myblock = NULL; static int bytesleft = 0;
	char *ptr; 

	if(size >= SMALLBLOCK)
	   {
	    if ((ptr = malloc (size)) == (char*)NULL)
		seperr ("alloc(sep):can't allocate %d bytes\n",size);
	    return (ptr);
	   }
	else
	  { 
	    if(bytesleft < size)
	       {
	        if ((myblock = malloc ((unsigned)SMALLBLOCK)) == (char*)NULL)
		   seperr ("alloc(sep):can't allocate %d bytes\n",size);
		bytesleft = SMALLBLOCK - size;
	       }
	   else
	       {
		bytesleft -= size;
	       }
	   ptr = myblock;
	   myblock += size;
	   return(ptr);
	  }
}
