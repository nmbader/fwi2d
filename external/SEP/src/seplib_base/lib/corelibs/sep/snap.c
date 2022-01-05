/*$
=head1 NAME

snap,slice - write to the screen an array (through Grey etc)

=head1 SYNOPSIS

C<int snap(name,n1,n2,matrix)>

=head1 INPUT PARAMETERS

=over 4

=item char* - name 

      tag of desired output file

=item int   - n1   

      size of slice in fast dimension

=item int   - n2   

      size of slice along slow axess

=item void* -  matrix

      matrix value of slice

=back

=head1 RETURN VALUES

=over 4

	0 = if successful

=back
 
=head1 DESCRIPTION

Slice, and snap provide a convenient method of output for
snapshots of intermediate computations or iterations in a form
suitable for processing by other seplib functions (such as Window,
Movie, etc.)  None of the inputs is modified.

=head1 COMMENTS
The input tag name is used to locate and update an associated
auxiliary history file.  The name of the output header is "name".
This may be modified by specifying `name=newtag' on the command line
when a program is invoked.  Input is an n1 by n2 matrix of floating
point values containing the snapshot.  

If more than one call to these subroutines is made, the snapshot is
appended to the file and n3 is updated in the history file.  (This
frees the user from the chore of counting the number of output frames,
a special convenience while debugging.)

The subroutine datapath() is used to generate a name for an output
data file.

slice() only: esize is the element size, typically 1 or 4.

snap() only: If the output header does not exist no snapshots are written.
This allows the user to decide whether or not to create snapshots
on a given run.

=head1 FILES

name	auxiliary header

=head1 SEE ALSO

seplib, datapath, srite

=head1 DIAGNOSTICS

Program execution is terminated with an appropriate error message
if the program attempts to change the snapshot size along the way.
Also for appropriate I/O related errors.

=head1 BUGS 
Prexisting snapshot files with the same name are overwritten
without qualms.

=head1 KEYWORDS 

snapshot slice output

=head1 LIBRARY

B<sep>

=cut


*/
/*

Author: jon
Revised: stew  3/19/83  changed method of determining whether to append or create 
		by keeping track of what names have already been snapped.
Revised: stew  3/21/83  improved core allocation scheme to be more efficient.
Revised: stew  3/26/83  call noheader()
Revised: stew  7/14/83  added seperr () and perror() diagnostics
Revised: stew  11/29/83 update rather than append to snap header
Revised: stew  4/29/84  Changed final ^D to \n. This parallels hclose() mod.
			Allow name=Hnewname as well as name=newname on cmdline.
Revised: stew  8/29/84  insure n3= stored in header
Revised: stew  2/10/85  assume name=Hnewname(or whatever) without modification
Revised: stew  2/21/85  fixed up previous fix correctly
Revised: stew  8/22/85  Time stamp output header
Revised: stew  2/9/86 	Noop if n1 or n2 <= 0.
Revised: joe   8/13/86  Removed "Hname" stuff from documentation
Modified 8/22/89  Dave Nichols fix suballoc to return double word aligned data.
Modified 9-16-90  Dave Nichols  made ansi-c and posix compatible
Revised: W. Bauske 03-26-91 added comments to the program for the next person
			    has to port this code and figure it out
Modified 8-94  Dave Nichols rewrote using new seplib internals

*/
#include <sepConfig.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

#include "sep_main_internal.h"
#include <sepcube.h>


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int snap (char *name,int n1,int n2,void* data)
_XFUNCPROTOEND 
#else
int snap (name,n1,n2,data)
char *name;
int  n1, n2;
void *data;
#endif
{
char snapname[1024];
struct stat statbuf;

if( !getch(name,"s",snapname) ) strcpy( snapname, name );

 /* stat the header file name */
 if(-1 == stat(snapname,&statbuf)) {
      return 0;
 }else{
      slice(name, 4, n1, n2, data );
 }
return 0;
}
