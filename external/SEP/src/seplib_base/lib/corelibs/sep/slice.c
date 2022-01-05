/*<
slice

USAGE
int snap(name,esize,n1,n2,matrix)

INPUT PARAMETERS
char* - name tag of desired output file
int   - esize size in bytes of one element of output
int   - n1   size of slice in fast dimension 
int   - n2   size of slice along slow axess
void* - matrix value of slice

RETURN VALUES
  0 = if successful


DESCRIPTION
Slice, and snap provide a convenient method of output for
snapshots of intermediate computations or iterations in a form
suitable for processing by other seplib functions (such as Window,
Movie, etc.)  None of the inputs is modified.

COMMENTS
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

=head1 DIAGNOSTICS

Program execution is terminated with an appropriate error message
if the program attempts to change the snapshot size along the way.
Also for appropriate I/O related errors.

=head1 BUGS 

Prexisting snapshot files with the same name are overwritten
without qualms.

=head1 SEE ALSO

L<sreed>, L<datapath>

=head1 KEYWORDS 

snapshot slice output
*/
/* slice()  --  Append a plane to a data file and update the History file.
	call slice("name",esize,n1,n2,data)
	(nt,nx,np) are synonyms for (n1,n2,n3) both on input and output.
	History file updated correctly if n2 or n1 or both =1

slice:
	getpar("name") to test for command line over-ride "name=newname"
	If name does not exist,
		create it.
	If name has zero length,
		Deposit  "title=name  n1=nt=  n2=nx=  in=datapathDname"  in name
		Create datapathDname
	Deduce n3 from length of file datapathDname.
	Append data to datapathDname.
	Update n3= in name.

  ---makefile example----

out name [name2 ...] [newname]: yourprocess
	Zero name [name2 ...]	 [newname]		#  cat /dev/null > name
	yourprocess [name=newname] > out

Author: jon
Revised: stew  3/19/83  changed method of determining whether to append or create 
		by keeping track of what names have already been sliceped.
Revised: stew  3/21/83  improved core allocation scheme to be more efficient.
Revised: stew  3/26/83  call noheader()
Revised: stew  7/14/83  added seperr () and perror() diagnostics
Revised: stew  11/29/83 update rather than append to slice header
Revised: stew  4/29/84  Changed final ^D to \n. This parallels hclose() mod.
			Allow name=Hnewname as well as name=newname on cmdline.
Revised: stew  8/29/84  insure n3= stored in header
Revised: stew  2/10/85  assume name=Hnewname(or whatever) without modification
Revised: stew  2/21/85  fixed up previous fix correctly
Revised: stew  8/22/85  Time stamp output header

Author: stew   12/12/85  Copied snap to slice for always create new header/data.
Revised: jon	1/19/86	Copied slice to slice, esize now an argument.
Revised: stew	2/19/86	Noop if n1 or n2 <= 0.
Modified 9/16/90  Dave Nichols make ansi-c and posix compatible.
Modified 8-94  Dave Nichols rewrote using new seplib internals

*/
#include <sepConfig.h>
#include <stdio.h>
#include "streamlist.h"
#include "sep_main_internal.h"
#include <sep_main_external.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int slice (char* name,int esize,int n1,int n2,void* data)
_XFUNCPROTOEND 
#else
int slice (name,esize,n1,n2,data)
char *name;
int  esize, n1, n2;
void *data;
#endif
{
    streaminf* info;
    int oldn1, oldn2, oldn3, n3, oldesize;
    int nwrite;

    /* open in inout mode */
    info = tag_info( name, TAG_INOUT );

    if( n1 == 1 ){
	n1=n2; n2=1;
    }

    if( !info->ready_out ){ /* first time through for this tag */

      if( info->headerbuf != 0 ){ /* header file already exists */
	if( !sepstrpar( info, "n1", "d", &oldn1 ) ) 
	  seperr("slice(): n1 required in existing header for tag \"%s\"\n",
		 info->tagname);
	if( !sepstrpar( info, "n2", "d", &oldn2 ) )
	  seperr("slice(): n2 required in existing header for tag \"%s\"\n",
		 info->tagname);
	if( !sepstrpar( info, "n3", "d", &oldn3 ) )
	  seperr("slice(): n3 required in existing header for tag \"%s\"\n",
		 info->tagname);
	if( !sepstrpar( info, "esize", "d", &oldesize ) )
	  seperr("slice(): esize required in existing header for tag \"%s\"\n");


       }else{
 	/* this is a new header file */
	oldn1 = n1; sepstrput( info, "n1", "d", &oldn1 );
	oldesize=esize; sepstrput( info, "esize", "d", &esize ); 
	if( n2 ==1 ){
	   /* we will add to the two axis */
	   oldn2 = 0;
	   oldn3=1; sepstrput( info, "n3", "d", &oldn3 ); 
        }else{
	   oldn2 = n2; sepstrput( info, "n2", "d", &oldn2 ); 
	   oldn3 = 0;
	}
       }

       /* get ready for output */
       sepstr_ready_out( info );

    }else{
	/* second and susequent times through, obtain, n1, n2, n3, esize */
	sepstrpar( info, "n1", "d", &oldn1 );  
	sepstrpar( info, "n2", "d", &oldn2 );
	sepstrpar( info, "n3", "d", &oldn3 );
	sepstrpar( info, "esize", "d", &oldesize );
    }

     /* check that old and new parameters match */
    if( oldn1 != n1 ) 
       seperr("slice(): n1=%d doesn't match earlier n1=%d for tag \"%s\"\n",
			     n1, oldn1, info->tagname );

    if( oldesize != esize ) 
       seperr("slice(): esize=%d doesn't match earlier esize=%d, tag \"%s\"\n",
			     esize, oldesize, info->tagname );

    if( n2 == 1  ){
       /* figure out new value for n2 */
       n2 = oldn2 +1;

       /* rewrite the last line of the header */
       sepstrputlast( info, "n2", "d", &n2 );

	nwrite = n1*esize;

    }else{

       if( oldn2 != n2 ) 
       seperr("slice(): n2=%d doesn't match earlier n2=%d for tag \"%s\"\n",
			     n2, oldn2, info->tagname );

       /* figure out new value for n3 */
       n3 = oldn3 +1;

       /* rewrite the last line of the header */
       sepstrputlast( info, "n3", "d", &n3 );

	nwrite = n1*n2*esize;
    }

    /* seek to the end of the datafile and write the new data */
    sseek( name, 0L, 2 );
    srite( name, data, nwrite );
		return(0);

}
