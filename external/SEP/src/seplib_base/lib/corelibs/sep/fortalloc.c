
/*
 * Interface for fortran dynamic allocation use by the sat preprocessor 
 * 
 * This documentation describes the use from a fortran program.
 * 
 * 
 * 1) 	integer function  fortalloc (  base, length )
 *      ---------------------------------------------
 * 
 * Input is a base address ( defined as character*1 base(1) ) and a length to 
 * allocate in bytes.
 * 
 * It returns an offset such that base(offset) is the address
 * of the start of newly allocated memory of the reqested length.
 * 
 * It returns 1 on error.
 * 
 * 
 * 2)	integer function fortfree( base(offset) )
 *      -----------------------------------------
 * 
 * Input is base(offset) where base was passed to the previous 
 * call and the offset is the value that it returned.
 * It frees the associated memory.
 * It returns 1 on error.
 *
 *============================================================================= 
 * 
 * 8-22-90 : Modified, Steve Cole
 *	     Force alignment on double word boundary. Caused alignment
 *	     errors on DEC3100.
 * 8-29-90 : Modified, Dave Nichols
 * 	     Fixed uninitialized pointer error introduced in previous fixes 
 * 7-20-97: Modified, Robert CLapp Moved stdlib ifdef include sepcube
 * 
 */

#include <sepConfig.h>
#include "sep_fortran_internal.h"
#include <stdlib.h>


#define WORDSIZE 8
/*<
 fortalloc

USAGE
int fortalloc (base,length)

INPUT PARAMETERS
base - int base addreess to begin allocation
length-int amount to allocate

RETURN VALUES
-1  = ERROR

DESCRIPTION
Fortran77 dynamic memory allocation. Mainly used by saw and sat

CATEGORY
Lib:Sep:Fortran

COMPILE LEVEL
DISTR
>*/
 
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int fortalloc( char *base, int len )
_XFUNCPROTOEND
#else
int fortalloc( base, len )
char *base;
int len;
#endif
{
	char *newmem;
	int offset;
	int allocsize;

  	if (len % WORDSIZE == 0) {
                allocsize = len;
        } else {
                allocsize = ( (len / WORDSIZE) + 1) * WORDSIZE;
        }

	if( ( newmem = malloc( allocsize ) ) == (char*)0 ) return 1;

#if (defined CRAY)
        offset = newmem - base + 8;
#else
        offset = newmem - base + 1;
#endif


	return offset;
}
/*<
 fortfree

USAGE
int fortfree (buf)

INPUT PARAMETERS
buf -  memory segment

DESCRIPTION
Frees Fortran77 dynamic 

CATEGORY
Lib:Sep:Fortran

COMPILE LEVEL
DISTR
>*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int fortfree( char *buf )
_XFUNCPROTOEND
#else
int fortfree( buf )
char *buf; 
#endif
{

/* note free doesn't return anything so if this isn't a valid address
   (i.e. not one that was previously alloc'd) free will silently do nothing 
   and your memory usage wil grow and grow and GROW !! */

	free( buf ); 

/* so in this implementation fortfree never returs an error code
   however sometime in the future it might */

	return 0;
}
