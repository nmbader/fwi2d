#ifndef STRFORMAT_H
#define STRFORMAT_H
#include <prototypes.h>

/* 
 *  Defines the formats possible for seplib files
 *  Used in the switch in sreed/srite. 
 */

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern int get_format_num( char* format_nam );
extern char* get_format_nam( int format_num );
extern int get_format_byte_length(const char* format_nam);
_XFUNCPROTOEND
#else
extern int get_format_num();
extern char* get_format_nam();
extern int get_format_byte_length();
#endif


#define FMT_NOT_KNOWN 	       -1
#define FMT_XDR_COMPLEX 		0
#define FMT_XDR_FLOAT 		0
#define FMT_XDR_INT 		1
#define FMT_XDR_BYTE 		2
#define FMT_NATIVE_BYTE 	3
#define FMT_NATIVE_FLOAT 	4
#define FMT_VPLOT 		5
#define FMT_TMC_FLOAT		6
#define FMT_TMC_BYTE		7

#define         NUM_FMT      8   /* FMT_NOT_KNOWN doesn't count */

#define IS_XDR_FMT( x )  ( x> -1 && x< 3 ) 

#define IS_TMC_FMT( x )  ( x> 5 && x< 8 ) 

#endif

