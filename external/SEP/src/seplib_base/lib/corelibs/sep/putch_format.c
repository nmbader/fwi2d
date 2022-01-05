/*
 * formatting routine for putch(), auxputch(), snap()
 *
 * Author: Stewart A. Levin  12/28/85  Cannibalized from putch_.c
 * Revised: Stewart A. Levin  1/7/86   Greater precision for floats
 * Revised: Stewart A. Levin  5/23/87  Added type="1" logical conversion
 * Revised: Joe Dellinger   Feb 26 92  Made type="1" write out "yes" or "no"
 *                                     instead of "0" or "1"
 * Revised 12-18-97 biondo  made type='l' equivalent to type='1' 
   for saw compatibility
 */
#include <sepConfig.h>
#include    <stdio.h>
#include "fastpar.h"
#include <sepcube.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void putch_format( char* string, char* name, char* type, MIXED ptr)
_XFUNCPROTOEND
#else
void putch_format(string,name,type,ptr)
char *string, *name, *type;
MIXED ptr;
#endif
{
	char c;
	int  integer;
	float real;
	double dubble;
        long long lg;

	c = *type;

	switch(c) {
	case 'm':
		lg = *ptr.m;
		sprintf (string,"\t\t%s=%lld\n",name,lg);
		break;
	case 'd':
	case 'i':
		integer = *ptr.i;
		sprintf (string,"\t\t%s=%d\n",name,integer);
		break;
	case 'l':
	case '1':
		/*
		 * Explicitly write out "yes" or "no"
		 * so people know it is a boolean and not an int.
		 */
		if (*ptr.i)
		    sprintf (string,"\t\t%s=%s\n",name,"yes");
		else
		    sprintf (string,"\t\t%s=%s\n",name,"no");
		break;
	case 'f':
	case 'r':
		real = *ptr.f;
		sprintf (string,"\t\t%s=%g\n",name,real);
		break;
	case 'g':
		dubble = *ptr.g;
		sprintf (string,"\t\t%s=%.15g\n",name,dubble);
		break;
	case 's':
		if(*ptr.s != '"')  /* Add quotes around string for insurance */
		    sprintf(string,"\t\t%s=\"%s\"\n",name,ptr.s);
		else
		    sprintf (string,"\t\t%s=%s\n",name,ptr.s);
		break;
	default:
		seperr("putch_format(): unknown conversion type %d\n%s\n%s\n",c,name,type);
	}
}
