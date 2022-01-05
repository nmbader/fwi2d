#define SET_SDOC 1 
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Sun Aug 16 15:04:02 PDT 1998

Purpose: 

*/	 

#include "superset_internal.h"
#include <superset.h> 

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
<
sep3d_print_info

USAGE
ierr= sep3d_print_info(char *sep3dname)

INPUT PARAMETERS
sep3dname -  char*   tag associated with sep3d structure



RETURN VALUES
0   =  if it works correctly
(check superset.h for other return values and there meaning)

DESCRIPTION
Prints information about the dataset
>
*/






#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_print_info(char *sep3dname)
_XFUNCPROTOEND
#else
int sep3d_print_info(sep3dname)
char *sep3dname;
#endif 
{
int i1,ierr;
sep_3d *info;


info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(info == SEPNULL)
  return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));


sep3d_print(info);
return(SUCCESS);
}
