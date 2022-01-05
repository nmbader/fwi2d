#define SET_SDOC 1
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Sun Aug 16 13:52:22 PDT 1998

Purpose: 

*/	 
#include <sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include<superset.h> 
#include"superset_internal.h" 

/*<
tag_info_sep3d

Usage
sep3d=tag_info_sep3d(sep3dname,usage)

Input Paramters
sep3dname   -   char*         pointer to sep_3d structure  
usage       -   usage_type    usage (INPUT,OUTPUT,SCRATCH) for structure
Output Paramters

Return Values

Description
Finds and initializes a tag

>*/




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
sep_3d* tag_info_sep3d(const char *name, enum usage_type type)
_XFUNCPROTOEND
#else
sep_3d* tag_info_sep3d(name,type) 
char *name ;
enum usage_type type;
#endif 
{ 
 sep_3d *curr;

    curr =sep3d_head();


		
    while( curr != SEPNULL ){
    if( strcmp(name,curr->name) == 0 ){
      return curr;
  }
  curr = curr->next;
    }


   if(type==INQUIRE) return(SEPNULL);

/* fell through so we must create it and put it at the end of the list*/
    curr = sep3d_new( name, type );
		sep3d_addend(curr);

    switch( curr->usage ){
     case INPUT:
        break;
     case OUTPUT:
        break;
     case SCRATCH:
/*				if(NULL==auxinout(name))*/
/*					seperr("trouble opening %s as tagsrc \n",name);*/
        break;
		 case INQUIRE:
				return(SEPNULL);
				break;
    }

    return curr;
}

/*  $Id: tagstream.c,v 1.2 2004/04/08 22:32:28 bob Exp $ */
