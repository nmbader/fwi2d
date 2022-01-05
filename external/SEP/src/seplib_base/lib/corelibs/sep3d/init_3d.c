/*$

=head1 NAME

init_3d - initialize sep3d i/o

=head1 SYNOPSIS

init_3d()

=head1 DESCRIPTION

Initiates sep3d dataset

=head1 SEE ALSO

L<sep_3d_close>

=head1 LIBRARY

B<sep3d>

=cut




>*/ 
/*
-------------------------------------------------

Author: Robert Clapp, ESMB 463, 7230253

Date Created:Sat Dec 13 18:49:36 PST 1997

Purpose: 

*/	 
#include <string.h>
#include <sep3d.h> 
#include "streamlist.h"
#include "./../include/sep_main_internal.h"

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void init_3d(void)
_XFUNCPROTOEND
#else
void init_3d()
#endif 

{ 
int hff_in,gff_in;
char hff_input[1024],gff_input[1024];
streaminf *info,*info_hff,*info_gff;


/*FIRST DO A GETCH TO GET EVERYTHING STARTED AND SEE IF WE HAVE A USER
OVERRIDE FOR HFF VALUE*/
memset(hff_input,'\0',sizeof(hff_input));
hff_in=hetch("hff_in","s",hff_input);
hff_in=getch("hff_in","s",hff_input);
memset(gff_input,'\0',sizeof(gff_input));
gff_in=getch("gff_in","s",gff_input);


/*SHOULD CHECK HERE TO SEE IF INPUT EXISTS*/
info = tag_info("in", TAG_INQUIRE);  /* get info on this tag */
if(info == SEPPOINTNULL) return ; /* No Input so this routine is meaningless*/

/*maybe we shouldn't exit if hff_in or gff_in is specified, think about it*/


/*GET HFF TAG  FROM HISTORY FILE IF NOT OVERWRITTEN*/
if(hff_in<=0) hff_in=hetch("hff","s",hff_input);

/*IF HFF DEFINED START UP TAG */
if(hff_in==1 && 0!=strcmp("-1",hff_input)){

	/*Check to see if headerformat file has already been specified, if
   so give an error, this should be the FIRST thing a sep3d calls */

	if((info->headerformatfile) != SEPSTRNULL) 
	  seperr("headerformatfile already set in structure init_3d should be the first call in a sep3d program \n");
	
	/*now copy the value of hff into the info structure */
/*  info->headerformatfile = (char *) malloc(1+(int)strlen(hff_input));*/
  info->headerformatfile = expandnm(hff_input,(char *) NULL);
/*  strcpy(info->headerformatfile, hff_input);*/

  /*finally start-up the tag */
	auxin(hff_input);
	info_hff = tag_info(hff_input, TAG_INQUIRE);  /* get info on this tag */
}

/*GET GFF TAG */
if(gff_in<=0) gff_in=hetch("gff","s",gff_input);
 
/*IF GFF TAG DEFINED START UP TAG*/
if(gff_in==1 && 0!=strncmp("-1",gff_input,2)){
 /*Check to see if gridformatfile file has already been specified, if
   so give an error, this should be the FIRST thing a sep3d calls */

  if((info->gridformatfile) != SEPSTRNULL)
    seperr("gridformatfile already set in structure init_3d should be the first call in a sep3d program \n");

  /*now copy the value of gff into the info structure */
  info->gridformatfile = expandnm(gff_input,(char *) NULL);
/*  info->gridformatfile = (char *) malloc(1+(int)strlen(gff_input));*/
/*  strcpy(info->gridformatfile, gff_input);*/

  /*start-up the tag */
  auxin(gff_input);
	info_gff = tag_info(gff_input, TAG_INQUIRE);  /* get info on this tag */

	/*unpipe the grid */
	if(0!=make_unpipe(gff_input)) seperr("trouble unpiping the grid \n");
	
} 
return;
} 


/*  $Id: init_3d.c,v 1.2 2004/04/08 22:32:27 bob Exp $ */
