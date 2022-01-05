/*$

=head1 NAME

auxpar - get parameter from auxilary file

=head1 SYNOPSIS

int auxpar (name, format, variable, tag)

=head1 INPUT PARAMETERS

=over 4

=item char* - name 

      parameter name

=item char* - type 

      format of variable (d,f,c,s,l)

=item char* - tag  

      name tag of history file

=back

=head1 OUTPUT PARAMETERS

=over 4

=item		 type* - variable 

         value of parameter

=back

=head1 RETURN VALUE

=over 4

=item		x = 

       number of values found matching given parameter last occurence

=back

=head1  DESCRIPTION

Read parameters to an auxilary history file

=head1 COMMENTS

       This function implements "hetch" for extracting values from
       expressions of the form `name=value' from an auxiliary input
       header. The fourth argument, tag, is used to determine the
       location of the auxiliary input header according the 
       following rules.

       Look for 'tag=header' on the command line.

       Use 'tag' as the auxiliary header name.

       When a keyword `par=filename' is encountered in the file, 
       that file (and any that it might in turn point to) is scanned 
       as well. Note auxpar() with tag="in" is now identical to hetch().

=head1 LIBRARY

B<sep>

=cut 
>*/

/*

SEE ALSO
      L<hetch>, L<auxin>

KEYWORDS   auxillary header parameter 

WHERE
       ./cube/libcube/auxpar.c

*/
/*
 * Author: Stewart A. Levin   4/29/84	From getch2()
 * Revised: stew 6/14/84  Needed to update _auxh definition
 * Revised: stew 9/6/87   Use <varargs.h> for portability
 * Revised: dave 9/17/90  Use stdarg for ANSI-C compilers
 */
#include <sepConfig.h>
#include <stdio.h>
#include "streamlist.h"
#include  <sep_pars_external.h>

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int auxpar( const char *name, const char *type, void* val, const char* tag )
_XFUNCPROTOEND 
#else
int auxpar( name, type,  val, tag )
char *name, *type, *tag;
char* val;
#endif
{
streaminf *info;

info = tag_info( tag, TAG_IN );

return sepstrpar(info,(char *) name,(char *) type,val); 
}
