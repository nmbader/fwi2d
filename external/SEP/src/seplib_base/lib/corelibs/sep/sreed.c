/*$

=head1 NAME

sreed - read in an array from a seplib file

=head1 SYNOPSIS

C<int sreed(tag, buffer, nbytes)>

=head1 INPUT PARAMETERS

=over 4

=item char*-  tag     

      tag of history file

=item int -   nbytes  

      number of bytes to be read

=back

=head1 OUTPUT PARAMETERS

=over 4

=item void*-  buffer  

      values to be read

=back

=head1 RETURN VALUES

=over 4

n = number of bytes read

=back

=head1 DESCRIPTION
Reads data from dataset attached to given tag

=head1 COMMENTS
The tag argument is either the string "in" or any tag appropriate for
use with auxin().  This means either an explicit filename or a command
line redirect parameter tag=filename.  Buffer is the location of
contiguous bytes into which the input will be placed. Unless an end
of file or I/O error is encountered, it is guaranteed that all nbytes
bytes will be read; this is true even for terminals and pipes. In any
event the number of characters read is returned.

sreed and srite perform conversions from machine independent external
data formats to the internal data representation. This is controlled by
the "data_format" keyword in the header file. The valid types are
"xdr_float", "xdr_integer", "xdr_byte" and "native". If the keyword is not
found the native format is assumed. 

If the internal representation is larger than the representation of 
the external data then you must be careful to make the buffer the 
correct size for the converted data.  e.g. on a cray system when 
reading xdr_float data, nbytes refers to the number of bytes to be read. 
This will be in 4-byte xdr_floats, but the internal storage will be in 
8-byte cray floats, so the buffer must be 2*nbytes long.

sreed_raw() just reads the raw bytes with no conversion. The buffer
size and nbytes are the same in this case.

If the returned value is 0, then
end-of-file has been reached.

=head1 SEE ALSO

seplib, file, L<srite>, L<auxclose>, L<auxpar>

=head1 DIAGNOSTICS
If an end of file is reached, the returned byte count may be less than
the requested amount.  The next call will return zero.  If the read
was otherwise unsuccessful the program will be terminated via
seperr().  Many conditions can generate an error: physical I/O errors,
bad buffer address, preposterous nbytes, file tag not that of an input
file.

=head1 BUGS

=head1 KEYWORDS

read input xdr

=head1 LIBRARY

B<sep>

=cut



*/

/*
Author Dave Nichols 8-94
	New version totally rewritten for new seplib.
Modified Dave Nichols 25-8-94
	Added seed2 for use by libsepattr.
Modified: Stew Levin  05-08-96
        Added explicit cast on last arg to xdr_vector calls
Modified: Robert Clapp 08-18-97 Added prototypes
Modified: Robert Clapp 6-1-99 Switched to GNU Prototypes
*/
#include <sepConfig.h>
#include <stdio.h>
#include <assert.h>

#include "streamlist.h"
#include <sepcube.h>
#include "sep_main_internal.h"

#if defined(HAVE_SYS_TYPES_H) || defined(__APPLE__)
#include <sys/types.h>
#endif
#if defined(HAVE_SYS_SOCKET_H) || defined(__APPLE__)
#include <sys/socket.h>
#endif

#if defined(HAVE_RPC_RPC_H) || defined(__APPLE__)
#include <rpc/rpc.h>
#endif

#if defined(HAVE_RPC_TYPES_H) || defined(__APPLE__)
#include <rpc/types.h>
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sreed_xdr( streaminf*, void*, int, int);
_XFUNCPROTOEND
#else
static int sreed_xdr();
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sreed(const char *tag, void *buf, int nbytes)
_XFUNCPROTOEND
#else
int sreed(tag,buf,nbytes)
char *tag;
char *buf;
int nbytes;
#endif
{
streaminf *info;



/*fprintf(stderr,"SREED %d \n",nbytes);*/

if( nbytes== 0 ) return 0;

assert( tag != 0 );
assert( buf != 0 );

info = tag_info( tag, TAG_IN );

assert( info != 0 );

/* check if this is the first I/O for this tag */
if( info->ioinf == 0 ){
    /* if it is then we open the dataset */
    (*info->open_func)( info, &(info->ioinf) );
     if( !info->valid ){seperr("sreed(): invalid input tag %s\n",tag);}
}

#ifdef TMC_DISK_ARRAY

if( is_tmc_array(buf) && info->iotype != CM_SDA_IO && 
			 info->iotype != CM_REG_IO  ){
     seperr("sreed(): tag %s, CM array and not setup for CM  I/O \n",tag);
}

if( IS_TMC_FMT(info->format_num) && info->iotype != CM_SDA_IO ){
     seperr("sreed(): tag %s, raw CM format and not setup for SDA I/O \n",tag);
}

#endif

/* check to see if the data is compressed, if so return error */

/* If native format is the same as xdr (SUN,IBM,HP etc.) we don't need
 * any conversions */


#ifdef WORDS_BIGENDIAN
     	return( (int)(*info->read_func)(info,info->ioinf, buf,nbytes));

#else

if( IS_XDR_FMT( info->format_num ) ){
	return( sreed_xdr(info,buf,nbytes,info->format_num));
}else{
     	return((int) (*info->read_func)(info,info->ioinf, buf,nbytes));
}
#endif

}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sreed_raw( char *tag, void *buf, int nbytes )
_XFUNCPROTOEND
#else
int sreed_raw( tag, buf, nbytes )
char* tag;
void* buf;
int nbytes;
#endif
{
streaminf* info;

assert( tag != 0 );
assert( buf != 0 );

info = tag_info( tag, TAG_IN );
assert( info != 0 );

/* check if this is the first I/O for this tag */
if( info->ioinf == 0 ){
    /* if it is then we open the dataset */
    (*info->open_func)( info, &(info->ioinf) );
    if( !info->valid ){seperr("sreed(): invalid input tag %s\n",tag);}
}

return((int) (*info->read_func)( info,info->ioinf,  buf, nbytes ));

}

_XFUNCPROTOBEGIN
size_t sreed_big(const char *tag, void *buf, size_t nbytes)
_XFUNCPROTOEND

{
  size_t done,rett,bigt;
  int doo,big,ret;
  done=0;
  bigt=2*1000*1000;
  big=bigt;
  rett=0;
  while(done < nbytes){
  
    if(nbytes-done < big) doo=nbytes-done;
    else doo=big;
    //fprintf(stderr,"sreed %lld %lld %lld \n",done,nbytes,doo);
    ret=sreed(tag,(((char *)buf)+done),doo);
    rett+=ret;
    if(ret!=doo) return rett;
    done+=ret;
  }
  return rett;
}
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sreed2(char *tag, void *buf, int nbytes, char* format)
_XFUNCPROTOEND
#else
int sreed2(tag,buf,nbytes,format)
char *tag;
char *buf;
int nbytes;
char *format;
#endif
{
streaminf *info;
int format_num;

if( nbytes== 0 ) return 0;

assert( tag != 0 );
assert( buf != 0 );

info = tag_info( tag, TAG_IN );

assert( info != 0 );

if( (format_num = get_format_num( format )) == -1 ){
    seperr( "sreed2(): Unknown format name \"%s\" for tag \"%s\" \n",
	     format,tag);
}

/* check if this is the first I/O for this tag */
if( info->ioinf == 0 ){
    /* if it is then we open the dataset */
    (*info->open_func)( info, &(info->ioinf) );
     if( !info->valid ){seperr("sreed(): invalid input tag %s\n",tag);}
}

#ifdef TMC_DISK_ARRAY

if( is_tmc_array(buf) && info->iotype != CM_SDA_IO && 
			 info->iotype != CM_REG_IO  ){
     seperr("sreed(): tag %s, CM array and not setup for CM  I/O \n",tag);
}

if( IS_TMC_FMT(info->format_num) && info->iotype != CM_SDA_IO ){
     seperr("sreed(): tag %s, raw CM format and not setup for SDA I/O \n",tag);
}

#endif

/* If native format is the same as xdr (SUN,IBM,HP etc.) we don't need
 * any conversions */

#ifdef WORDS_BIGENDIAN
     	return((int) (*info->read_func)(info,info->ioinf, buf,nbytes));

#else


if( IS_XDR_FMT( format_num ) ){
	return( sreed_xdr(info,buf,nbytes,format_num));
}else{
     	return((int) (*info->read_func)(info,info->ioinf, buf,nbytes));
}
#endif

}

#define XDRSIZE SEP_BUFSIZ
char inxdrbuf[XDRSIZE]; 

/* read "nbytes" converted bytes into the buffer "buf" 
  The length of the buffer should be large enough to hold the result which
  might be longer than nbytes (e.g. xdr_float I/O on a cray )
 */

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
static int sreed_xdr( streaminf *info, void *buf, int nbytes, int format_num )
_XFUNCPROTOEND
#else
static int sreed_xdr( info, buf, nbytes, format_num )
streaminf* info;
void* buf;
int nbytes;
int format_num;
#endif
{
  unsigned int nread, total,size;
	long num;
	char* cpos;
	static XDR* inxdr=0;

	/* the memory xdr stream used for decoding data */

	if( inxdr == 0 ) {
	    /* initialize memory xdr stream first time through */
	    inxdr= (XDR*) malloc( sizeof(XDR));
            xdrmem_create( inxdr, inxdrbuf, XDRSIZE, XDR_DECODE );
        }

	nread = nbytes;
	total = 0;
	cpos = buf;


	do{
	   size = MIN(  nread-total, XDRSIZE );
	   xdr_setpos( inxdr, 0 );
	   num = (*info->read_func)(info,info->ioinf,  (inxdrbuf), size );

	   /* eof */
	   if( num == 0 ) return total;

	   switch( format_num ){
    	     case( FMT_XDR_FLOAT ):
#if defined(CRAY)   /* use optimised routines on the cray */
	        { int ierr, numel, bitoff, type;
		  numel = num/4; bitoff=0; type=2;
		  ierr = IEG2CRAY( &type, &numel, (float*)inxdrbuf, &bitoff, 
				   (float*)cpos );   
		  if( ierr != 0 ){
		   seperr("xdr_reed(): cray convert error %d, tag \"%s\" \n",
					ierr,info->tagname);
                  }
		  
		}
#else
    		if( xdr_vector( inxdr,cpos,(int)(num/4),sizeof(float),(xdrproc_t) xdr_float) ==
							FALSE ){
		seperr("xdr_reed(): xdr error, tag \"%s\" \n",info->tagname);
		}
#endif
	        cpos += num/4*sizeof(float);
		break;
    	     case( FMT_XDR_INT ):
    		if( xdr_vector( inxdr,(char *) (&(cpos[0])),(int)(num/4),sizeof(int),(xdrproc_t) xdr_int) == 
							FALSE ){
		seperr("xdr_reed(): xdr error, tag \"%s\" \n",info->tagname);
		}	
	        cpos += num/4*sizeof(int);
		break;
    	     case( FMT_XDR_BYTE ):
    		if( xdr_bytes( inxdr, &cpos, (u_int *) &num, (u_int) num ) == FALSE ){
		seperr("xdr_reed(): xdr error, tag \"%s\" \n",info->tagname);
		}	
	        cpos += num;
		break;
	     default:
		seperr("xdr_reed(): tag \"%s\" I don't grok this format\n",
			info->tagname);
		break;
	   }
	   total += num;
	
	}while(total<nread  );

	return total;
}
