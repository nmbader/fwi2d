/*$

=head1 NAME

srite  - write an array to seplib tag


=head1 SYNOPSIS

C<int srite(tag, buffer, nbytes)>

=head1 INPUT PARAMETERS

=over 4

=item char*-  tag     

      tag of history file

=item int -   nbytes  

      number of bytes to be written

=item void*-  buffer  

      values to be written

=back

=head1 DESCRIPTION 

Writes out data to dataset attached to given tag

=head1 RETURN VALUES

=over 4

n= number of bytes written

=back

=head1 COMMENTS

The tag argument is either the string "out" or any tag appropriate for
use with auxout().  This means either an explicit filename or a
command line redirect parameter tag=filename.  Buffer is the location
of contiguous bytes from which the output will be written.  Unless an
end of file or I/O error is encountered, it is guaranteed that all
nbytes bytes will be written; this is true even for terminals and
pipes.  In any event the number of characters written is returned.

Sreed and srite now perform conversions from machine independent
external data formats to the internal data representation. For srite
this defaults to "xdr_float". This may be overridden by
"data_format=..." on the command line or (for tag="out" only) by a
call to C<set_format(tag, "...")>  The valid types are "xdr_float",
"xdr_integer", "xdr_byte", "native_float", "native_byte" and "vplot". 

If the internal representation is larger than the representation of
the external data then you must be careful to make nbytes the size
of the converted data.  e.g. on a cray system when
writing xdr_float data, nbytes refers to the number of bytes to be written.
This will be in 4-byte xdr_floats, but the internal storage will be in
8-byte cray floats, so the buffer must be 2*nbytes long.

=head1 SEE ALSO 

L<sreed>, L<sseek>, L<ssize>, L<auxclose> 

=head1 DIAGNOSTICS 

If an end of file is reached, the returned byte count may be less than
the requested amount.  The next call will return zero.  If the write
was otherwise unsuccessful the program will be terminated via
seperr().  Many conditions can generate an error: physical I/O errors,
bad buffer address, preposterous nbytes, tag not that of an output
file.

=head1 BUGS 

=head1 KEYWORDS 

write xdr output

=head1 LIBRARY

B<sep>

=cut



*/

/* Author Dave Nichols 9/94
	New version using the new seplib IO architecture 
Modified Martin Karrenbach 25-8-94
	Added srite2 for use by libsepattr.
Modified: Stew Levin  05/08/96
        Added explicit cast on last arg of xdr_vector calls
Modified: Robert CLapp 07/18/97
				Added casts and prototyping
Modified: Robert CLapp 06/01/97
				Changed to  GNU prototyping style
 */

#include <sepConfig.h>

#include <stdio.h>

#include "streamlist.h"
#include "sep_main_internal.h"
#include <sepcube.h>

#include <assert.h>


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
static int srite_xdr( streaminf*, void*,int,int);
_XFUNCPROTOEND 
#else
static int srite_xdr();
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int srite(const char *tag, void *buf, int nbytes)
_XFUNCPROTOEND 
#else
int srite(tag,buf,nbytes)
char *tag;
char *buf;
int nbytes;
#endif
{

streaminf *info;

assert( tag != 0 );
assert( buf != 0 );

if( nbytes== 0 ) return 0;

info = tag_info( tag, TAG_OUT );

assert( info != 0 );

/* open the dataset if it isn't already open */
if( info->ioinf == 0 ){
	(*info->open_func)(info, &(info->ioinf) );
	if( !info->valid ){seperr("srite(): invalid output tag %s\n",tag);}

}

if( ! info->ready_out ) sepstr_ready_out(info);

#ifdef TMC_DISK_ARRAY
if( is_tmc_array(buf) && info->iotype != CM_SDA_IO &&
                         info->iotype != CM_REG_IO ) {
     seperr("srite(): tag %s, CM array and not setup for CM I/O \n",tag);
}
if( IS_TMC_FMT(info->format_num) && info->iotype != CM_SDA_IO ){
     seperr("srite(): tag %s, Not setup for CM  I/O and raw CM format\n",tag);
}
#endif


#ifdef WORDS_BIGENDIAN

    return( (*info->write_func)(info,info->ioinf,buf,nbytes));

#else /* we may have to convert the data to xdr_formats */

if( IS_XDR_FMT( info->format_num ) ){
	return( srite_xdr(info,buf,nbytes,info->format_num));
}else{
     	return((int) (*info->write_func)(info,info->ioinf,buf,nbytes));
}

#endif

}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int srite_raw( char *tag, void *buf, int nbytes )
_XFUNCPROTOEND
#else
int srite_raw( tag, buf, nbytes )
char *tag;
void *buf;
int nbytes;
#endif
{
streaminf *info;

assert( tag != 0 );
assert( buf != 0 );

if( nbytes== 0 ) return 0;

info = tag_info( tag, TAG_OUT );

assert( info != 0 );

/* open the dataset if it isn't already open */
if( info->ioinf == 0 ){
	(*info->open_func)(info,&(info->ioinf) );
	if( !info->valid ){seperr("srite(): invalid output tag %s\n",tag);}
}

if( ! info->ready_out ) sepstr_ready_out(info);

return((int)(*info->write_func)(info,info->ioinf,buf,nbytes));
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
int srite2(char *tag, void *buf, int nbytes, char* format)
_XFUNCPROTOEND 
#else
int srite2(tag,buf,nbytes,format)
char *tag;
char *buf;
int nbytes;
char* format;
#endif
{

streaminf *info;
int format_num;

assert( tag != 0 );
assert( buf != 0 );

if( nbytes== 0 ) return 0;

info = tag_info( tag, TAG_OUT );

assert( info != 0 );


if( (format_num = get_format_num( format )) == -1 ){
    seperr( "srite2(): Unknown format name \"%s\" for tag \"%s\" \n",
	     format,tag);
}

info->format_num=format_num;

/* open the dataset if it isn't already open */
if( info->ioinf == 0 ){
	(*info->open_func)(info, &(info->ioinf) );
	if( !info->valid ){seperr("srite(): invalid output tag %s\n",tag);}

}

if( ! info->ready_out ) sepstr_ready_out(info);

#ifdef TMC_DISK_ARRAY
if( is_tmc_array(buf) && info->iotype != CM_SDA_IO &&
                         info->iotype != CM_REG_IO ) {
     seperr("srite(): tag %s, CM array and not setup for CM I/O \n",tag);
}
if( IS_TMC_FMT(info->format_num) && info->iotype != CM_SDA_IO ){
     seperr("srite(): tag %s, Not setup for CM  I/O and raw CM format\n",tag);
}
#endif


#ifdef WORDS_BIGENDIAN

    return((int) (*info->write_func)(info,info->ioinf,buf,nbytes));

#else /* we may have to convert the data to xdr_formats */

if( IS_XDR_FMT( format_num ) ){
	return( srite_xdr(info,buf,nbytes,format_num));
}else{
     	return((int) (*info->write_func)(info,info->ioinf,buf,nbytes));
}

#endif

}



#define XDRSIZE SEP_BUFSIZ
char outxdrbuf[XDRSIZE];

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN 
static int srite_xdr( streaminf *info, void *buf, int nbytes, int format_num)
_XFUNCPROTOEND 
#else
static int srite_xdr( info, buf, nbytes, format_num)
streaminf* info;
void* buf;
int nbytes;
int format_num;
#endif
{
        unsigned int nwrite, total,size;
	char *cpos;
  int ierr;

	/* the memory xdr stream used for decoding data */
	static XDR* outxdr=0;

	if( outxdr == 0 ) {
	    /* initialize memory xdr stream first time through */
	    outxdr= (XDR*) alloc( sizeof(XDR));
            xdrmem_create( outxdr, outxdrbuf, XDRSIZE, XDR_ENCODE );
        }

	nwrite = nbytes;
	total = 0;
	cpos= buf;

	do{
	   size = MIN(  nwrite-total, XDRSIZE-4 ); /*REMOVE THE -4 WHEN LINUX XDR_BYTES WORKS*/
	   xdr_setpos( outxdr, 0 );

	   switch( format_num ){
    	     /* case( FMT_XDR_COMPLEX ): */ 
    	     case( FMT_XDR_FLOAT ):
#if defined(CRAY)
                { int ierr, numel, bitoff, type;
                  numel = size/4; bitoff=0; type=2;
                  ierr = CRAY2IEG( &type, &numel, (float*)outxdrbuf, &bitoff, 
                                   (float*)cpos );   
                  if( ierr != 0 ){
                   seperr("xdr_rite(): cray convert error %d, tag \"%s\" \n",
                                        ierr,info->tagname);
                  }
                 
                }
#else
    		if( xdr_vector( outxdr,cpos,size/4,sizeof(float),(xdrproc_t) xdr_float) ==
							FALSE ){
		seperr("xdr_rite(): xdr error, tag \"%s\" \n",info->tagname);
		}
#endif
	        cpos += size/4*sizeof(float);
		break;
    	     case( FMT_XDR_INT ):
    		if( xdr_vector( outxdr,cpos,size/4,sizeof(int),(xdrproc_t) xdr_int) == 
							FALSE ){
		seperr("xdr_rite(): xdr error, tag \"%s\" \n",info->tagname);
		}	
	        cpos += size/4*sizeof(int);
		break;
    	     case( FMT_XDR_BYTE ):
    		if( xdr_bytes( outxdr, &cpos, &size, size ) == FALSE ){
		seperr("xdr_rite(): xdr error, tag \"%s\" \n",info->tagname);
		}	
	        cpos += size;
		break;
	     default:
		seperr("xdr_rite(): tag \"%s\" I don't grok this format\n",
			info->tagname);
		break;
	   }

     ierr=(int) (*info->write_func)(info,info->ioinf, (outxdrbuf), size );
	   if(size != ierr)
		seperr("xdr_rite(): I/O error on output for tag \"%s\" -- ierr=%d \n",
				info->tagname,ierr);
	   total +=size;
	
	}while(total<nwrite );

	return total;
}


_XFUNCPROTOBEGIN
size_t srite_big(const char *tag, void *buf, size_t nbytes)
_XFUNCPROTOEND

{
  size_t done,rett,bigt;
  int doo,big,ret;
  done=0;
  bigt=2*1000*1000*1000;
  big=bigt;
  rett=0;
  while(done < nbytes){
    if(nbytes-done < big) doo=nbytes-done;
    else doo=big;
    ret=srite(tag,(((char *) buf)+done),doo);
    rett+=ret;
    if(ret!=doo) return rett;
    done+=ret;
  }
  return rett;
}
