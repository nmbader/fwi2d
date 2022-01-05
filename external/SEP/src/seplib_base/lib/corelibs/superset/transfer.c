#define SET_SDOC 1
#include<sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include<superset.h>
#include<superset_internal.h>
#ifdef SEP_MPI
#include<mpi.h>
#endif
#ifndef MAX
#define MAX(a,b) ( ((a)>(b)) ? (a):(b) )
#endif
#ifndef MIN
#define MIN(a,b) ( ((a)<(b)) ? (a):(b) )
#endif

extern int sep_thread_num(void);


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_broadcast_headers(char *sep3dname, int ifrom)
_XFUNCPROTOEND
#endif
{
sep_3d *info;
#ifdef SEP_MPI
MPI_Status status;
MPI_Comm comm;
#endif
char usage[256],data[256],file[256];
int nsz,i,k,iam,impi;
int n,nbuf;
float o,d;
char temp1[256],temp2[256],temp3[256];
char *buf;



#ifdef SEP_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &impi);
if(ifrom==impi){
  info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
  if(info == SEPNULL)
    return (sepwarn(INVALID_STRUC,"tag:%s 2invalid struc\n",sep3dname));
  
  nbuf=SEP3D_nconvert(info);

}
else{
  info = tag_info_sep3d(sep3dname, SCRATCH);  /* get info on this tag */
  if(info!=SEPNULL) SEP3D_clean(info);
  info = tag_info_sep3d(sep3dname, SCRATCH);  /* get info on this tag */
  if(info == SEPNULL)
        return (sepwarn(INVALID_STRUC,"tag:%s 1invalid struc\n",sep3dname));
}
  MPI_Bcast( &nbuf, 1, MPI_INT,  ifrom,MPI_COMM_WORLD); 
  sep_mpi_stop();
  buf=(char*) malloc(sizeof(char)*(nbuf));
if(ifrom==impi){
    nbuf=0;
   if(0!=SEP3D_convert(info,&nbuf,buf))
     return(sepwarn(NOT_MET,"trouble converting tag  %s \n",info->name));
}
  MPI_Bcast((void*)buf,nbuf,MPI_BYTE, ifrom, MPI_COMM_WORLD);
sep_mpi_stop();
if(ifrom!=impi){
    nbuf=0;
   if(0!=SEP3D_unconvert(info,&nbuf,buf))
     return(sepwarn(NOT_MET,"trouble converting tag  %s \n",info->name));
free(buf);
}
if(ifrom==impi) free(buf);
#else
/*return(sepwarn(FAIL_OTHER,"Distribute header requested with a non-mpi version of SEPLIB \n"));*/
#endif


return(SUCCESS);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_pass_headers(int ithread,char *sep3dname, int ifrom, int ito)
_XFUNCPROTOEND
#endif
{
sep_3d *info;
#ifdef SEP_MPI
MPI_Status status;
#endif
char usage[256],data[256],file[256];
int nsz,i,k;
int n,impi,nbuf;
float o,d;
char temp1[256],temp2[256],temp3[256];
char *buf;


if(ifrom==ito)  return(SUCCESS);
if(ithread!=ifrom && ithread!=ito) return(SUCCESS);
#ifdef SEP_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &impi);
  info = tag_info_sep3d(sep3dname, INQUIRE);  /* get info on this tag */
if(ifrom==impi){
  if(info == SEPNULL)
    return (sepwarn(INVALID_STRUC,"tag:%s  invalid struc\n",sep3dname));
  
  nbuf=SEP3D_nconvert(info);
    MPI_Send( &nbuf, 1, MPI_INT,ito, 2321, MPI_COMM_WORLD);

}
else{
  if(info!=SEPNULL) SEP3D_clean(info);
  info = tag_info_sep3d(sep3dname, SCRATCH);  /* get info on this tag */
   MPI_Recv( &nbuf, 1, MPI_INT, ifrom, 2321,MPI_COMM_WORLD,&status);
}
  buf=(char*) malloc(sizeof(char)*nbuf);
if(ifrom==impi){
    nbuf=0;
   if(0!=SEP3D_convert(info,&nbuf,buf))
     return(sepwarn(NOT_MET,"trouble converting tag  %s \n",info->name));
    MPI_Send( buf, nbuf, MPI_BYTE, ito,2322,MPI_COMM_WORLD);
}
else{
   MPI_Recv( buf, nbuf, MPI_BYTE, ifrom, 2322,MPI_COMM_WORLD,&status);
   nbuf=0;
   if(0!=SEP3D_unconvert(info,&nbuf,buf))
     return(sepwarn(NOT_MET,"trouble unconverting tag %s \n",sep3dname));
}
free(buf);
#else
if(ifrom!=ito)
return(sepwarn(FAIL_OTHER,"Pass header requested with a non-mpi version of SEPLIB \n"));
#endif

return(SUCCESS);
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_broadcast_ints(int *val, int nsz, int ifrom)
_XFUNCPROTOEND
{
int impi;
#endif
#ifdef SEP_MPI
MPI_Status status;
#endif

#ifdef SEP_MPI
 MPI_Bcast(val,nsz,MPI_INT, ifrom, MPI_COMM_WORLD);
#endif
return(SUCCESS);
}




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_send_int(int ithread,int *val, int ifrom, int ito)
_XFUNCPROTOEND
{
int impi;
#endif
#ifdef SEP_MPI
MPI_Status status;
#endif

if(ifrom==ito)  return(SUCCESS);
if(ithread!=ifrom && ithread!=ito) return(SUCCESS);
#ifdef SEP_MPI
  if(ifrom==ithread) MPI_Send( val, 1, MPI_INT, ito, 2329,MPI_COMM_WORLD);
  else MPI_Recv( val, 1, MPI_INT, ifrom, 2329,MPI_COMM_WORLD,&status);
#else
if(ifrom!=ito)
return(sepwarn(FAIL_OTHER,"Pass integer requested with a non-mpi version of SEPLIB \n"));
#endif
return(SUCCESS);
}






#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_broadcast_data(int ifrom,  int esize, int nelem, char *data)
_XFUNCPROTOEND
{
sep_3d *info;
#endif
#ifdef SEP_MPI
MPI_Status status;
#endif
int chunk,done;


done=0;
while(done < nelem*esize){
  chunk=MIN(nelem*esize-done,1000000);
#ifdef SEP_MPI
   MPI_Bcast( (data+done), chunk, MPI_BYTE, ifrom, MPI_COMM_WORLD);
#else
/*return(sepwarn(FAIL_OTHER,"Broadcast data requested with a non-mpi version of SEPLIB \n"));*/
#endif
done=done+chunk;
sep_mpi_stop();
}
return(SUCCESS);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int sep3d_pass_data(int ithread, int ifrom, int ito, int esize, int nelem, char *data)
_XFUNCPROTOEND
{
sep_3d *info;
#endif
#ifdef SEP_MPI
MPI_Status status;
#endif
int done,chunk;

if(ifrom==ito)  return(SUCCESS);
if(ithread!=ifrom && ithread!=ito) return(SUCCESS);
done=0;
while(done < nelem*esize){
  chunk=MIN(nelem*esize-done,1000000);
#ifdef SEP_MPI
  if(ifrom==ithread) MPI_Send( (data+done), chunk, MPI_BYTE, ito, 2329,MPI_COMM_WORLD);
  else MPI_Recv( (data+done), chunk, MPI_BYTE, ifrom, 2329,MPI_COMM_WORLD,&status);
#else
if(ifrom!=ito)
return(sepwarn(FAIL_OTHER,"Pass data requested with a non-mpi version of SEPLIB \n"));
#endif
  done+=chunk;
}

return(SUCCESS);
}

int sep3d_pass_string(char *string, int ifrom,int ito){
int nlen;
#ifdef SEP_MPI
MPI_Status status;

 if(sep_thread_num()==ifrom){
    nlen=(int)(strlen(string)+1);
    MPI_Send( &nlen, 1, MPI_INT,ito, 2321, MPI_COMM_WORLD);
    MPI_Send( string, nlen, MPI_CHAR,ito, 2322, MPI_COMM_WORLD);
 } 
 else if(sep_thread_num()==ito){
   MPI_Recv( &nlen, 1, MPI_INT, ifrom, 2321,MPI_COMM_WORLD,&status);
   MPI_Recv( string, nlen, MPI_CHAR, ifrom, 2322,MPI_COMM_WORLD,&status);
 }
#endif
return(SUCCESS);
}

int sep3d_broadcast_string(char *string, int ifrom){
int nlen;

if(sep_thread_num()==ifrom) nlen=(int)(strlen(string)+1);
#ifdef SEP_MPI
 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Bcast(&nlen,1,MPI_INT, ifrom, MPI_COMM_WORLD);
 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Bcast(string,nlen,MPI_INT, ifrom, MPI_COMM_WORLD);
 MPI_Barrier(MPI_COMM_WORLD);
#endif


return(SUCCESS);
}

void sep_mpi_stop(void){

#ifdef SEP_MPI
MPI_Barrier(MPI_COMM_WORLD);
#endif


}
int sep_num_thread(void){
int i;
#ifdef SEP_MPI
MPI_Comm_size(MPI_COMM_WORLD, &i);
return(i);
#else
return(1);
#endif
}


int sep_thread_num(void){
int i;
#ifdef SEP_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &i);
   return(i);
#else
return(0);
#endif
}

