#include <sepConfig.h>
#include <stdio.h>
#if defined(HAVE_SYS_SOCKET_H)
#include <sys/socket.h>
#endif
                                                                                 
#include <seplib.h>
#include <sep_main_external.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <netinet/in.h>

#define BUFFSIZE 1024
int sep_send_msg(char *status,char *extra) {
  int sock,port;
  struct sockaddr_in server;
  char buffer[BUFFSIZE];
  unsigned int echolen;
  int received = 0,seperrit=1;
  char jobid[1024];
  char master_ip[1024];
  char mach_label[1024];
  int gip,gport,gmach,gjob,ierr;

  if(0==strcmp(status,"error"))  seperrit=0;
  if(0==getch("sep.par_job","d",&gip))  gip=1;
  if(gip==0) return(0);

  gjob=getch("sep.jobid","s",jobid);
  gmach=getch("sep.mach_label","s",mach_label);
  gip=getch("sep.master_ip","s",master_ip);
  gport=getch("sep.master_port","d",&port);

  if(gjob ==1 || gmach ==1 || gip ==1 || gport ==1){
    if(gjob !=1 || gmach !=1 || gip !=1 || gport !=1)
      if(seperrit==1)
      seperr("when running in socket mode must provide jobid,mach_label,master_ip, and master_port %d %d %d %d",gjob,gmach,gip,gport);
      else return(1);




  /* Create the TCP socket */
  if ((sock = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
    if(seperrit==1) seperr("Failed to create socket \n");
    return(1);
  }
 
   sprintf(buffer,"%s:%s:%s:%s", status,jobid,mach_label,extra);

  /* Construct the server sockaddr_in structure */
  memset(&server, 0, sizeof(server));       /* Clear struct */
       server.sin_family = AF_INET;              /* Internet/IP */
  server.sin_addr.s_addr = inet_addr(master_ip);  /* IP address */
  server.sin_port = htons(port);       /* server port */
  /* Establish connection */
  if (connect(sock,
              (struct sockaddr *) &server,
              sizeof(server)) < 0) {
    if(seperrit==1) seperr("Failed to connect with server");
    else return(1);
  }
      

  /* Send the word to the server */
/*       echolen = strlen(buffer);*/
  echolen=strlen(buffer);
  if(echolen > BUFFSIZE-5) seperr("to long of message ");
  strcat(buffer,"EnDiT");
  ierr=send(sock, buffer, BUFFSIZE, 0);
  if (ierr != BUFFSIZE) {
    fprintf(stderr,"Mismatch in number of sent bytes %d !=%d \n",ierr,BUFFSIZE);
    if(seperrit==1) seperr("trouble sending message %s \n",buffer);
    else return(1);
  }
  ierr=recv(sock, buffer, 5, 0) ;
  if (ierr != 5) {
    fprintf(stderr,"Mismatch in number of received bytess %d!=5 ",ierr);
    if(seperrit==1)seperr("trouble received message %s",buffer);
    else return(1);
  }
 close(sock);
 } 
 return(0);
}
int sep_begin_prog(){
 int begin=0;
 getch("sep.begin_par","d",&begin);
 if(begin==1) return(sep_send_msg("record","none"));
 return(0);
}
int sep_end_prog(){
 int end;
 end=0; 
 getch("sep.end_par","d",&end);
 if(end==1) return(sep_send_msg("finish","none"));
 return(0);
}
int sep_progress_prog(char *string){
 return(sep_send_msg("progress",string));
 return(0);
}
int sep_exit_prog(){
 return(sep_send_msg("exit","none"));
 return(0);
}
int sep_err_prog(){
 return(sep_send_msg("error","none"));
 return(0);
}
int sep_prog_stat(char *msg, int i, int n,int isend)
{
char temp_ch[1204];
sprintf(temp_ch,"%s finished=%d of %d",msg,i,n);
fprintf(stderr,"%s\n",temp_ch);
if(isend==1) return(sep_progress_prog(temp_ch));
else return(0);
}
