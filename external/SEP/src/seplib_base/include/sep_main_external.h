#ifndef SEP_MAIN_EXTERNAL_H
#define SEP_MAIN_EXTERNAL_H
#include <prototypes.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern char* get_format( char *tag );
extern char *expandnm(char*,char*);
extern void grab_history(const char *tag,char *buf, int nmax,int *nsize);
extern char *datapath(char*);
extern int auxclose(const char *);
extern int auxhclose(const char *);
extern int aux_unlink(const char *);
extern FILE *auxin(const char*);
extern int fauxin(const char*);
extern FILE *auxout(const char*);
extern int copy_history( const char*, const char*);
extern FILE *auxinout(const char*);
extern FILE *auxsockout(const char*);
extern FILE *auxtmp(const char*);
extern FILE *auxscr(const char*);
extern int doc(char*)  ;
extern void sep_add_doc_line (const char* line);
extern int hclose(void)  ;
extern int isapipe(int);
extern int make_unpipe(char*);
extern int noheader(void);
extern int redin(void);
extern int redout(void);
extern int sepwarn(int wrn, const char *format, ... );
extern int sep_prog( char *arg);
extern int separg(int iarg, char *arg);
extern int seperr( const char *format, ... );
extern int sreed_window_new(const char*,const int,const int*,const int*,const int*,const int*,const int,void*);
extern int srite_window_new(const char*,const int,const int*,const int*,const int*,const int*,const int,void*);
extern int sreed_window(const char*,const int*,const int*,const int*,const int*,const int*,const int,void*);
extern int srite_window(const char*,const int*,const int*,const int*,const int*,const int*,const int,void*);
extern int sep_window (const int,const char*,const int*, const int*,const int*,const int*,const int*,const int, void *);
extern int snap( char*, int, int, void* );
extern void set_format( char*, char* );
extern int slice (char*,int,int,int,void*);
extern int sreed_raw(char*,void*,int);
extern int sreed2(char*,void*,int,char*);
extern int srite(const char*, void*, int);
extern size_t sreed_big(const char*, void*,size_t);
extern size_t srite_big(const char*, void*, size_t);
extern int srite_raw(char*,void*,int);
extern int srite2(char*,void*,int,char*);
extern int sreed(const char*, void*, int);
extern long long sseek_l_l(const char*,long long, int);
extern int sseek(const char*,int, int);
extern int sseekable(const char*);
extern int sseek_block(const char*,int,int,int);
double sseek_block_d(const char *,int,int,int);
extern int file_position(const char*,int,int*,int*);
extern int ssize (char*);
extern int ssize_block (char* ,int);
extern void debug_tag (char*);
extern FILE *input(void);
extern FILE *output(void);
extern int hcount(int);
extern FILE *sep_head(void);
extern char *alloc(size_t);
extern int seploc(void*);
extern void initpar_f(char*);
extern void c2h(int*,int*, int,int*);
extern void h2c(int, int*, int,int*);
extern int evaluate_expression(char *exp,int get_val(char *,double*),int nvals ,double *result);
extern double ssize_info (char * tag);
extern int sep_send_msg(char *status,char *extra) ;
extern int sep_end_prog();
extern int sep_exit_prog();
extern int sep_progress_prog(char *fd);
extern int sep_prog_stat(char *fd,int i, int n,int send);
extern int sep_err_prog();
extern int sep_begin_prog();


_XFUNCPROTOEND
#else
extern char *expandnm();
extern char *datapath();
extern int auxclose();
extern int auxhclose();
extern FILE *auxin();
extern FILE *auxout();
extern int copy_history();
extern int aux_unlink();
extern FILE *auxinout();
extern FILE *auxsockout();
extern FILE *auxscr();
extern FILE *auxtmp();
extern int doc();
extern int hclose();
extern int isapipe();
extern int make_unpipe();
extern int noheader();
extern int redin();
extern int redout();
extern int sepwarn();
extern int seperr();
extern int sreed_window();
extern int srite_window();
extern int sep_window ();
extern int snap();
extern void set_format();
extern int slice ();
extern int sreed_raw();
extern int sreed2();
extern int sreed();
extern int srite_raw();
extern int srite2();
extern int srite();
extern int sseek();
extern int sseek_block();
extern dobule sseek_block_d();
extern int ssize ();
extern int ssize_block();
extern FILE *input();
extern FILE *output();
extern int hcount();
extern FILE *sep_head();
extern char *alloc();
extern int seploc();
extern void c2h();
extern void h2c();
extern int file_position();
extern int evaluate_expression();
extern void debug_tag ();
extern int sseekable();
#endif


#ifndef MAX
#define MAX(a,b) ( ((a)>(b)) ? (a):(b) )
#endif
#ifndef MIN
#define MIN(a,b) ( ((a)<(b)) ? (a):(b) )
#endif

#ifdef __cplusplus
}
#endif

#endif
