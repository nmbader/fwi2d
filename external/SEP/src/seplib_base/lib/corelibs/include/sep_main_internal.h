#include <prototypes.h>
#include <stdio.h>
#include <stdlib.h>


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern char *envhead(char*);
extern int findnm( int, char*, int );
extern int fullnm( char*, int );
extern int isatape(int);
extern int isclosed(FILE*);
extern int fdordinary(int);
extern int isordinary(char*);
extern char *maketitle( char*);
extern double sseek_block_d(const char*,int,int,int);
void mkrandom_string(char *string_in, char *string_out);
void sep_headername(char *tag, char *headername);
extern char* sep_tail(char*);
extern char* get_format_name( int num);
int get_data_name(char *tag,char *string_out);


/*fortran junk to be rmoved some day*/
int srite2_f(const char *tag, void *v, const int sz,const char *typ);
int srite2_i(const char *tag, void *v, const int sz,const char *typi);
int srite2_c(const char *tag, void *v, const int sz,const char *typ);
int sreed2_f(const char *tag, void *v, const int sz,const char *typ);
int sreed2_i(const char *tag, void *v, const int sz,const char *typ);
int sreed2_c(const char *tag, void *v, const int sz,const char *typ);

int srite_f(const char *tag, void *v, const int sz);
int srite_i(const char *tag, void *v, const int sz);
int srite_c(const char *tag, void *v, const int sz);
int sreed_f(const char *tag, void *v, const int sz);
int sreed_i(const char *tag, void *v, const int sz);
int sreed_c(const char *tag, void *v, const int sz);

int srite_window_i(const char *tag,  int*, int *ng, int *nw, int *fw, int *jw, int sz,
  void *val);
int srite_window_f(const char *tag,  int*, int *ng, int *nw, int *fw, int *jw, int sz,
  void *val);
int srite_window_c(const char *tag,  int*, int *ng, int *nw, int *fw, int *jw, int sz,
  void *val);
int sreed_window_i(const char *tag,  int*, int *ng, int *nw, int *fw, int *jw, int sz,
  void *val);
int sreed_window_f(const char *tag,  int*, int *ng, int *nw, int *fw, int *jw, int sz,
  void *val);
int sreed_window_c(const char *tag,  int*, int *ng, int *nw, int *fw, int *jw, int sz,
  void *val);

int tetch_f_f(const char *arg, const char *typ, void *val);
int tetch_i_f(const char *arg, const char *typ, void *val);
int tetch_s_f(const char *arg, const char *typ, void *val);
int tetch_l_f(const char *arg, const char *typ, void *val);

int putch_f_f(const char *arg, const char *typ, void *val);
int putch_i_f(const char *arg, const char *typ, void *val);
int putch_s_f(const char *arg, const char *typ, void *val);
int putch_l_f(const char *arg, const char *typ, void *val);

int hetch_f_f_a(const char *arg, const char *typ, void *val);
int hetch_i_f_a(const char *arg, const char *typ, void *val);
int hetch_f_f(const char *arg, const char *typ, void *val);
int hetch_i_f(const char *arg, const char *typ, void *val);
int hetch_s_f(const char *arg, const char *typ, void *val);
int hetch_l_f(const char *arg, const char *typ, void *val);

int getch_f_f_a(const char *arg, const char *typ, void *val);
int getch_i_f_a(const char *arg, const char *typ, void *val);
int getch_f_f(const char *arg, const char *typ, void *val);
int getch_i_f(const char *arg, const char *typ, void *val);
int getch_s_f(const char *arg, const char *typ, void *val);
int getch_l_f(const char *arg, const char *typ, void *val);

int fetch_f_f_a(const char *arg, const char *typ, void *val);
int fetch_i_f_a(const char *arg, const char *typ, void *val);
int fetch_f_f(const char *arg, const char *typ, void *val);
int fetch_i_f(const char *arg, const char *typ, void *val);
int fetch_s_f(const char *arg, const char *typ, void *val);
int fetch_l_f(const char *arg, const char *typ, void *val);

int auxpar_f_f_a(const char *arg, const char *typ, void *val,const char *tag);
int auxpar_i_f_a(const char *arg, const char *typ, void *val,const char *tag);
int auxpar_f_f(const char *arg, const char *typ, void *val,const char *tag);
int auxpar_i_f(const char *arg, const char *typ, void *val,const char *tag);
int auxpar_s_f(const char *arg, const char *typ, void *val,const char *tag);
int auxpar_l_f(const char *arg, const char *typ, void *val,const char *tag);

int auxputch_f_f_a(const char *arg, const char *typ, void *val,const char *tag);
int auxputch_i_f_a(const char *arg, const char *typ, void *val,const char *tag);
int auxputch_f_f(const char *arg, const char *typ, void *val,const char *tag);
int auxputch_i_f(const char *arg, const char *typ, void *val,const char *tag);
int auxputch_s_f(const char *arg, const char *typ, void *val,const char *tag);
int auxputch_l_f(const char *arg, const char *typ, void *val,const char *tag);

extern int socklisten(  int , int );
extern int opensock2( char* , char* );
extern int opensock1( char*, int);

_XFUNCPROTOEND
#else /*END OF NO PROTO */
extern char *envhead();
extern int findnm();
extern int fullnm();
extern int isatape();
extern int isclosed();
extern int fdordinary();
extern int isordinary();
extern char *maketitle();


extern int socklisten();
extern int opensock2();
extern int opensock1();
#endif /*END OF NO PROTO */

/*  $Id: sep_main_internal.h,v 1.3 2004/04/08 22:32:27 bob Exp $ */
