#include <sep_main_internal.h>
#include <sep_main_external.h>
#include <sep_pars_external.h>

#define CONSTDECL const
#define UNCONST (char *)
/*this file is only necessary because gfortran has a bug in older versions*/
int sreed2_f(CONSTDECL char *tag,void *v, int sz,CONSTDECL char *typ){
   return sreed2(UNCONST tag,v,sz,UNCONST typ);
}
int sreed2_c(CONSTDECL char *tag,void *v, const int sz,CONSTDECL char *typ){
  return sreed2(UNCONST tag,v,sz,UNCONST typ);
}
int sreed2_i(CONSTDECL char *tag,void *v, const int sz,CONSTDECL char *typ){
  return  sreed2(UNCONST tag,v,sz,UNCONST typ);
}
int srite2_f(CONSTDECL char *tag,void *v, const int sz,CONSTDECL char *typ){
  return srite2(UNCONST tag,v,sz,UNCONST typ);
}
int srite2_c(CONSTDECL char *tag,void *v, const int sz,CONSTDECL char *typ){
  return  srite2(UNCONST tag,v,sz,UNCONST typ);
}
int srite2_i(CONSTDECL char *tag,void *v, const int sz,CONSTDECL char *typ){
  return  srite2(UNCONST tag,v,sz,UNCONST typ);
}
long long sreedll_i(CONSTDECL char *tag,void *v, const long long sz){
  return sreed(UNCONST tag,v,sz);
}
long long sreedll_c(CONSTDECL char *tag,void *v, const long long sz){
  return sreed(UNCONST tag,v,sz);
}
long long sreedll_x(CONSTDECL char *tag,void *v, const long long sz){
  return sreed(UNCONST tag,v,sz);
}
long long sreedll_f(CONSTDECL char *tag,void *v, const long long sz){
  return sreed(UNCONST tag,v,sz);
}
long long sritell_x(CONSTDECL char *tag,void *v, const long long sz){
  return srite(UNCONST tag,v,sz);
}
long long sritell_c(CONSTDECL char *tag,void *v, const long long sz){
  return srite(UNCONST tag,v,sz);
}
long long sritell_i(CONSTDECL char *tag,void *v, const long long sz){
  return srite(UNCONST tag,v,sz);
}
long long sritell_f(CONSTDECL char *tag,void *v, const long long sz){
  return srite(UNCONST tag,v,sz);
}
int sreed_f(CONSTDECL char *tag,void *v, const int sz){
  return sreed(UNCONST tag,v,sz);
}
int sreed_c(CONSTDECL char *tag,void *v, const int sz){
  return sreed(UNCONST tag,v,sz);
}
int sreed_i(CONSTDECL char *tag,void *v, const int sz){
  return  sreed(UNCONST tag,v,sz);
}
int srite_f(CONSTDECL char *tag,void *v, const int sz){
  return  srite(UNCONST tag,v,sz);
}
int srite_c(CONSTDECL char *tag,void *v, const int sz){
  return srite(UNCONST tag,v,sz);
}
int srite_i(CONSTDECL char *tag,void *v, const int sz){
  return  srite(UNCONST tag,v,sz);
}
int srite_window_f(CONSTDECL char *tag, int *nd, int *ng, int *nw, int *fw,int *jw, int sz, void *val){
   return  srite_window(UNCONST tag,nd,ng,nw,fw,jw,sz,val);
}
int srite_window_c(CONSTDECL char *tag, int *nd, int *ng, int *nw, int *fw,int *jw, int sz, void *val){
return   srite_window(UNCONST tag,nd,ng,nw,fw,jw,sz,val);
}
int srite_window_i(CONSTDECL char *tag, int *nd, int *ng, int *nw, int *fw,int *jw, int sz, void *val){
  return srite_window(UNCONST tag,nd,ng,nw,fw,jw,sz,val);
}

int sreed_window_i(CONSTDECL char *tag, int *nd, int *ng, int *nw, int *fw,int *jw, int sz, void *val){
  return  srite_window(UNCONST tag,nd,ng,nw,fw,jw,sz,val);
}
int sreed_window_f(CONSTDECL char *tag, int *nd, int *ng, int *nw, int *fw,int *jw, int sz, void *val){
return   srite_window(UNCONST tag,nd,ng,nw,fw,jw,sz,val);
}
int sreed_window_c(CONSTDECL char *tag, int *nd, int *ng, int *nw, int *fw,int *jw, int sz, void *val){
  return  srite_window(UNCONST tag,nd,ng,nw,fw,jw,sz,val);
}
int tetch_f_f(CONSTDECL char *arg, const char *typ, void *val){return tetch(UNCONST arg,UNCONST typ,val);}
int tetch_g_f(CONSTDECL char *arg, const char *typ, void *val){return tetch(UNCONST arg,UNCONST typ,val);}
int tetch_i_f(CONSTDECL char *arg, const char *typ, void *val){return tetch(UNCONST arg,UNCONST typ,val);}
int tetch_s_f(CONSTDECL char *arg, const char *typ, void *val){return tetch(UNCONST arg,UNCONST typ,val);}
int tetch_l_f(CONSTDECL char *arg, const char *typ, void *val){return tetch(UNCONST arg,UNCONST typ,val);}

int fetch_f_f_a(CONSTDECL char *arg, const char *typ, void *val){return fetch(UNCONST arg,UNCONST typ,val);}
int fetch_i_f_a(CONSTDECL char *arg, const char *typ, void *val){return fetch(UNCONST arg,UNCONST typ,val);}
int fetch_f_f(CONSTDECL char *arg, const char *typ, void *val){return fetch(UNCONST arg,UNCONST typ,val);}
int fetch_g_f(CONSTDECL char *arg, const char *typ, void *val){return fetch(UNCONST arg,UNCONST typ,val);}
int fetch_l_f(CONSTDECL char *arg, const char *typ, void *val){return fetch(UNCONST arg,UNCONST typ,val);}
int fetch_s_f(CONSTDECL char *arg, const char *typ, void *val){return fetch(UNCONST arg,UNCONST typ,val);}
int fetch_i_f(CONSTDECL char *arg, const char *typ, void *val){return fetch(UNCONST arg,UNCONST typ,val);}

int getch_f_f_a(CONSTDECL char *arg, const char *typ, void *val){return getch(UNCONST arg,UNCONST typ,val);}
int getch_g_f_a(CONSTDECL char *arg, const char *typ, void *val){return getch(UNCONST arg,UNCONST typ,val);}
int getch_i_f_a(CONSTDECL char *arg, const char *typ, void *val){return getch(UNCONST arg,UNCONST typ,val);}
int getch_f_f(CONSTDECL char *arg, const char *typ, void *val){return getch(UNCONST arg,UNCONST typ,val);}
int getch_g_f(CONSTDECL char *arg, const char *typ, void *val){return getch(UNCONST arg,UNCONST typ,val);}
int getch_i_f(CONSTDECL char *arg, const char *typ, void *val){return getch(UNCONST arg,UNCONST typ,val);}
int getch_s_f(CONSTDECL char *arg, const char *typ, void *val){return getch(UNCONST arg,UNCONST typ,val);}
int getch_l_f(CONSTDECL char *arg, const char *typ, void *val){return getch(UNCONST arg,UNCONST typ,val);}

int hetch_f_f_a(CONSTDECL char *arg, const char *typ, void *val){return hetch(UNCONST arg,UNCONST typ,val);}
int hetch_i_f_a(CONSTDECL char *arg, const char *typ, void *val){return hetch(UNCONST arg,UNCONST typ,val);}
int hetch_l_f_a(CONSTDECL char *arg, const char *typ, void *val){return hetch(UNCONST arg,UNCONST typ,val);}
int hetch_g_f_a(CONSTDECL char *arg, const char *typ, void *val){return hetch(UNCONST arg,UNCONST typ,val);}
int hetch_f_f(CONSTDECL char *arg, const char *typ, void *val){return hetch(UNCONST arg,UNCONST typ,val);}
int hetch_l_f(CONSTDECL char *arg, const char *typ, void *val){return hetch(UNCONST arg,UNCONST typ,val);}
int hetch_g_f(CONSTDECL char *arg, const char *typ, void *val){return hetch(UNCONST arg,UNCONST typ,val);}
int hetch_i_f(CONSTDECL char *arg, const char *typ, void *val){return hetch(UNCONST arg,UNCONST typ,val);}
int hetch_s_f(CONSTDECL char *arg, const char *typ, void *val){return hetch(UNCONST arg,UNCONST typ,val);}

int putch_f_f_a(CONSTDECL char *arg, const char *typ, void *val){return putch(UNCONST arg,UNCONST typ,val);}
int putch_i_f_a(CONSTDECL char *arg, const char *typ, void *val){return putch(UNCONST arg,UNCONST typ,val);}
int putch_f_f(CONSTDECL char *arg, const char *typ, void *val){return putch(UNCONST arg,UNCONST typ,val);}
int putch_g_f(CONSTDECL char *arg, const char *typ, void *val){return putch(UNCONST arg,UNCONST typ,val);}
int putch_i_f(CONSTDECL char *arg, const char *typ, void *val){return putch(UNCONST arg,UNCONST typ,val);}
int putch_s_f(CONSTDECL char *arg, const char *typ, void *val){return putch(UNCONST arg,UNCONST typ,val);}
int putch_l_f(CONSTDECL char *arg, const char *typ, void *val){return putch(UNCONST arg,UNCONST typ,val);}
 
int auxpar_f_f_a(const char *arg, const char *typ, void *val,const char *fle){return auxpar(UNCONST arg,typ,val,fle);}
int auxpar_g_f_a(const char *arg, const char *typ, void *val,const char *fle){return auxpar(UNCONST arg,typ,val,fle);}
int auxpar_i_f_a(const char *arg, const char *typ, void *val,const char *fle){return auxpar(UNCONST arg,typ,val,fle);}
int auxpar_f_f(const char *arg, const char *typ, void *val,const char *fle){return auxpar(UNCONST arg,typ,val,fle);}
int auxpar_g_f(const char *arg, const char *typ, void *val,const char *fle){return auxpar(UNCONST arg,typ,val,fle);}
int auxpar_i_f(const char *arg, const char *typ, void *val,const char *fle){return auxpar(UNCONST arg,typ,val,fle);}
int auxpar_l_f(const char *arg, const char *typ, void *val,const char *fle){return auxpar(UNCONST arg,typ,val,fle);}
int auxpar_s_f(const char *arg, const char *typ, void *val,const char *fle){return auxpar(UNCONST arg,typ,val,fle);}

int auxputch_f_f_a(const char *arg, const char *typ, void *val,const char *fle){return auxputch(UNCONST arg,typ,val,fle);}
int auxputch_g_f_a(const char *arg, const char *typ, void *val,const char *fle){return auxputch(UNCONST arg,typ,val,fle);}
int auxputch_i_f_a(const char *arg, const char *typ, void *val,const char *fle){return auxputch(UNCONST arg,typ,val,fle);}
int auxputch_f_f(const char *arg, const char *typ, void *val,const char *fle){return auxputch(UNCONST arg,typ,val,fle);}
int auxputch_g_f(const char *arg, const char *typ, void *val,const char *fle){return auxputch(UNCONST arg,typ,val,fle);}
int auxputch_i_f(const char *arg, const char *typ, void *val,const char *fle){return auxputch(UNCONST arg,typ,val,fle);}
int auxputch_s_f(const char *arg, const char *typ, void *val,const char *fle){return auxputch(UNCONST arg,typ,val,fle);}
int auxputch_l_f(const char *arg, const char *typ, void *val,const char *fle){return auxputch(UNCONST arg,typ,val,fle);}

