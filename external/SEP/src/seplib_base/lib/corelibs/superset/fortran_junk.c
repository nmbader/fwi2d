#include <sep_fortran_internal.h>
#include <sep_main_external.h>
#include "sep3d.h"
#include "superset.h"

#define UNCONST (char *)

int sep3d_grab_header_vals_sf(const char *tag,char *tag2, void *vals){
  return sep3d_grab_header_vals_s(UNCONST tag,tag2,vals);
}
int sep3d_grab_header_vals_si(const char *tag,char *tag2, void *vals){
  return sep3d_grab_header_vals_s(UNCONST tag,tag2,vals);
}

int sep3d_grab_header_vals_if(const char *tag,int tag2, void *vals){
  return sep3d_grab_header_vals_i(UNCONST tag,tag2,vals);
}
int sep3d_grab_header_vals_ii(const char *tag,int tag2, void *vals){
  return sep3d_grab_header_vals_i(UNCONST tag,tag2,vals);
}

int sep3d_set_header_vals_sf(const char *tag,char *tag2, void *vals){
  return sep3d_set_header_vals_s(UNCONST tag,tag2,vals);
}
int sep3d_set_header_vals_si(const char *tag,char *tag2, void *vals){
  return sep3d_set_header_vals_s(UNCONST tag,tag2,vals);
}

int sep3d_set_header_vals_if(const char *tag,int tag2, void *vals){
  return sep3d_set_header_vals_i(UNCONST tag,tag2,vals);
}
int sep3d_set_header_vals_ii(const char *tag,int tag2, void *vals){
  return sep3d_set_header_vals_i(UNCONST tag,tag2,vals);
}



int sep3d_ritecf(const char *tag, char *t2, int *n, int *f, int *j, void *v, int i1, int i2, int i3,int i4){
  return sep3d_rite(UNCONST tag,t2,n,f,j,v,i1,i2,i3,i4);
}
int sep3d_riteff(const char *tag, char *t2, int *n, int *f,int *j, void *v, int i1, int i2, int i3,int i4){
  return sep3d_rite(UNCONST tag,t2,n,f,j,v,i1,i2,i3,i4);
}
int sep3d_riteif(const char *tag, char *t2, int *n, int *f,int *j, void *v, int i1, int i2, int i3,int i4){
  return sep3d_rite(UNCONST tag,t2,n,f,j,v,i1,i2,i3,i4);
}
int sep3d_read_datacf(const char *tag, char *t2, int i1, int i2, int i3, void *v){
  return sep3d_read_data(UNCONST tag,t2,i1,i2,i3,v);
}
int sep3d_read_dataff(const char *tag, char *t2, int i1, int i2, int i3, void *v){
  return sep3d_read_data(UNCONST tag,t2,i1,i2,i3,v);
}
int sep3d_read_dataif(const char *tag, char *t2, int i1, int i2, int i3, void *v){
  return sep3d_read_data(UNCONST tag,t2,i1,i2,i3,v);
}

