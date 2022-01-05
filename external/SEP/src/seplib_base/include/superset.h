#ifndef SUPERSET_H
#define SUPERSET_H dfad
#include <sepConfig.h>
#include<sep3d.h>

#ifdef __cplusplus
extern "C" {
#endif

enum{
 /*****************/
 /* ERROR RETURNS */
/*****************/
SUCCESS       = 0,/* successul exection of function */
FAIL_OTHER    =-1,/* failed because of lower level call */
EXISTS        = 1,/* tag has already been declared with other state */
INVALID_USAGE = 2, /* usage passed not acceptable */ 
INVALID_STRUC = 3, /* usage passed not acceptable */ 
INVALID_DATA  = 4, /* usage passed not acceptable */ 
NOT_MET       = 5  /* requirements to run routine not met */
};



#ifndef SEPH2C
#define SEPH2C(index,block,n) (int)(((long long)index/(long long)block)%((long long)n))
#endif
#ifndef SEPC2H
#define SEPC2H(index,block)  ((long long)(index)*((long long)(block)))
#endif

#ifdef SU_SUPPORT
enum{
BEG_TR_DEFAULT=-1,/*the default value for beg_trace */
END_TR_DEFAULT=-2,/*the default value for end_trace */
BEG_TR_INIT   =-3,/*the default value for beg_trace */
EXTRA_KEY_DEFAULT=-1 /*the default value for  extra_keys */
};
#else
#ifndef SET_SDOC
char *sdoc[] = {NULL};
#define SET_SDOC 1
#endif
#endif


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern int sep3d_set_sect_threads(char *sep3dname, int *sect); /*private*/
extern int sep3d_grab_axis_sect_ithread(char *sep3dname,int isect, int *axis);
extern int sep3d_grab_naxis_ithread(char *sep3dname,int isect, int *naxis);
extern int sep3d_read_data(char *tag, char *sep3dname ,int nwind,int fwind,int jwind, char *array);
int sep3d_read_header_section(char *struct_big,char *struct_in,int thread_from, char *struct_out, int thread_to, int isect, int *nwind, int *fwind, int *jwind,int *nh);

extern int sep3d_return_local_size(char *sep3din,int *nin, int *fin, int *jin, int *nout, int *fout, int *jout);
extern int sep3d_read_data_section(char *sect_big, char *sect_in, int thread_from, char *sect_out, int thread_to,
int isect, int nwind, int fwind, int jwind,char *data);
extern int sep3d_read_header_section_list(char *sep3din,int ithread, int *nwind, int *fwind, int *jwind,int nsect, int *sect_list,int *nh);
extern int sep3d_read_header_sections(char *sep3din, int ithread, int *nwind, int *fwind, int *jwind,int *nh);
extern int sep3d_read_data_section_list(char *sep3din, int ithread, int nwind, int fwind, int jwind,int nsect, int *sect_list,char *data);
extern int sep3d_read_data_sections(char *sep3din, int ithread, int nwind, int fwind, int jwind,char *data);
extern int sep3d_set_coord(char *sep3dname,int *coords);
extern int sep3d_grab_coord(char *sep3dname,int *coords);
extern int sep3d_grab_inorder(const char *, int *);
extern int sep3d_set_inorder(const char *);
extern int sep3d_unset_inorder(const char *);
extern int sep3d_read_header_all(char *sep3din,int *nwind, int *fwind, int *jwind, int *nh);
extern int sep3d_read_header_local(char *sep3din,int *nwind, int *fwind, int *jwind, int *nh);
extern int sep3d_rite_all(char *tag, char *sep3dname, int *nwind, int *fwind,
  int *jwind,char *data, int write_data, int write_headers,
   int write_grid,int ntraces);
extern int sep3d_rite_local(char *tag, char *sep3dname,  int *nwind, int *fwind,
  int *jwind,char *data, int write_data, int write_headers,
   int write_grid,int ntraces);
extern int sep3d_rite_headers_section(char *tag, char *sep3dname, int *nwind, int *fwind, int *jwind,char *data, int ntraces, int write_data, int write_headers,
   int write_grid,int nsect, int *sect_list);
extern int sep3d_set_nsect_mpi(char*,int);
extern int sep3d_grab_nsect_mpi(char*,int*);
extern int sep3d_set_mpi_axis(char*,int);
extern int sep3d_grab_mpi_axis(char*,int*);
extern int sep3d_set_distrib_type(char*,char*);
extern int sep3d_grab_distrib_type(char*,char*);
extern int sep3d_set_tag0(char*,char*);
extern int sep3d_grab_tag0(char*,char*);
extern int sep3d_del_mpi(char*);

extern int sep3d_set_thread_tag(char*,int,char*);
extern int sep3d_grab_thread_tag(char*,int,char*);

extern int sep3d_grab_mpi_sections(char*,int*,int*,int*);
extern int sep3d_set_mpi_sections(char*,int,int*,int*);

extern int sep3d_struct_init(char*, char*, char*);
extern int sep3d_par_init(char*, char *);
extern int sep3d_tag_init(char*, char*, char*);
extern int sep3d_set_ndims(char*, int);
extern int sep3d_grab_ndims(char*, int*);
extern int sep3d_grab_nh(char*, int*);
extern int sep3d_set_nkeys(char*, int);
extern int sep3d_grab_nkeys(char*, int*);
extern int sep3d_grab_key(char*, int, char*, char*, char*);
extern int sep3d_set_key(char*, int, char*, char*, char*);
extern int sep3d_set_axis(char*, int,int,float,float,char*,char*);
extern int sep3d_set_axis2(float,float);
extern int sep3d_close_tag(char*,char*);
extern int sep3d_grab_axis(char*, int,int*,float*,float*,char*,char*);
extern int sep3d_conform(char *sep3dname1,char  *sep3dname2);
extern int sep3d_ge_space(char *sep3dname1,char  *sep3dname2);
extern int sep3d_grab_nod(char*, int*,float*,float*);
extern int sep3d_grab_axes(char*, int*,float*,float*,char**,char**);
extern int sep3d_set_axes(char*, int*,float*,float*,char**,char**);
extern int sep3d_grab_keys(char*,char**,char**,char**);
extern int sep3d_set_keys(char*, char**,char**,char**);
extern int sep3d_grab_usage(char*, char*);
extern int sep3d_grab_data_type(char*, char*);
extern int sep3d_set_data_type(char*, char*);
extern int sep3d_print_info(char*);
extern int sep3d_set_big_ntraces(char*,int,int,int);
extern int sep3d_set_ntraces(char*,int ntraces);
extern int sep3d_grab_ntraces(char*,int *ntraces);
extern int sep3d_grab_key_index(char*,char*,int*);
extern int sep3d_grab_axis_index(char*,char*,int*);
extern int sep3d_grab_drn(char*,int*);
extern int sep3d_set_drn(char*,int*);
extern int sep3d_reed_axis(char*,int, int*, float*, float*, char*, char*);
extern int sep3d_rite_axis(char*,int, int, float, float, char*, char*);
extern int sep3d_read_data_from_grid(char*, char*,int*,int*,int*,int,char*);
extern int sep3d_rite(char*,char*,int*,int*,int*,void*,int,int,int,int);
extern int sep3d_rite_list(char*,int ,int,int,int,int,int,int*,char*);
extern int sep3d_grab_header_vals_s(char*,char*,int*);  
extern int sep3d_grab_header_vals_i(char*,int,int*);  
extern int sep3d_set_header_vals_s(char*,char*,int*);  
extern int sep3d_set_header_vals_i(char*,int,int*);  
extern int sep3d_read_headers(char*,char*,int*,int*,int*,int*);
extern int sep3d_set_nh(char*,int);
extern int sep3d_grab_trace_numbers(char *, int*);
extern int sep3d_set_rite_status(char*,int,int);
extern int sep3d_rite_ntraces(char*,char*);
extern int sep3d_rite_format(char*,char*);
extern int sep3d_set_ntraces(char*,int);
extern int sep3d_grab_ntraces(char*,int*);
extern int sep3d_clear_headers(char*);
extern int sep3d_grab_file_type(char*,char*);
extern int sep3d_set_file_type(char*,char*);
extern int sep3d_read_list(char*,int,int,int,int,int,int,int*,char*);
extern int sep3d_status(char *);
extern int* sep3d_grab_header_pointer(char*);
extern int sep3d_set_header_pointer(char*,int,int*);
extern int sep3d_nullify_header_pointer(char*);
extern int sep3d_copy_headers(char*,char*);
extern int sep3d_grab_header_block(char*,void*);
extern int sep3d_clear_grid(char*);
extern int sep3d_set_ng(char*,int);
extern int sep3d_grab_ng(char*,int*);
extern int* sep3d_grab_grid_pointer(char*);
extern int sep3d_set_grid_pointer(char*,int,int*);
extern int sep3d_set_grid_vals(char*,int*);
extern int sep3d_nullify_grid_pointer(char*);
extern int sep3d_copy_grid(char*,char*);
extern int sep3d_grab_grid_block(char*,int*);
extern int sep3d_rite_file_stat(char*,int,int);
extern int sep3d_change_dims(char*,int,int*);
extern int sep3d_delete(char*);
extern int sep_get_non_default_axes(char *tag,int *naxes,int istart);
extern int sep3d_drn_set(char*);
extern int sep3d_write_dff(char *tag, char *sep3dname);
extern int sep3d_copy_struct(char *sep3din, char *sep3dout);
extern int fget_distrib_format_tag(char * tag_history, char *tag_distrib);
extern int sep_get_distrib_format_tag(char* tag_history, char** tag_distrib);
extern int sep_put_distrib_format_tag(char* tag_history, char* tag_dff);
extern int sep3d_sync(char *tag, char *sep3dname);
extern int sep3d_have_dff(char *tag);
extern int sep3d_copy_coords(char *sep3in,char *sep3dout);
extern int sep3d_broadcast_data(int ifrom,  int esize, int nelem, char *data);
extern void sep_mpi_stop(void);




_XFUNCPROTOEND
#endif
#ifdef __cplusplus
}
#endif
#endif
/*  $Id: superset.h,v 1.3 2004/07/08 18:15:32 bob Exp $ */
