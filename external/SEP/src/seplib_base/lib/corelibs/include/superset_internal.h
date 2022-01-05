#ifndef SUPERSET_INT_I  /*SUPERSET_INT_I*/
#define SUPERSET_INT_I aaa
#include<sep3d.h>
#include<../include/superset_struct.h>

#define SEPNULL ((void *) 0)
#if NeedFunctionPrototypes  /*begin proto*/
_XFUNCPROTOBEGIN
sep_3d* tag_info_sep3d(const char *name, enum usage_type type);
void sep3d_print( sep_3d *info );
sep_3d* sep3d_head(void);
void sep3d_addstart( sep_3d *curr );
void sep3d_addend( sep_3d *curr );
sep_3d *sep3d_new( const char *tag, enum usage_type usage );
#ifdef SU_SUPPORT /*begin SU_SUPPORT*/
#include <su.h>
#include <segy.h>
#include<header.h>
#include<hdr.h>
extern int finish_susep(void);
void both_initargs(int, char **);
extern int tgettr(char *, segy *);
extern int susep_read_trace_block(char*,int,int);
extern int susep_rite_trace_block(char*);
extern int tgettra(char *, segy *, int);
extern int tputttra(char *, segy *);
extern int susep_input_init(char *);
extern int susep_output_init(char*,segy*);
#endif /*end SU_SUPPORT*/
extern int sep3d_check_add_keys(char*);
extern int calc_additional_headers(sep_3d*,int*,int*,int*);
extern int sep3d_del_axes(char *);
extern int sep3d_set_wind(char *sep3dname,int *nwind, int *fwind, int *jwind);
extern int sep3d_grab_wind(char *sep3dname,int *nwind, int *fwind, int *jwind);
extern int sep3d_del_mpi_info(char *sep3dname);
extern int sep3d_pass_headers(int ithread,char *sep3name, int ifrom, int ito);

extern int sep3d_send_int(int ithread,int *val, int from, int to);
extern int sep3d_pass_data(int ithread,int ifrom, int ito, int esize, int nelem, char *data);
extern int sep3d_move_part(char *tag,sep_3d *info, int in_all, int out_all, int nmove, int *move_list,char **tags_in,int *thread_from);
extern int sep3d_header_copy(sep_3d *input, sep_3d *output);
extern int sep3d_drn_copy(char *input, char *output);
extern int sep3d_grid_copy(sep_3d *input, sep_3d *output);

extern int sep3d_grab_impi(char *tag,int *impi);
extern int sep3d_set_header_block(char *tag,void *block);
extern int sep3d_grab_nmpi(char *tag,int *nmpi);
extern int sep3d_find_unique_names(int nlist, char **list, int *norig, char**orig, int *mapping,int istart);
extern int sep3d_pass_string(char *string,int from,int ito);
extern int sep3d_broadcast_string(char *string,int from);
extern int sep3d_set_coord_vals(char *sep3dname,int *coords);
extern int sep3d_grab_coord_vals(char *sep3dname,int *coords);
extern int sep3d_set_coordh(char *sep3dname,long long *coords);
extern int sep3d_grab_coordh(char *sep3dname,long long *coords);
extern int sep3d_nconvert(char *sep3dname);
extern int sep3d_convert(char *sep3dname,int *off, char *buf);
extern int sep3d_unconvert(char *sep3dname,int *off,char *buf);
extern int sep3d_buf_size(char *sep3dname,int *nwind, int *fwind, int *jwind, int *n);





extern int sep3d_rite_local(char *tag, char *sep3dname, int *nwind, int *fwind,
  int *jwind,char *data, int ntraces,int write_data, int write_headers,
   int write_grid);
extern int SEP3D_nconvert(sep_3d  *info);
extern int SEP3D_unconvert(sep_3d  *info,int *offset,char *buf);
extern int SEP3D_convert(sep_3d  *info,int *offset,char *buf);
extern int sep3d_broadcast_ints(int *val, int nsz, int ifrom);
extern int sep3d_broadcast_headers(char *sep3dname, int ifrom);
int sep3d_col_headers(char *sep3din, char *sep3dout, int ithread,int *nh);
int sep3d_dist_headers(char *sep3din, char *sep3dout, int ithread);
int sep3d_col_data(char *sep3din, char *sep3dout, int ithread, char *data_in, char *data_out);
int sep3d_distribute_data(char *sep3din, char *sep3dout, int ithread, char *data_in, float *data_out);
int SEP3D_transfer_data(sep_3d *info_from, sep_3d *info_to, char *data_in, char *data_out);
int sep3d_trans_data(char *sep3din, char *sep3dout,  char *data_in, char *data_out);
int sep3d_com_data(char *sep3din, char *sep3dout, int ithread, char *data_in, char *data_out);

                                                                                







extern int SEP3D_get_esize(sep_3d *info);





/*internal sep3d functions */
extern int SEP3D_del_axes(sep_3d *info);
extern int SEP3D_set_ndims(sep_3d *info,int ndims);
extern int SEP3D_grab_ndims(sep_3d *info, int *ndims);
extern int SEP3D_set_axis(sep_3d *info, int axis, int n, float o, float d,char *label, char *unit);
extern int SEP3D_grab_axis(sep_3d *info, int axis, int *n, float *o, float *d,char *label, char *unit);
extern int SEP3D_grab_nod(sep_3d *info, int *n, float *o, float *d);
extern int SEP3D_grab_axes(sep_3d *info, int *n, float *o, float *d,char **label, char **unit);
extern int SEP3D_set_axes(sep_3d *info, int *n, float *o, float *d,char **label, char **unit);
extern int SEP3D_grab_usage(sep_3d *info,char *usage);
extern int SEP3d_grab_file_type(sep_3d *info,char *file_type);
extern int SEP3D_set_file_type(sep_3d *info,char *file_type);
extern int SEP3D_grab_data_type(sep_3d *info,char *data_type);
extern int SEP3D_set_data_type(sep_3d *info,char *data_type);
extern int SEP3D_grab_ntraces(sep_3d *info,long long *ntraces);
extern int SEP3D_set_ntraces(sep_3d *info,long long ntraces);
extern int SEP3D_alloc_coord(sep_3d *info, int ncoord);
extern int sep3d_alloc_coord(char *tag, int ncoord);
extern int SEP3D_delete_coord(sep_3d *info);
extern int SEP3D_coord_copy(sep_3d *input, sep_3d *output);
extern void sep_mpi_stop(void);
void SEP3D_del( sep_3d *curr );
int SEP3D_clean( sep_3d *curr );





extern int SEP3D_wind_coords(sep_3d *info,int *nsz);
extern int SEP3D_rite(char *tag, sep_3d *info, void *data, int write_data,int write_headers, int write_grid);
extern long long SEP3D_count_ntraces(sep_3d *info);
extern int sep3d_count_ntraces(char *tag);
extern int sep3d_set_inorder(const char *sep3dname);
extern int sep3d_unset_inorder(const char *sep3dname);
extern int sep3d_grab_inorder(const char *sep3dname,int *inorder);
extern int sep3d_fileexist(char *tag ,int ilocal);
extern int SEP3D_change_dims(sep_3d *info,int naxes, int *axisindex);
extern int sep3d_rite_format_ith(char  *tagout,char *sep3dname, int ith);
extern int SEP3D_copy_space(sep_3d *in, sep_3d *out);
extern int sep3d_grab_long_ntraces(char  *tag, long long*);
extern int sep3d_set_long_ntraces(char  *tag, long long);



                                                                                        




_XFUNCPROTOEND
#endif
#endif  /*end SUPERSET_INT_I*/
