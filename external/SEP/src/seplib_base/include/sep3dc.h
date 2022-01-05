#ifndef SEP3DC_HEAD
#define SEP3DC_HEAD 1
#include<superset.h>
#ifndef SEPNULL2
#define SEPNULL2 ((void *) 0)
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct _sep3d{
char tag[256]    ;  /*sep_3d name */
int      *n;         /*size of dataset */
float   *o,*d;     /*other axes descriptors for regular portion*/
char   **label,**unit;  /*other axes descriptors for regular portion*/
int     nkeys;        /*number of keys in dataset*/
int     ndims;          /*number of dimensions in dataset */
int     drn;        /*data record number for the dataset*/
long long     ntraces;        /*number of traces in dataset*/
char   **keyname,**keytype,**keyfmt ; /*key name, type, and format */
char usage[256],data_format[256],file_format[256];
int *nwind,*fwind, *jwind; /*window parameters */
char tag0[256];
char first[10];
};
typedef struct _sep3d sep3d;

int init_sep3d_struct(sep3d sep3din,sep3d* sep3dout,char *usage);
int sep3d_grab_sep3d(char *sep3dt, sep3d *sep3dc);
int sep3dc_set_drn(sep3d *sep3din,int *drn);
int sep3dc_grab_drn(sep3d *sep3din,int *drn);
int sep3d_set_sep3d(sep3d* sep3dc);
int init_sep3d_tag(char *tag,sep3d* sep3dc,char *usage);
int init_sep3d_par(sep3d *sep3dc,char *usage,char *data_type,char *file_type,int ndim,int nkeys);
int valid_structure(sep3d *sep3dc);
void sep3d_initialize(sep3d *sep3dc);
int sep3d_key_allocate(sep3d *sep3dc,int nkey);
int sep3d_axes_allocate(sep3d *sep3dc,int ndim);
void sep3d_clean(sep3d *sep3dc);
int sep3dc_set_header_block(sep3d *sep3dc,void * block);
int sep3dc_grab_header_block(sep3d *sep3dc,void * block);
int sep3dc_delete(sep3d *sep3dc);
int sep3dc_grab_headers(char *tag, sep3d *sep3dc, int *nh, int *nwind, int *fwind, int *jwind);
int sep3dc_read_data(char *tag, sep3d *sep3dc, char *data,int nt, int ft, int jt);
int  sep3dc_write_description(char *tag, sep3d *sep3dc);
int sep3dc_set_write_status(sep3d *sep3dc,int data,int headers);
int sep3d_rite_num_traces(char *tag,sep3d *sep3dc);
int sep3dc_write_data(char *tag,sep3d *sep3dc,char *data,int *nwind,int *fwind,
int *jwind,int nh,int write_headers,int write_grid);
int print_sep3dc(sep3d *sep3dc);
int sep3dc_ndims(sep3d *sep3dc);
int sep3dc_conform(sep3d *s1, sep3d *s2);
int sep3dc_ge_space(sep3d *s1, sep3d *s2);
int sep3dc_key_index(sep3d *sep3dc, char *name);
int sep3dc_axis_index(sep3d *sep3dc, char *name);
int sep3dc_grab_key_vals(sep3d *sep3dc, char *name, float *header);
int sep3dc_grab_key_vali(sep3d *sep3dc, int num, float *header);
int sep3dc_set_key_vals(sep3d *sep3dc, char *name, float *header);
int sep3dc_set_vali(sep3d *sep3dc, int num, float *header);
int sep3dc_set_grid_values(sep3d *sep3dc, int *vals);
int sep3dc_grab_grid_values(sep3d *sep3dc, int *vals);
int sep3dc_grid_copy(sep3d *s1, sep3d *s2);
int sep3dc__copy(sep3d *s1, sep3d *s2);
int sep3dc_drn_copy(sep3d *s1, sep3d *s2);
int sep3dc_header_copy(sep3d *s1, sep3d *s2);
int sep3dc_set_number_headers(sep3d *sep3dc, int num);
int sep3dc_grab_inorder(sep3d *sep3dc,int *order);
int sep3dc_inorder(sep3d *sep3dc);
int sep3dc_reshape(sep3d *sep3dc,int ndim, int *n);
int sep3dc_collect_headers(sep3d *sep3din, sep3d *sep3dout,int *nh);
int sep3dc_collect_data(sep3d *sep3din, sep3d *sep3dout,char *data_in, char *data_out);
int sep3dc_distribute_headers(sep3d *sep3din, sep3d *sep3dout);
int sep3dc_distribute_data(sep3d *sep3din, sep3d *sep3dout,char *data_in, char *data_out);
int  sep3dc_rite_file_stat(sep3d *sep3din,int idat,int ihead);
int sep3dc_coord_copy(sep3d *s1, sep3d *s2);
int sep3dc_update_ntraces(sep3d *sep3dc);
int  sep3dc_own(sep3d *sep3din,int i);
int sep3dc_close_tag(char *tag, sep3d *sep3din);
int  sep3dc_get_esize(sep3d *sep3din);
int  sep3dc_write_local_description(char *tag, int ith,sep3d *sep3dc);
int sep3dc_broadcast_headers(sep3d *sep3dc, int ifrom);
int sep3dc_broadcast_data(sep3d *sep3dc, int ntr, int from, char *data);
int sep3dc_pass_headers(sep3d *sep3dc, int ifrom,int ito);
int sep3dc_compress_data(sep3d *sep3din, sep3d  *sep3dout,char *in, char *out);
int sep3dc_grab_nh(sep3d *sep3dc, int *nh);                                
int sep3dc_clear(sep3d *sep3dc);
int  sep3dc_grab_coord_vals(sep3d *s1, int *coords);
int  sep3dc_set_coord_vals(sep3d *s1, int *coords);
int  sep3dc_grab_coordh(sep3d *s1, long long *coords);
int  sep3dc_set_coordh(sep3d *s1, long long *coords);
int  sep3dc_set_ncoord(sep3d *s1, int ncoord);

#ifdef __cplusplus
}
#endif

#endif
