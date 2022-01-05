#ifndef _SEP3D_H
#define _SEP3D_H

/* Include file <sep3d.h>  This file contains function definitions for
 * the sep3d routines. It should be included in all the c-files in
 * that library.
 *
 * I should probably end up adding this to include.dir/sepcube.h
 *
 */

#include <sepcube.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern int sep_put_number_keys(const char* tag, int* nkeys);
extern int sep_get_number_keys(const char* tag, int* nkeys);
extern int sep_tag_is_pipe(char *tag);
extern void sep_3d_close(void);
extern void init_3d(void);


extern int sep_put_key(const char* tag, const char* key_name, const char* key_type, const char* key_fmt, int* key_index);
extern int sep_get_key_index(const char* tag, char* key_name, int* key_index);
extern int sep_get_key_name(const char* tag, int* key_index, char* key_name);
extern int sep_get_key_type(const char* tag, int* key_index, char* key_type);
extern int sep_get_key_fmt(char* tag, int* key_index, char* key_fmt);
extern int sep_get_key_first_byte(char* tag, int* key_index, int* key_byte);
extern int sep_get_header_bytes(char* tag, int* n_bytes);
extern int sep_copy_grid(char* in_tag, char* out_tag);

extern int sep_put_header_format_tag(const char* tag_history, char* tag_header);
extern int sep_get_header_format_tag(const char* tag_history, char** tag_header);
extern int fget_header_format_tag(char* tag_history, char* tag_header);
extern int fget_grid_format_tag(char* tag_history, char* tag_grid);
extern int sep_copy_hff(char* tag_in, char* tag_out);
extern int sep_copy_gff(char* tag_in, char* tag_out);
extern int sep_set_no_headers(char* tag_history);
extern int sep_set_regular_grid(char* tag_history);
extern int sep_set_no_grid(char* tag_history);

extern int sep_put_grid_format_tag(const char *tag_history, char *tag_grid);
extern int sep_get_grid_format_tag(const char *tag_history, char **tag_grid);

extern int sep_put_data_axis_par(const char *tag_history,int *i_axis,int *n,float *o,float *d,const char *label);
extern int sep_get_data_axis_par(const char *tag_history ,int *i_axis,int *n,float *o,float *d,char *label);
extern int sep_get_number_data_axes(const char *tag_histor,int *n_axes);

extern int sep_put_header_axis_par(const char *tag_history,int *i_axis,int *n,float *o,float *d,const char *label );
extern int sep_get_header_axis_par(const char *tag_history,int *i_axis,int *n,float *o,float *d,char *label);
extern int sep_get_number_header_axes(const char *tag_history,int *n_axis);

extern int sep_put_grid_axis_par(const char *tag_history,int *i_axis,int *n,float *o,float *d,const char *label );
extern int sep_get_grid_axis_par(const char *tag_history,int *i_axis,int *n,float *o,float *d,char *label);
extern int sep_get_number_grid_axes(const char *tag_history,int *n_axis);

extern int sep_put_val_by_index(char *tag_history, int *record_number , int *key_index , int *n_values, void *values);
extern int sep_get_val_by_index(char *tag_history, int *record_number, int *key_index, int *n_values, void *values);
extern int sep_put_val_by_name(char *tag_history, int *record_number, char *key_name, int *n_values, void *values);
extern int sep_get_val_by_name(char *tag_history, int *record_number, char *key_name, int *n_values, void *values);
extern int sep_put_val_headers(char *tag_history, int *record_number, int *n_headers, void *header_values);
extern int sep_get_val_headers(const char *tag_history, int *record_number, int *n_headers, void *header_values);
extern int sep_insert_val_by_index(char *tag_history, int *key_index, int *n_values, void *values, void *all_values);
extern int sep_extract_val_by_index(char *tag_history, int *key_index, int *n_values, void *values, void *all_values);

extern int sep_put_grid_window(char *tag_history, int *n_dim_grid, int *n_grid, int *n_wind, int *f_wind, int *j_wind, int *header_numbers);
extern int sep_get_grid_window(char *tag_history, int *n_dim_grid, int *n_grid, int *n_wind, int *f_wind, int *j_wind, int *header_numbers);

extern int sep_copy_data_pointer(char *tag_in, char *tag_out);
extern int sep_copy_header_keys(char *tag_in, char *tag_out);
extern int sep_reorder_data(char *tag_in,char *tag_out,int n2h,int tsize,int *order);
extern int sep_reorder_data_fast(char *tag_in,char *tag_out,int n2h,int tsize,int *order,int megabyte);

_XFUNCPROTOEND
#else
extern void sep_3d_close();
extern void init_3d();
extern int sep_put_number_keys();
extern int sep_get_number_keys();
extern int stdout_pipe();

extern int sep_put_key();
extern int sep_get_key_index();
extern int sep_get_key_name();
extern int sep_get_key_type();
extern int sep_get_key_fmt();
extern int sep_get_key_first_byte();
extern int sep_get_header_bytes();

extern int sep_put_header_format_tag();
extern int sep_get_header_format_tag();
extern int fget_header_format_tag();
extern int fget_grid_format_tag();
extern int sep_copy_hff();
extern int sep_copy_gff();
extern int sep_set_no_headers();
extern int sep_set_regular_grid();
extern int sep_set_no_grid();

extern int sep_put_data_axis_par();
extern int sep_get_data_axis_par();
extern int sep_get_number_data_axes();
extern int sep_put_header_axis_par();
extern int sep_get_header_axis_par();
extern int sep_get_number_header_axes();
extern int sep_put_grid_axis_par();
extern int sep_get_grid_axis_par();
extern int sep_get_number_grid_axes();


extern int sep_put_val_by_index();
extern int sep_get_val_by_index();
extern int sep_put_val_by_name();
extern int sep_get_val_by_name();
extern int sep_put_val_headers();
extern int sep_get_val_headers();
extern int sep_insert_val_by_index();
extern int sep_extract_val_by_index();

extern int sep_put_grid_window();
extern int sep_get_grid_window();

extern int sep_copy_data_pointer();
extern int sep_copy_header_keys();
extern int sep_reorder_data();
#endif

#ifdef __cplusplus
}
#endif
#endif/* _SEP3D_H */
/*  $Id: sep3d.h,v 1.1.1.1 2004/03/25 06:37:22 cvs Exp $ */
