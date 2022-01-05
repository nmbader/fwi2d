 
#include "cfortran.h"
#include <sep3d.h>

/*
  Modified: Paul hargrove 11/06/96, Added alternate FORTRAN entry points
             for all functions that contain an undescore for compatibilty
             with g77's behavior of placing a double underscore at the
             end of an external name that contains underscores.
*/

/* "Normal" stubs */
FCALLSCFUN2(INT,sep_put_number_keys,SEP_PUT_NUMBER_KEYS,sep_put_number_keys,STRING,PINT)
FCALLSCFUN2(INT,sep_get_number_keys,SEP_GET_NUMBER_KEYS,sep_get_number_keys,STRING,PINT)

FCALLSCFUN5(INT,sep_put_key,SEP_PUT_KEY,sep_put_key,STRING,STRING,STRING,STRING,PINT)
FCALLSCFUN3(INT,sep_get_key_index,SEP_GET_KEY_INDEX,sep_get_key_index,STRING,STRING,PINT)
FCALLSCFUN3(INT,sep_get_key_name,SEP_GET_KEY_NAME,sep_get_key_name,STRING,PINT,PSTRING)
FCALLSCFUN3(INT,sep_get_key_type,SEP_GET_KEY_TYPE,sep_get_key_type,STRING,PINT,PSTRING)
FCALLSCFUN3(INT,sep_get_key_fmt,SEP_GET_KEY_FMT,sep_get_key_fmt,STRING,PINT,PSTRING)
FCALLSCFUN3(INT,sep_get_key_first_byte,SEP_GET_KEY_FIRST_BYTE,sep_get_key_first_byte,STRING,PINT,PINT)
FCALLSCFUN2(INT,sep_get_header_bytes,SEP_GET_HEADER_BYTES,sep_get_header_bytes,STRING,PINT)

FCALLSCFUN2(INT,sep_put_header_format_tag,SEP_PUT_HEADER_FORMAT_TAG,sep_put_header_format_tag,STRING,STRING)
FCALLSCFUN2(INT,fget_header_format_tag,SEP_GET_HEADER_FORMAT_TAG,sep_get_header_format_tag,STRING,PSTRING)

FCALLSCFUN2(INT,sep_put_grid_format_tag,SEP_PUT_GRID_FORMAT_TAG,sep_put_grid_format_tag,STRING,STRING)
FCALLSCFUN2(INT,fget_grid_format_tag,SEP_GET_GRID_FORMAT_TAG,sep_get_grid_format_tag,STRING,PSTRING)

FCALLSCFUN6(INT,sep_put_data_axis_par,SEP_PUT_DATA_AXIS_PAR,sep_put_data_axis_par,STRING,PINT,PINT,PFLOAT,PFLOAT,STRING)
FCALLSCFUN6(INT,sep_get_data_axis_par,SEP_GET_DATA_AXIS_PAR,sep_get_data_axis_par,STRING,PINT,PINT,PFLOAT,PFLOAT,PSTRING)
FCALLSCFUN2(INT,sep_get_number_data_axes,SEP_GET_NUMBER_DATA_AXES,sep_get_number_data_axes,STRING,PINT)

FCALLSCFUN6(INT,sep_put_header_axis_par,SEP_PUT_HEADER_AXIS_PAR,sep_put_header_axis_par,STRING,PINT,PINT,PFLOAT,PFLOAT,STRING)
FCALLSCFUN6(INT,sep_get_header_axis_par,SEP_GET_HEADER_AXIS_PAR,sep_get_header_axis_par,STRING,PINT,PINT,PFLOAT,PFLOAT,PSTRING)
FCALLSCFUN2(INT,sep_get_number_header_axes,SEP_GET_NUMBER_HEADER_AXES,sep_get_number_header_axes,STRING,PINT)

FCALLSCFUN6(INT,sep_put_grid_axis_par,SEP_PUT_GRID_AXIS_PAR,sep_put_grid_axis_par,STRING,PINT,PINT,PFLOAT,PFLOAT,STRING)
FCALLSCFUN6(INT,sep_get_grid_axis_par,SEP_GET_GRID_AXIS_PAR,sep_get_grid_axis_par,STRING,PINT,PINT,PFLOAT,PFLOAT,PSTRING)
FCALLSCFUN2(INT,sep_get_number_grid_axes,SEP_GET_NUMBER_GRID_AXES,sep_get_number_grid_axes,STRING,PINT)

FCALLSCFUN5(INT,sep_put_val_by_index,SEP_PUT_VAL_BY_INDEX,sep_put_val_by_index,STRING,PINT,PINT,PINT,PVOID)
FCALLSCFUN5(INT,sep_get_val_by_index,SEP_GET_VAL_BY_INDEX,sep_get_val_by_index,STRING,PINT,PINT,PINT,PVOID)
FCALLSCFUN5(INT,sep_put_val_by_name,SEP_PUT_VAL_BY_NAME,sep_put_val_by_name,STRING,PINT,STRING,PINT,PVOID)
FCALLSCFUN5(INT,sep_get_val_by_name,SEP_GET_VAL_BY_NAME,sep_get_val_by_name,STRING,PINT,STRING,PINT,PVOID)
FCALLSCFUN5(INT,sep_extract_val_by_index,SEP_EXTRACT_VAL_BY_INDEX,sep_extract_val_by_index,STRING,PINT,PINT,PVOID,PVOID)
FCALLSCFUN5(INT,sep_insert_val_by_index,SEP_INSERT_VAL_BY_INDEX,sep_insert_val_by_index,STRING,PINT,PINT,PVOID,PVOID)
FCALLSCFUN4(INT,sep_put_val_headers,SEP_PUT_VAL_HEADERS,sep_put_val_headers,STRING,PINT,PINT,PVOID)
FCALLSCFUN4(INT,sep_get_val_headers,SEP_GET_VAL_HEADERS,sep_get_val_headers,STRING,PINT,PINT,PVOID)

FCALLSCFUN7(INT,sep_put_grid_window,SEP_PUT_GRID_WINDOW,sep_put_grid_window,STRING,PINT,PINT,PINT,PINT,PINT,PINT)
FCALLSCFUN7(INT,sep_get_grid_window,SEP_GET_GRID_WINDOW,sep_get_grid_window,STRING,PINT,PINT,PINT,PINT,PINT,PINT)
FCALLSCFUN2(INT,sep_copy_data_pointer,SEP_COPY_DATA_POINTER,sep_copy_data_pointer,STRING,STRING)
FCALLSCFUN2(INT,sep_copy_header_keys,SEP_COPY_HEADER_KEYS,sep_copy_header_keys,STRING,STRING)
FCALLSCFUN2(INT,sep_copy_grid,SEP_COPY_GRID,sep_copy_grid,STRING,STRING)
FCALLSCFUN2(INT,sep_copy_hff,SEP_COPY_HFF,sep_copy_hff,STRING,STRING)
FCALLSCFUN2(INT,sep_copy_gff,SEP_COPY_GFF,sep_copy_gff,STRING,STRING)
FCALLSCFUN1(INT,sep_set_no_headers,SEP_SET_NO_HEADERS,sep_set_no_headers,STRING)
FCALLSCFUN1(INT,sep_set_regular_grid,SEP_SET_REGULAR_GRID,sep_set_regular_grid,STRING)
FCALLSCFUN1(INT,sep_set_no_grid,SEP_SET_NO_GRID,sep_set_no_grid,STRING)
FCALLSCFUN6(INT,sep_reorder_data_fast,SEP_REORDER_DATA_FAST,sep_reorder_data_fast,STRING,STRING,INT,INT,PINT,INT)
FCALLSCFUN5(INT,sep_reorder_data,SEP_REORDER_DATA,sep_reorder_data,STRING,STRING,INT,INT,PINT)


/* stubs w/ an extra underscore for g77 compatibility */
FCALLSCFUN2(INT,sep_put_number_keys,SEP_PUT_NUMBER_KEYS_,sep_put_number_keys_,STRING,PINT)
FCALLSCFUN2(INT,sep_get_number_keys,SEP_GET_NUMBER_KEYS_,sep_get_number_keys_,STRING,PINT)

FCALLSCFUN5(INT,sep_put_key,SEP_PUT_KEY_,sep_put_key_,STRING,STRING,STRING,STRING,PINT)
FCALLSCFUN3(INT,sep_get_key_index,SEP_GET_KEY_INDEX_,sep_get_key_index_,STRING,STRING,PINT)
FCALLSCFUN3(INT,sep_get_key_name,SEP_GET_KEY_NAME_,sep_get_key_name_,STRING,PINT,PSTRING)
FCALLSCFUN3(INT,sep_get_key_type,SEP_GET_KEY_TYPE_,sep_get_key_type_,STRING,PINT,PSTRING)
FCALLSCFUN3(INT,sep_get_key_fmt,SEP_GET_KEY_FMT_,sep_get_key_fmt_,STRING,PINT,PSTRING)
FCALLSCFUN3(INT,sep_get_key_first_byte,SEP_GET_KEY_FIRST_BYTE_,sep_get_key_first_byte_,STRING,PINT,PINT)
FCALLSCFUN2(INT,sep_get_header_bytes,SEP_GET_HEADER_BYTES_,sep_get_header_bytes_,STRING,PINT)

FCALLSCFUN2(INT,sep_put_header_format_tag,SEP_PUT_HEADER_FORMAT_TAG_,sep_put_header_format_tag_,STRING,STRING)
FCALLSCFUN2(INT,fget_header_format_tag,SEP_GET_HEADER_FORMAT_TAG_,sep_get_header_format_tag_,STRING,PSTRING)

FCALLSCFUN2(INT,sep_put_grid_format_tag,SEP_PUT_GRID_FORMAT_TAG_,sep_put_grid_format_tag_,STRING,STRING)
FCALLSCFUN2(INT,sep_copy_grid,SEP_COPY_GRID_,sep_copy_grid_,STRING,STRING)
FCALLSCFUN2(INT,fget_grid_format_tag,SEP_GET_GRID_FORMAT_TAG_,sep_get_grid_format_tag_,STRING,PSTRING)

FCALLSCFUN6(INT,sep_put_data_axis_par,SEP_PUT_DATA_AXIS_PAR_,sep_put_data_axis_par_,STRING,PINT,PINT,PFLOAT,PFLOAT,STRING)
FCALLSCFUN6(INT,sep_get_data_axis_par,SEP_GET_DATA_AXIS_PAR_,sep_get_data_axis_par_,STRING,PINT,PINT,PFLOAT,PFLOAT,PSTRING)
FCALLSCFUN2(INT,sep_get_number_data_axes,SEP_GET_NUMBER_DATA_AXES_,sep_get_number_data_axes_,STRING,PINT)

FCALLSCFUN6(INT,sep_put_header_axis_par,SEP_PUT_HEADER_AXIS_PAR_,sep_put_header_axis_par_,STRING,PINT,PINT,PFLOAT,PFLOAT,STRING)
FCALLSCFUN6(INT,sep_get_header_axis_par,SEP_GET_HEADER_AXIS_PAR_,sep_get_header_axis_par_,STRING,PINT,PINT,PFLOAT,PFLOAT,PSTRING)
FCALLSCFUN2(INT,sep_get_number_header_axes,SEP_GET_NUMBER_HEADER_AXES_,sep_get_number_header_axes_,STRING,PINT)

FCALLSCFUN6(INT,sep_put_grid_axis_par,SEP_PUT_GRID_AXIS_PAR_,sep_put_grid_axis_par_,STRING,PINT,PINT,PFLOAT,PFLOAT,STRING)
FCALLSCFUN6(INT,sep_get_grid_axis_par,SEP_GET_GRID_AXIS_PAR_,sep_get_grid_axis_par_,STRING,PINT,PINT,PFLOAT,PFLOAT,PSTRING)
FCALLSCFUN2(INT,sep_get_number_grid_axes,SEP_GET_NUMBER_GRID_AXES_,sep_get_number_grid_axes_,STRING,PINT)

FCALLSCFUN5(INT,sep_put_val_by_index,SEP_PUT_VAL_BY_INDEX_,sep_put_val_by_index_,STRING,PINT,PINT,PINT,PVOID)
FCALLSCFUN5(INT,sep_get_val_by_index,SEP_GET_VAL_BY_INDEX_,sep_get_val_by_index_,STRING,PINT,PINT,PINT,PVOID)
FCALLSCFUN5(INT,sep_put_val_by_name,SEP_PUT_VAL_BY_NAME_,sep_put_val_by_name_,STRING,PINT,STRING,PINT,PVOID)
FCALLSCFUN5(INT,sep_get_val_by_name,SEP_GET_VAL_BY_NAME_,sep_get_val_by_name_,STRING,PINT,STRING,PINT,PVOID)
FCALLSCFUN5(INT,sep_extract_val_by_index,SEP_EXTRACT_VAL_BY_INDEX_,sep_extract_val_by_index_,STRING,PINT,PINT,PVOID,PVOID)
FCALLSCFUN5(INT,sep_insert_val_by_index,SEP_INSERT_VAL_BY_INDEX_,sep_insert_val_by_index_,STRING,PINT,PINT,PVOID,PVOID)

FCALLSCFUN4(INT,sep_put_val_headers,SEP_PUT_VAL_HEADERS_,sep_put_val_headers_,STRING,PINT,PINT,PVOID)
FCALLSCFUN4(INT,sep_get_val_headers,SEP_GET_VAL_HEADERS_,sep_get_val_headers_,STRING,PINT,PINT,PVOID)

FCALLSCFUN7(INT,sep_put_grid_window,SEP_PUT_GRID_WINDOW_,sep_put_grid_window_,STRING,PINT,PINT,PINT,PINT,PINT,PINT)
FCALLSCFUN7(INT,sep_get_grid_window,SEP_GET_GRID_WINDOW_,sep_get_grid_window_,STRING,PINT,PINT,PINT,PINT,PINT,PINT)
FCALLSCFUN2(INT,sep_copy_data_pointer,SEP_COPY_DATA_POINTER_,sep_copy_data_pointer_,STRING,STRING)
FCALLSCFUN2(INT,sep_copy_header_keys,SEP_COPY_HEADER_KEYS_,sep_copy_header_keys_,STRING,STRING)
FCALLSCFUN2(INT,sep_copy_hff,SEP_COPY_HFF_,sep_copy_hff_,STRING,STRING)
FCALLSCFUN2(INT,sep_copy_gff,SEP_COPY_GFF_,sep_copy_gff_,STRING,STRING)
FCALLSCFUN1(INT,sep_set_no_headers,SEP_SET_NO_HEADERS_,sep_set_no_headers_,STRING)
FCALLSCFUN1(INT,sep_set_regular_grid,SEP_SET_REGULAR_GRID_,sep_set_regular_grid_,STRING)
FCALLSCFUN1(INT,sep_set_no_grid,SEP_SET_NO_GRID_,sep_set_no_grid_,STRING)
FCALLSCFUN5(INT,sep_reorder_data,SEP_REORDER_DATA_,sep_reorder_data_,STRING,STRING,INT,INT,PINT)

FCALLSCFUN1(INT,sep_tag_is_pipe,SEP_TAG_IS_PIPE,sep_tag_is_pipe,STRING)
FCALLSCFUN1(INT,sep_tag_is_pipe,SEP_TAG_IS_PIPE_,sep_tag_is_pipe_,STRING)

FCALLSCSUB0(init_3d,INIT_3D,init_3d)
FCALLSCSUB0(init_3d,INIT_3D_,init_3d_)
FCALLSCSUB0(sep_3d_close,SEP_3D_CLOSE,sep_3d_close)
FCALLSCSUB0(sep_3d_close,SEP_3D_CLOSE_,sep_3d_close_)

