#ifndef STREAMIN3d_H
#define STREAMIN3d_H adfasd
#include <sep_file_types.h>

#include <stdio.h>
#include<prototypes.h>

#define MAX_EXTRA_KEYS 15
#ifndef SU_NKEYS
#define SU_NKEYS 80
#endif
#define SEP_3D_STRING_LEN 128 /* all string entries in _sep_3d struct are
                                 Fortran CHARACTER*128 */


/* types of stream , the way a file is opened */
enum usage_type { INPUT, OUTPUT, SCRATCH, INQUIRE };
enum data_type { FLOAT, INTEGER, COMPLEX, BYTE, UNKNOWN};
enum file_type { REGULAR, HEADER, GRID,UNSPECIFIED};

typedef long long sep_coord_type;








struct _sep_3d;  /* declare it as a type that can be pointed at */


struct _sep_3d{
struct _sep_3d  *next; /*next sep3d structure */
struct _sep_3d  *prev; /*previous sep3d structure*/
char *name    ;  /*structure name */
enum data_type data_format;/*data type (float, int , or complex) */
enum file_type file_format; /*regular, header, or grid*/
enum usage_type usage;    /*data type (in, out, scratch) */

int     ndims;          /*number of dimensions in dataset */
int      *n;         /*size of dataset */
float   *o,*d;     /*other axes descriptors for regular portion*/
char   **label,**unit;  /*other axes descriptors for regular portion*/
int *nwind,*fwind,*jwind;  /*portion of the headers in memory*/

int     nkeys;        /*number of keys in dataset*/
char   **keyname,**keytype,**keyfmt ; /*key name, type, and format */

int    nh              ;/*number of headers being currently stored */
int  *headers                  ; /*headers buffer */
int     *drn;        /*data record number for the dataset*/
long long int ntraces_wrote;   /* number of traces written*/

/*
int    ng              ;  size of the grid currently being stored 
int  *grid             ;  grid values 
*/

int     in_order;        /*in_order whether or not tag is synched*/
int     ncoord;         /*coordinate array size (usually nh)*/
sep_coord_type     *coord;      /*coordinate  (in helix space) of trace */

long long  int  ntraces;        /*number of traces in dataset*/
int     nextra_keys;    /* number of keys before on the fly math */
int     nkeys_in;    /* number of keys before on the fly math */
char   **exp ;   /*buffer for header math on the fly */
int  wrote_data,wrote_headers; /*wheter or not we have written data/headers*/
int   su_input; /* Data is SU and an input file (could be reshaped ) */
 
#ifdef SU_SUPPORT
/*now  for the su io stuff */
int  key_index[SU_NKEYS+MAX_EXTRA_KEYS]; /*index relate su keys to seplib keys */
float  float_fract[SU_NKEYS+MAX_EXTRA_KEYS];/*conversion between su->seplib  (0->int)*/

int nmem ; /* The number of traces that we will buffer in su format */
char *su_tr_block,*su_hdr_block; /*blocks to hold buffered header and traces*/

int extra_keys;  /*whether there are some extra (non-su) keys in the dataset */
int beg_trace, end_trace;/*the begining and ending trace currently stored */
int last_trace; /* last_trace read or written */
int nextra_mapping; /*extra or different mapping exists */
char extra_name[MAX_EXTRA_KEYS][257]; /*offset of extra key */
int extra_offset[MAX_EXTRA_KEYS]; /*offset of extra key */
char extra_type[MAX_EXTRA_KEYS][5]; /*offset of extra key */
int      su_ndims;         /*size of dataset in su */
int      *n_su;         /*size of dataset in su */
float   *o_su,*d_su;     /*sampling of dataset in su*/
#endif
};
typedef struct _sep_3d sep_3d;


#endif
/*  $Id: superset_struct.h,v 1.2 2004/04/08 22:32:27 bob Exp $ */
