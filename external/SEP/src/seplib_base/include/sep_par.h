#include<sep3dc.h>
#ifndef MPI_SEP_LIB_H
#define MPI_SEP_LIB_H 1

#ifdef __cplusplus
extern "C" {
#endif

struct _file_locks;
                                                                                
struct _file_locks{
struct _file_locks  *next; /*next sep3d structure */
char *tag;
char *lock_file;
int  fd;
};

typedef struct _file_locks sep_file_lock;
struct _distrib{
  int nsect;  /*number of sections we are dealing with*/
  char **tag_sect; /*the section names we are dealing with*/
  int *isect;  /*what section corresponds to what  portion of the divided axis/axes*/
  int *sect_thread; /*the thread that each section belongs to*/
  int *dff_axis;  /*the axis/exes we are spread along*/
  int **axis_beg; /*the begining and ending portion for each section we are dealing with [nsect][naxis]*/
  int **axis_end; /*the ending of the above*/
  int nown; /*how many sections we own*/
  int *nblock; /*the number of blocks along each split axis*/
  int *noverlap; /*the amount of overlap in these blocks*/
  int *iown; /*the list of sections i own*/
  int *ilocal_pt; /*a pointer of size nsect, -1 if i don't the own section otherwise the local section number*/
  int ndim; /*the number of axes we are split along*/
};
typedef struct _distrib distrib;

int read_local_datasets( sep3d *sects, sep3d *data, distrib *spread, int *nwind, int *fwind, float *buf1, float *buf2,float *buf3);
int get_distrib_info( distrib *spread, int output,char *prefix);
int calc_output_sections(sep3d *big, distrib *spread);
int calc_input_sections(sep3d *big, distrib *spread,sep3d *structs);

int combo_data(int impi, int nmpi,float *data,float *buf, int n123);
int create_bounds(int isect, distrib *spread, int ndims,int *ngrid,int *nwind, int *fwind, int *b, int *e, int *s);
int patch_boundary(sep3d *data,distrib *spread,int *nwind,  int isect, int *nout, int *fout, int *fsect, float *buf2);
int create_global_file(distrib *spread,sep3d *sects, sep3d *combo);
void create_3d(sep3d *data, distrib *spread,int *nsect,int idim,  int isect,
int *nuse, int *nblock, int *ibeg0,int *ibeg1, int *iend0, int *iend1);
int local_window(int isect, sep3d *data, distrib *spread, int *nwind,
  int *fwind, int *nsect, int *fsect, int  *nout, int *fout);





#if NeedFunctionPrototypes
int mpi_sep_send_args(int,int,int);

int sep_mpi_tag_sum(int ithread,char *tag_from, int ithread_to, char *tag_to, int bs, int nb,int add);
int sep_mpi_tag_sum_index(int ithread,char *tag_from, int ithread_to, char *tag_to, int bs, int nb,int add,int index);
int sep_mpi_tag_bcast_index(int ithread_from,char *tag_from, int ithread, char *tag_to, int bs, int nb,int index);
int sep_mpi_tag_bcast(int ithread_from,char *tag_from, int ithread, char *tag_to, int bs, int nb);
void sep_close_lock(char *tag);
void sep_open_lock(char *tag);

#ifdef SEP_MPI
#include<mpi.h>
MPI_Comm *sep_mpi_return_com(int index);
int sep_mpi_init();
int sep_mpi_auto_io_init();
int sep_mpi_io_tags(char *tagin, char *tagout);
int sep_mpi_proc_begin(char *tagin, char *tagout);
int sep_mpi_do_proc();
int sep_mpi_inttag(char *tag);
int sep_mpi_outtag(char *tag);
int sep_mpi_in_out_tags(char *tagin, char *tag_out);
int sep_mpi_finish();
int sep_mpi_distrib_in_tag(char *tag_in, char  *io_tag, int iaxis, int nblock,char *distrib_type);
int sep_mpi_distrib_out_tag(char *tag_in,  char  *io_tag,int iaxis, int nblock,char *distrib_type);
int  sep_mpi_collect_all(char *tag);
int  sep_mpi_distrib_all(char *tag);

#endif

int sep_mpi_thread_portion(char*,int);
int sep_mpi_update_block(char *tag,int nblock,int input);
int sep_mpi_begin_io(char *tag,int nblock);
int sep_mpi_tag_distribute(int ithread_from,char *tag_from, int ithread, char *tag_to, int bs, int nb, int *send_to);

int sep_mpi_tag_combine_index(int ithread,char *tag_from, int ithread_to, char *tag_to, int bs, int nb, int *send_from,int add,int index);
int sep_mpi_tag_combine_big(int ithread,char *tag_from, int ithread_to, char *tag_to, int bs, int nb, int nb2, int *send_from,int add);
int sep_mpi_tag_combine(int ithread,char *tag_from, int ithread_to, char *tag_to, int bs, int nb, int *send_from,int add);

int mpi_sep_receive_args();
int mpi_sep_set_dead(int nproc, int *dead_list);
int mpi_sep_loop_init(int,int,int);
int mpi_sep_loop_close();
int mpi_sep_next_proc(int,int);
int mpi_sep_get_proc_num(int);
int mpi_sep_clean();
int mpi_sep_get_loop_num(int);
int mpi_sep_finish_iproc(int);
int mpi_sep_finish_iloop(int);
int mpi_sep_check_dead();
int mpi_sep_print_status();
int mpi_sep_next_group(int, int*);
int mpi_sep_valid_nproc();
int mpi_sep_valid_processors(int*);
int mpi_sep_jobs_running();
int mpi_sep_open_lock(char *fname);
int mpi_sep_close_lock(char *fname,int fd);
int mpi_sep_check_lock_file(char *fname);
int mpi_sep_create_group(int ngrp, int *use);
int sep_num_thread(void);
int sep_thread_num(void);
#else
void mpi_open_lock();
void mpi_close_lock();
int mpi_sep_set_dead();
int mpi_sep_send_args();
int mpi_sep_receive_args();
int mpi_sep_loop_init();
int mpi_sep_loop_close();
int mpi_sep_next_proc();
int mpi_sep_get_loop_num();
int mpi_sep_get_proc_num();
int mpi_sep_finish_iproc();
int mpi_sep_finish_iloop();
int mpi_sep_check_dead();
int mpi_sep_print_status();
int mpi_sep_next_group();
int mpi_sep_valid_nproc();
int mpi_sep_valid_processors();
int mpi_sep_jobs_running();
int sep_num_thread();
int sep_thread_num();
#endif

#ifdef __cplusplus
}
#endif

#endif

