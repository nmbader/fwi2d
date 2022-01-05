#include<seplib.h>
#ifndef MPI_SEP_LIB_H
#define MPI_SEP_LIB_H 1

#ifdef __cplusplus
extern "C" {
#endif

#if NeedFunctionPrototypes
int mpi_sep_send_args(int,int,int);

int sep_mpi_tag_sum(int ithread,char *tag_from, int ithread_to, char *tag_to, int bs, int nb,int add);
int sep_mpi_tag_sum_index(int ithread,char *tag_from, int ithread_to, char *tag_to, int bs, int nb,int add,int index);
int sep_mpi_tag_bcast_index(int ithread_from,char *tag_from, int ithread, char *tag_to, int bs, int nb,int index);
int sep_mpi_tag_bcast(int ithread_from,char *tag_from, int ithread, char *tag_to, int bs, int nb);

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
int mpi_sep_open_lock(char *fname);
int mpi_sep_close_lock(char *fname,int fd);
int mpi_sep_check_lock_file(char *fname);
#endif

int sep_mpi_thread_portion(char*,int);
int sep_mpi_update_block(char *tag,int nblock,int input);
int sep_mpi_begin_io(char *tag,int nblock);
int sep_mpi_tag_distribute(int ithread_from,char *tag_from, int ithread, char *tag_to, int bs, int nb, int *send_to);

int sep_mpi_tag_combine_index(int ithread,char *tag_from, int ithread_to, char *tag_to, int bs, int nb, int *send_from,int add,int index);
int sep_mpi_tag_combine_big(int ithread,char *tag_from, int ithread_to, char *tag_to, int bs, int nb, int nb2, int *send_from,int add);
int sep_mpi_tag_combine(int ithread,char *tag_from, int ithread_to, char *tag_to, int bs, int nb, int *send_from,int add);

int mpi_sep_pass_args();
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
int sep_open_lock(char *tag);
int sep_close_lock(char *tag);
int mpi_sep_open_lock(char *fname);
int mpi_sep_close_lock(char *fname,int fd);
int mpi_sep_check_lock_file(char *fname);
int mpi_sep_create_group(int ngrp, int *use);
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
#endif

#ifdef __cplusplus
}
#endif

#endif

