#ifndef SEPAUX_H
#define SEPAUX_H
/*this should be changed - don't need to pull in all of the includes that
this implies */
#include<prototypes.h> 

#ifdef __cplusplus
extern "C" {
#endif

#define MAXHEAD 100 /* Maximum number of headers */

/* convenient structures for storing key header information */

struct _sep123
{
    int             n;
    float           o;
    float           d;
};

struct _sep123def
{
    int             n;
    int             o;
    int             d;
};

struct _sep
{
    int             esize;
    int       transp;
    char           *header;
    char            tag[100];
    char            in[100];
    struct _sep123  standard[4];
    struct _sep123def set[4];
};



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern int arrayreadf(char*, char*, int, float*,float);
extern int arrayreadi(char*, char*, int, int*,int);
extern int compare_c(char*, char*);
extern int parcat(char*,int,char*);
extern int sgainpar (char*,float*,int*,int*,int*,int*,float*,float*,float*,float*,float*,float*,float*,float*,int*,int*);
extern float cent (float,float*,int);
extern void pqueue_init (int);
extern void pqueue_close (void);
extern void pqueue_insert (float*);
extern float* pqueue_extract (void);
extern int time_window(int,float*,float*,char*);
extern int init_multihdr(int,char**,struct _sep**);
extern void struct_dump(struct _sep*);
extern void process_one(struct _sep *);
extern int  pad_it(char *,char *,int,int*,int*,int*,int,int);
extern int pad_portion(char*,char*,int,int,int,int*,int*,int*,int,int,int*,int*,int grab_input(char *,int*,int*,int));
extern int oc_patch_put(char *tag, int ipatch, int esize,char *buf);
extern int oc_patch_get(char *tag, int ipatch, int esize,char *buf);
extern void oc_patch_init(int ndim,int *nin, int *nwind, int *npatch, int *npatches);
extern void oc_patch_clean();
extern void oc_patch_bound(int ipatch,int esize, float *buf);
extern int init_loop_calc(int ndims,int *n,char *string,int max_size);
extern int do_sep_loop(char *name, int *nwind, int *fwind);
_XFUNCPROTOEND
#else
extern int arrayreadf();
extern int arrayreadi();
extern int compare_c();
extern int parcat();
extern int sgainpar ();
extern float cent();
extern void pqueue_init ();
extern void pqueue_close ();
extern void pqueue_insert ();
extern float* pqueue_extract ();
extern int time_window();
extern int init_multihdr();
extern void struct_dump();
extern void process_one();
extern int  pad_it();
extern int pad_portion();
extern int init_loop_calc();
extern int do_sep_loop();
#endif
#ifdef __cplusplus
}
#endif
#endif
