#ifndef FASTPAR_H
#define FASTPAR_H
#include <prototypes.h>
/* added timestamp 3-8-86  stew */
typedef struct hash_dummy {
	struct hash_dummy *next;
	char *tag; int tlen;
	char *val; int vlen;
	int timestamp;
	} hash_item;

typedef union { int *i; float *f; double *g ; char *s; long long *m;} MIXED;



#if NeedFunctionPrototypes/* { */
_XFUNCPROTOBEGIN 
extern hash_item** new_queue( int );
extern int getpar_decode(hash_item **,int ,const char*,char*,MIXED);
extern void putch_format(char*,char*,char*,MIXED);
extern hash_item *getpar_hash_lookup(register hash_item **,register int,
register char*,register int);
extern void getpar_push_input(register char*,register int);
extern void getpar_scan(register hash_item **,register int);
extern void getpar_string_store(hash_item **,int,register char *);
extern void getpar_hash_store(hash_item**,int,register char*,register char*,
     register int,int);
extern void getpar_stack_par( char*);
_XFUNCPROTOEND 
#else/* } !NeedFunctionPrototypes { */
extern hash_item** new_queue();
extern int getpar_decode();
extern void putch_format();
extern hash_item *getpar_hash_lookup();
extern void getpar_push_input();
extern void getpar_scan();
extern void getpar_string_store();
extern void getpar_hash_store();
extern void getpar_stack_par();
#endif/* } */
#endif/* FASTPAR_H */
/*  $Id: fastpar.h,v 1.1.1.1 2004/03/25 06:37:24 cvs Exp $ */
