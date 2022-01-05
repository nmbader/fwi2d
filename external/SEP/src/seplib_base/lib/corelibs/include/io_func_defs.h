#ifndef IO_FUNC_DEFS_H
#define IO_FUNC_DEFS_H 1

#include <sys/types.h>
#include "../include/sep_file_types.h"

/* typedef for the offset type goes here */

/* six functions define how to access a dataset */

#if defined(__STDC__) || defined(__stdc__)
typedef void (*PF_OPEN)( struct _streamh* , void** );
typedef void (*PF_CLOSE)( struct _streamh*, void* );
typedef ssize_t (*PF_READ)( struct _streamh*, void*, void*, size_t  );
typedef ssize_t (*PF_WRITE)( struct _streamh*, void*, void*, size_t  );
typedef sep_file_size_t (*PF_SEEK)( struct _streamh*, void*, sep_off_t, int );
typedef sep_file_size_t (*PF_SIZE)( struct _streamh*, void*);
#else
typedef void (*PF_OPEN)();
typedef void (*PF_CLOSE)();
typedef ssize_t (*PF_READ)();
typedef ssize_t (*PF_WRITE)();
typedef sep_file_size_t (*PF_SEEK)();
typedef sep_file_size_t (*PF_SIZE)();
#endif

#endif




