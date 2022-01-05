#ifndef STREAMINF_H
#define STREAMINF_H adfasd

#include <stdio.h>
#include<prototypes.h>

#include "../include/fastpar.h"
#include "strformats.h"


/*#define SEPPOINTNULL ((void *) 0)*/
/*#define SEPSTRNULL ((char *) 0)*/
#define SEPPOINTNULL 0
#define SEPSTRNULL 0
/* maximum number of files in a seplib dataset */
#define MAX_SEP_FILE 1000

/* default size of multi part file sections (In Mbytes)  */
#define SEP_DIR_SIZE  9999000
#ifdef FILE_2GB
#define SEP_FILE_SIZE 1999
#else
#define SEP_FILE_SIZE 999999
#endif
#ifndef SEP_BUFSIZ
#define SEP_BUFSIZ 4096
#endif

#define EOT 004

/* types of stream , the way a file is opened */
enum strtype { STREAMIN, STREAMOUT, STREAMINOUT, STREAMSOCKOUT, STREAMSCR};  

/* types of I/O , the underlying mechanism */
enum IOtyp { FD_IO, FILE_IO, CM_SDA_IO, CM_REG_IO, MULTI_FD_IO };

/* Reqests types for turning a tag into an info structure */
enum tagtype { TAG_OUT, TAG_IN, TAG_INOUT, TAG_SOCKET, TAG_INQUIRE, TAG_SCR}; 

struct _streamh;  /* declare it as a type that can be pointed at */

#include "io_func_defs.h"

struct _streamh {
	struct _streamh  *next;
	struct _streamh  *prev;
	char *tagname;
	enum strtype entrytype;      /* STREAMIN, STREAMOUT, STREAMSCR
                                        STREAMINOUT, STREAMSOCKOUT */

        int valid;   /* is stream valid ( dataset OK etc. ) */

	char *headername;    /* name of header (file, pipe or socket ) */
	FILE *headfile;      /* file pointer associated with the header */

	char *headerformatfile; /* name of associated header that describes 
                                  the properties of the Data Records       */
	char *gridformatfile; /* name of associated Grdi Format File */


	/* for input datasets only */
	char *headerbuf;     /* buffer containing a copy of input header */
	int hdrlen;          /* length of the input header */
	hash_item **hetch_queue; /* hashed queue of par in this header */
	int hqlen;           /* length of queue */
	hash_item **tetch_queue; /* hashed tetch queue for this header */
	int tqlen;           /* length of queue */

        int backseekable; /* can we backward seek on this input stream */

	/* for output datasets only */
	int header_title;  /* has the title line been written in the header */
	int header_outname;/* has the output dataset name been written ? */
	int ready_out;     /* Is the tag ready for data output ? */

	/* for all datasets */

	char *dataname;      /* name of data file (filename or stream) */

	enum IOtyp iotype;   /* type of IO for this dataset */

	void *ioinf;         /* pointer to something that the different
				I/O routines can use, the actual thing
				pointed at is IOtyp dependent */

	/* pointers to the functions to perform I/O on this dataset */
	PF_OPEN open_func;
	PF_CLOSE close_func;
	PF_READ read_func;
	PF_WRITE write_func;
	PF_SEEK seek_func;
	PF_SIZE size_func;

	/* bogus entries for backwards compatibility */
	FILE* streamfile;
	int streamfd;

	/* format of the data, see strformats.h for what they mean */
	int format_num; 
	int isapipe;  /* is this a pipe, socket or other IPC mechanism */
	int sockfd;   /* aux socket for IPC with the other end of stream */

	/* for Sep 90 data sets */
	int n_key; /*Number of keys in headers*/
	int *key_bytes; /*pointers to number of bytes in each key*/
};
typedef struct _streamh streaminf;

struct _sep_self_doc {
  char *line;
	struct _sep_self_doc  *next;
};
typedef struct _sep_self_doc sep_self_doc;

#endif /*END STREAMINF_H*/
