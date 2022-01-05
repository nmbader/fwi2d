#include "./../include/sepstream.h"

/* utility functions used in handling the list of sepstream functions */


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern void init_fd_io(streaminf *);
extern char* expand_info( char* , streaminf* );
extern void init_file_io(streaminf*);
extern void init_io(streaminf*);
extern void init_multifd_io(streaminf*);
extern void open_instream(streaminf*);
extern void sepstr_in_head( streaminf*);
extern void readhdr(streaminf*);
extern void open_inoutstream( streaminf*);
extern void sepstr_inout_head( streaminf*);
extern void open_outstream( streaminf *);
extern void sepstr_out_head( streaminf*);
extern void open_scrstream( streaminf*);
extern void sepstr_scr_head(streaminf*);
extern void open_socketstream( streaminf*);
extern void sepstr_socket_head(streaminf*);
extern void syncin(streaminf*);
extern void syncout(streaminf*);
extern streaminf* sepstr_head(void);
extern void sepstr_del(streaminf*);
extern void sepstr_addstart( streaminf*);
extern void sepstr_addend( streaminf * );
extern void print_streaminf( streaminf* );
extern void sepstr_copyh(streaminf*,streaminf*);
extern void outname(streaminf*);
extern streaminf *tag_info(const char*, enum tagtype);
extern streaminf *fd_info( int);
extern sep_file_size_t fsize( int );
extern streaminf *sepstr_new(char*,enum tagtype);
extern void sepstr_hclose(streaminf*);
extern void sepstr_ready_out(streaminf*);
extern int sepstrpar(streaminf*,char*,char*,void*);
extern int sepstrput(streaminf*,char*,char*,void*);
extern int sepstrputlast(streaminf*, char*, char*, void*);
extern void write_title(streaminf*);
_XFUNCPROTOEND
#else
extern streaminf * sepstr_head();
extern void init_fd_io();
extern char* expand_info();
extern void init_file_io();
extern void init_io();
extern void init_multifd_io();
extern void open_instream();
extern void sepstr_in_head();
extern void readhdr();
extern void open_inoutstream();
extern void sepstr_inout_head();
extern void open_outstream();
extern void sepstr_out_head();
extern void sepstr_copyh();
extern void outname();
extern void open_scrstream();
extern void sepstr_scr_head();
extern void open_socketstream();
extern void sepstr_socket_head();
extern void syncin();
extern void syncout();
extern void sepstr_del();
extern void sepstr_addstart();
extern void sepstr_addend();
extern void print_streaminf();
extern streaminf * tag_info();
extern streaminf * fd_info();
extern sep_file_size_t fsize();
extern streaminf *sepstr_new();
extern void sepstr_hclose();
extern void sepstr_ready_out();
extern int sepstrpar();
extern int sepstrput();
extern int sepstrputlast();
extern void write_title();
#endif


#ifdef BOOGSA
extern void write_outname(streaminf*);
extern void sepstr_del(streaminf*);
extern int sepstrput( streaminf *, char *, char *, void* );
extern int sepstrpar( streaminf *, char*, char *, void * );
extern void readhdr(  streaminf *);
extern void print_streaminf( streaminf* );
extern int fdordinary(int);
extern int isordinary(char*);
extern int sepstrputlast( streaminf*,char*,char*,void*);
extern void sepstr_in_head();
extern void sepstr_out_head();
extern void sepstr_inout_head();
extern void write_title();
extern void write_outname();
extern void sepstr_ready_out();
extern void sepstr_hclose();
extern void sepstr_del();
extern void sepstr_copyh();
extern int sepstrput( );
extern int sepstrpar( );
extern void readhdr();
extern void syncin();
extern void syncout();
extern streaminf * tag_info();
extern streaminf * fd_info();
extern void print_streaminf();
extern char* expand_info();
extern sep_file_size_t fsize();
extern void open_outstream();
extern void init_fd_io();
extern int fdordinary();
extern int isordinary();
extern void init_io();
extern void init_fd_io();
extern void init_file_io();
extern void init_multifd_io();
extern void open_instream();
extern int sepstrputlast();
#endif /* END ANSI CHECK */
/*  $Id: streamlist.h,v 1.1.1.1 2004/03/25 06:37:24 cvs Exp $ */
