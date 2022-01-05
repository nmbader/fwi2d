#ifndef VPLOT_H
#define VPLOT_H
#include <stdio.h>


#ifndef NO_VPLOT_DEFINES
/*
 * Weird backwards-compatible units
 */
#define RPERIN 		600.0f	/* vplot units per inch */
#define HATCHPERIN	100.0f	/* Hatch units per inch */
#define TXPERIN 	33.0f	/* Text units per inch */
#define FATPERIN	200.0f	/* Fatness units per inch */
/*
 * Height in inches of "standard" device, standard style
 */
#define STANDARD_HEIGHT 10.24f
/*
 * Height in inches of "standard" device, rotated style
 */
#define ROTATED_HEIGHT	 7.5f
/*
 * Aspect ratio of the default window (height/width)
 */
#define SCREEN_RATIO 0.75f
#define VP_MAX 54.6f		/* absolute maximum x or y coordinate in inches */

/*
 * text alignment enumerations
 */
/* horizontal */
#define	TH_NORMAL	0
#define	TH_LEFT		1
#define	TH_CENTER	2
#define	TH_RIGHT	3
#define TH_SYMBOL	4

/* vertical */
#define	TV_NORMAL	0
#define TV_BOTTOM	1
#define TV_BASE		2
#define TV_HALF		3
#define TV_CAP		4
#define TV_TOP		5
#define TV_SYMBOL	6

/*
 * text precision enumerations
 */
#define STRING	0
#define CHAR	1
#define STROKE	2
/* leave it what it already was */
#define NO_CHANGE -1

/*
 * text overlay enumerations
 */
#define OVLY_NORMAL	0
#define OVLY_BOX	1
#define OVLY_SHADE	2
#define OVLY_SHADE_BOX	3

/*
 * colors
 */
#define BLACK    0
#define BLUE     1
#define RED      2
#define PURPLE   3
#define GREEN    4
#define CYAN     5
#define YELLOW   6
#define WHITE    7

/*
 * Coordinate Origin
 */
#define STANDARD	0
#define ROTATED		1
#define ABSOLUTE 	3

/*
 * Fonts
 */

#define PEN		0
#define ROMANS		1
#define ROMAND		2
#define ROMANC		3
#define ROMANT		4
#define ITALICC		5
#define ITALICT		6
#define SCRIPTS		7
#define SCRIPTC		8
#define GREEKS		9
#define GREEKC		10
#define CYRILC		11
#define GOTHGBT		12
#define GOTHGRT		13
#define GOTHITT		14
#define MATH		15
#define MISC		16

/*
 * vplot metafile op-codes
 */

#define VP_SETSTYLE		'S'

#define VP_MOVE			'm'
#define VP_DRAW			'd'
#define VP_PLINE	    	'L'
#define VP_PMARK	   	'M'
#define VP_TEXT			'T'
#define VP_GTEXT		'G'
#define VP_AREA			'A'
#define VP_OLDAREA		'a'
#define VP_BYTE_RASTER		'R'
#define VP_BIT_RASTER		'r'
#define VP_SHORT_RASTER		'B'

#define VP_ERASE		'e'
#define VP_BREAK		'b'
#define VP_PURGE		'p'
#define VP_NOOP			'n'

#define VP_ORIGIN		'o'
#define VP_WINDOW		'w'
#define VP_UORIGSCL             'u'

#define VP_FAT			'f'
#define VP_SETDASH		's'
#define VP_COLOR		'c'
#define VP_SET_COLOR_TABLE	'C'
#define VP_TXALIGN		'J'
#define VP_TXFONTPREC		'F'
#define VP_PATLOAD		'l'
#define VP_OVERLAY		'v'

#define VP_MESSAGE		'z'
#define VP_BEGIN_GROUP		'['
#define VP_END_GROUP		']'

/* Hopefully now dead primitives */
#define VP_OLDTEXT		't'


#endif/* NO_VPLOT_DEFINES */

#ifdef __cplusplus
extern "C" {
#endif

struct txalign {
	int hor;
	int ver;
};
extern short geth(register FILE*);
extern float getf(register FILE *);
extern int name_to_coltab(char *colname, int nocol, float *red, float *green, float *blue);
extern int puth (register int w, register FILE *iop);
extern int putf (float w, FILE *iop);
extern int vp_arc(float x, float y, float r, float angstart, float angend);
extern int vp_arc_g(double x, double y, double r, double angstart, double angend);
extern int vp_uarc(float x, float y, float r, float angstart, float angend);
extern int vp_uarc_g(double x, double y, double r, double angstart, double angend);
extern int vp_area(float *xp, float *yp,int  lp, int fat,int  xmask,int  ymask);
extern int vp_area_g(double *xp, double *yp,int  lp, int fat,int  xmask,int  ymask);
extern int vp_uarea(float *xp, float *yp,int  lp, int fat,int  xmask,int  ymask);
extern int vp_uarea_g(double *xp, double *yp,int  lp, int fat,int  xmask,int  ymask);
extern int vp_arrow (float x0,float  y0,float  x,float  y,float r);
extern int vp_arrow_g (double x0,double  y0,double  x,double  y,double r);
extern int vp_uarrow (float x0,float  y0,float  x,float  y,float r);
extern int vp_uarrow_g (double x0,double  y0,double  x,double  y,double r);
extern int vp_bgroup(char *string);
extern int vp_bgroup_g(char *string);
extern int vp_break(void);
extern int vp_break_g(void);
extern int vp_circle(float, float, float);
extern int vp_circle_g(double, double, double);
extern int vp_ucircle(float, float, float);
extern int vp_ucircle_g(double, double, double);
extern int vp_clip(float, float, float, float);
extern int vp_clip_g(double, double, double, double);
extern int vp_uclip(float, float, float, float);
extern int vp_uclip_g(double, double, double, double);
extern int vp_color(int col);
extern int vp_color_g(int col);
extern int vp_coltab(int,float,float,float);
extern int vp_coltab_g(int,float,float,float);
extern int vp_dash(float,float,float,float);
extern int vp_dash_g(double,double,double,double);
extern int vp_draw(float,float);
extern int vp_draw_g(double,double);
extern int vp_udraw(float,float);
extern int vp_udraw_g(double,double);
extern int vp_egroup(void);
extern int vp_egroup_g(void);
extern int vp_endplot(void);
extern int vp_endplot_g(void);
extern int vp_erase(void);
extern int vp_erase_g(void);
extern int vp_fat(int FATNESS);
extern int vp_fat_g(int FATNESS);
extern int vp_file(char *filename);
extern int vp_file_g(char *filename);
extern int vp_filep(FILE *filepntr);
extern int vp_filep_g(FILE *filepntr);
extern int vp_fill(float*,float*,int);
extern int vp_fill_g(double*,double*,int);
extern int vp_ufill(float*,float*,int);
extern int vp_ufill_g(double*,double*,int);
extern int vp_gtext(float X,float Y,float  XPATH,float  YPATH,float  XUP,float  YUP,char *string);
extern int vp_gtext_g(double X,double Y,double  XPATH,double  YPATH,double  XUP,double  YUP,char *string);
extern int vp_hatchload(int ANGLE,int  NUMHATCH, int IHATCH, int *hatcharray);
extern int vp_hatchload_g(int ANGLE,int  NUMHATCH, int IHATCH, int *hatcharray);
extern int vp_message(char*); 
extern int vp_message_g(char*); 
extern int vp_move(float,float);
extern int vp_move_g(double,double);
extern int vp_umove(float,float);
extern int vp_umove_g(double,double);
extern int vp_orig(float,float);
extern int vp_orig_g(double,double);
extern int vp_uorig(float,float);
extern int vp_uorig_g(double,double);
extern int vp_patload(int PPI,int  NX,int  NY,int  IPAT,int  *colarray);
extern int vp_patload_g(int PPI,int  NX,int  NY,int  IPAT,int  *colarray);
extern int vp_pendn(float,float); 
extern int vp_pendn_g(double,double); 
extern int vp_upendn(float,float); 
extern int vp_upendn_g(double,double); 
extern int vp_penup(void);
extern int vp_penup_g(void);
extern int vp_pline(float *XP, float *yp,int LP);
extern int vp_pline_g(double *XP, double *yp,int LP);
extern int vp_upline(float *XP, float *yp,int LP);
extern int vp_upline_g(double *XP, double *yp,int LP);
extern int vp_plot_init(void);
extern int vp_plot_init_g(void);
extern int vp_plot(float,float,int);
extern int vp_plot_g(double,double,int);
extern int vp_pmark(int NPTS,int  MTYPE, int MSIZE, float *xp,float *yp);
extern int vp_pmark_g(int NPTS,int  MTYPE, int MSIZE, double *xp,double *yp);
extern int vp_upmark(int NPTS,int  MTYPE, int MSIZE, float *xp,float *yp);
extern int vp_upmark_g(int NPTS,int  MTYPE, int MSIZE, double *xp,double *yp);
extern int vp_purge(void);
extern int vp_purge_g(void);
extern int vp_rascol16tab(int nreserve,char *colname);
extern int vp_rascol16tab_g(int nreserve,char *colname);
extern int vp_rascoltab(int nreserve,char *colname);
extern int vp_rascoltab_g(int nreserve,char *colname);
extern int vp_rastershort (unsigned short *array, int BLAST, int BIT, int OFFSET, int XPIX, int YPIX, float XLL, float YLL,float  PPI,float *xur,float *yur, int ORIENT, int INVERT);
extern int vp_rastershort_g (unsigned short *array, int BLAST, int BIT, int OFFSET, int XPIX, int YPIX, double XLL, double YLL,double  PPI,double *xur,double *yur, int ORIENT, int INVERT);
extern int vp_raster (unsigned char *array, int BLAST, int BIT, int OFFSET, int XPIX, int YPIX, float XLL, float YLL,float  PPI,float *xur,float *yur, int ORIENT, int INVERT);
extern int vp_raster_g (unsigned char *array, int BLAST, int BIT, int OFFSET, int XPIX, int YPIX, double XLL, double YLL,double  PPI,double *xur,double *yur, int ORIENT, int INVERT);
extern int vp_scale(float,float);
extern int vp_scale_g(double,double);
extern int vp_setdash(float *dashp,float *gapp, int LP);
extern int vp_setdash_g(double *dashp,double *gapp, int LP);
extern int vp_stretch(float XMIN,float YMIN,float  XMAX,float YMAX);
extern int vp_stretch_g(double XMIN,double YMIN,double  XMAX,double YMAX);
extern int vp_style(int);
extern int vp_style_g(int);
extern int vp_text(float X,float  Y,int  SIZE,int  ORIENT, char *string);
extern int vp_text_g(double X,double Y,int  SIZE,int  ORIENT, char *string);
extern int vp_tfont(int,int,int);
extern int vp_tfont_g(int,int,int);
extern int vp_tjust(int,int);
extern int vp_tjust_g(int,int);
extern int vp_ucircle(float,float,float);
extern int vp_ucircle_g(double,double,double);
extern int vp_uclip(float,float,float,float);
extern int vp_uclip_g(double,double,double,double);
extern int vp_gtext(float X,float Y, float XPATH,float YPATH,float XUP,float YUP, char *string);
extern int vp_gtext_g(double X,double Y, double XPATH,double YPATH,double XUP,double YUP, char *string);
extern int vp_ugtext(float X,float Y, float XPATH,float YPATH,float XUP,float YUP, char *string);
extern int vp_ugtext_g(double X,double Y, double XPATH,double YPATH,double XUP,double YUP, char *string);
extern int vp_uplot(float X, float Y,int  DOWN);
extern int vp_uplot_g(double X, double Y,int  DOWN);
extern int vp_uraster (unsigned char *array, int BLAST,int  BIT,int  OFFSET, int  XPIX,int  YPIX,float  XLL,float  YLL,float  PPI,float *xur,float *yur, int ORIENT, int INVERT);
extern int vp_uraster_g (unsigned char *array, int BLAST,int  BIT,int  OFFSET, int  XPIX,int  YPIX,double  XLL,double  YLL,double  PPI,double *xur,double *yur, int ORIENT, int INVERT);
extern int vp_utext(float X,float Y,int SIZE,int ORIENT, char *string);
extern int vp_utext_g(double X,double Y,int SIZE,int ORIENT, char *string);
extern int vp_where(float*,float*);
extern int vp_where_g(double*,double*);
extern void p_pout (float xp,float  yp,int  down, FILE *plt);
extern int vp_fixpc(void);
extern int vp_fixpc_g(void);

/* handful of internal filter utility routines */
extern void wlimit(int,int,int*,int*);
extern void dithline (unsigned char *,unsigned char *, int,int,int);

#ifdef __cplusplus
}
#endif


#endif
