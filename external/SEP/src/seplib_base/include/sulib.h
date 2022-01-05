#ifndef SEP_SULIB_H
#define SEP_SULIB_H yada
/*this should be changed - don't need to pull in all of the includes that
this implies */
#include<prototypes.h> 
#include<sepmath.h>
#include <math.h>
#ifndef PI
static double sxftEkd=3.14159265358979323846264338327950288419716939937510;
#define PI sxftEkd
#endif
#ifndef pi
static double snftEkd=3.14159265358979323846264338327950288419716939937510;
#define pi snftEkd
#endif

#ifndef __cplusplus  /* C++ has its own complex arithmetic class */
#ifndef OLDcomplex
typedef struct { float re, im;} d0u1m2m3y4cmplx;
#define OLDcomplex d0u1m2m3y4cmplx
#endif
#endif

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef NINT
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef dcomplex
typedef struct _dcomplexStruct { /* double-precision complex number */
  double re,im;
} dcomplex;
#endif/* dcomplex */





#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
/*extern int npfa (int);
extern int npfao (int , int );
*/
extern int npfar (int );
extern int npfaro (int , int );
extern void pfacc(int , int, OLDcomplex*);
extern void pfacr (int, int, OLDcomplex*,float*);
extern void pfarc (int, int, float*, float*);
extern void pfamcc (int, int, int, int, int kt, OLDcomplex*);
extern void pfa2cc (int, int, int, int, OLDcomplex*);
extern void pfa2cr (int, int, int, int, OLDcomplex*, float*);
extern void pfa2rc (int, int, int, int, float*, OLDcomplex*);
_XFUNCPROTOEND
#else
extern int npfa ();
extern int npfao ();
extern int npfar ();
extern int npfaro ();
extern void pfacc();
extern void pfacr ();
extern void pfarc ();
extern void pfamcc ();
extern void pfa2cc ();
extern void pfa2cr ();
extern void pfa2rc ();
#endif


#include "prototypes.h"
#include<stdio.h>
#include<stdlib.h>
#include<rpc/types.h>
#include<rpc/xdr.h>

#ifndef HDR_H
#define HDR_H
static struct {
	char *key;	char *type;	int offs;
} hdr[] = {
	{   "tracl",		"i",		0},
	{   "tracr",		"i",		4},
	{    "fldr",		"i",		8},
	{   "tracf",		"i",		12},
	{      "ep",		"i",		16},
	{     "cdp",		"i",		20},
	{    "cdpt",		"i",		24},
	{    "trid",		"h",		28},
	{     "nvs",		"h",		30},
	{     "nhs",		"h",		32},
	{    "duse",		"h",		34},
	{  "offset",		"i",		36},
	{   "gelev",		"i",		40},
	{   "selev",		"i",		44},
	{  "sdepth",		"i",		48},
	{    "gdel",		"i",		52},
	{    "sdel",		"i",		56},
	{   "swdep",		"i",		60},
	{   "gwdep",		"i",		64},
	{  "scalel",		"h",		68},
	{  "scalco",		"h",		70},
	{      "sx",		"i",		72},
	{      "sy",		"i",		76},
	{      "gx",		"i",		80},
	{      "gy",		"i",		84},
	{  "counit",		"h",		88},
	{   "wevel",		"h",		90},
	{  "swevel",		"h",		92},
	{     "sut",		"h",		94},
	{     "gut",		"h",		96},
	{   "sstat",		"h",		98},
	{   "gstat",		"h",		100},
	{   "tstat",		"h",		102},
	{    "laga",		"h",		104},
	{    "lagb",		"h",		106},
	{   "delrt",		"h",		108},
	{    "muts",		"h",		110},
	{    "mute",		"h",		112},
	{      "ns",		"u",		114},
	{      "dt",		"u",		116},
	{    "gain",		"h",		118},
	{     "igc",		"h",		120},
	{     "igi",		"h",		122},
	{    "corr",		"h",		124},
	{     "sfs",		"h",		126},
	{     "sfe",		"h",		128},
	{    "slen",		"h",		130},
	{    "styp",		"h",		132},
	{    "stas",		"h",		134},
	{    "stae",		"h",		136},
	{   "tatyp",		"h",		138},
	{   "afilf",		"h",		140},
	{   "afils",		"h",		142},
	{  "nofilf",		"h",		144},
	{  "nofils",		"h",		146},
	{     "lcf",		"h",		148},
	{     "hcf",		"h",		150},
	{     "lcs",		"h",		152},
	{     "hcs",		"h",		154},
	{    "year",		"h",		156},
	{     "day",		"h",		158},
	{    "hour",		"h",		160},
	{  "minute",		"h",		162},
	{     "sec",		"h",		164},
	{  "timbas",		"h",		166},
	{    "trwf",		"h",		168},
	{  "grnors",		"h",		170},
	{  "grnofr",		"h",		172},
	{  "grnlof",		"h",		174},
	{    "gaps",		"h",		176},
	{   "otrav",		"h",		178},
	{      "d1",		"f",		180},
	{      "f1",		"f",		184},
	{      "d2",		"f",		188},
	{      "f2",		"f",		192},
	{  "ungpow",		"f",		196},
	{ "unscale",		"f",		200},
	{     "ntr",		"i",		204},
	{    "mark",		"h",		208},
	{"shortpad",		"h",		210},
};
#endif
/*
 * header.h - include file for segy sizes
 * THIS HEADER FILE IS GENERATED AUTOMATICALLY - 
 * see the makefile in this directory
 */

#ifndef HEADER_H
#define HEADER_H

#define SU_NKEYS	80	/* Number of key header words */
#define HDRBYTES	240	/* Bytes in the trace header */
#define	MAXSEGY		131312

#endif
typedef union {
        char s[8];
        short h;
        unsigned short u;
        long l;
        unsigned long v;
        int i;
        unsigned int p;
        float f;
        double d;
} Value;

/* Copyright (c) Colorado School of Mines, 1997.*/
/* All rights reserved.                       */

/* Copyright (c) Colorado School of Mines, 1994.*/
/* All rights reserved.                       */

/* segy.h - include file for SEGY traces
 *
 * declarations for:
 *	typedef struct {} segy - the trace identification header
 *	typedef struct {} bhed - binary header
 *
 * Note:
 *	If header words are added, run the makefile in this directory
 *	to recreate hdr.h.
 *
 * Reference:
 *	K. M. Barry, D. A. Cavers and C. W. Kneale, "Special Report:
 *		Recommended Standards for Digital Tape Formats",
 *		Geophysics, vol. 40, no. 2 (April 1975), P. 344-352.
 *	
 * $Author: cvs $
 * $Source: /sepwww2/bob/SEP_CVS/seplib/seplib-6.3.5/seplib_base/include/sulib.h,v $
 * $Revision: 1.1.1.1 $ ; $Date: 2004/03/25 06:37:24 $
 */ 

#include "prototypes.h"
#ifndef SEGY_H
#define SEGY_H

#define SU_NFLTS	32768	/* Arbitrary limit on data array size	*/


/* TYPEDEFS */
typedef struct {	/* segy - trace identification header */

	int tracl;	/* trace sequence number within line */

	int tracr;	/* trace sequence number within reel */

	int fldr;	/* field record number */

	int tracf;	/* trace number within field record */

	int ep;	/* energy source point number */

	int cdp;	/* CDP ensemble number */

	int cdpt;	/* trace number within CDP ensemble */

	short trid;	/* trace identification code:
			1 = seismic data
			2 = dead
			3 = dummy
			4 = time break
			5 = uphole
			6 = sweep
			7 = timing
			8 = water break
			9---, N = optional use (N = 32,767)

			Following are CWP id flags:

			 9 = autocorrelation

			10 = Fourier transformed - no packing
			     xr[0],xi[0], ..., xr[N-1],xi[N-1]

			11 = Fourier transformed - unpacked Nyquist
			     xr[0],xi[0],...,xr[N/2],xi[N/2]

			12 = Fourier transformed - packed Nyquist
	 		     even N:
			     xr[0],xr[N/2],xr[1],xi[1], ...,
				xr[N/2 -1],xi[N/2 -1]
				(note the exceptional second entry)
			     odd N:
			     xr[0],xr[(N-1)/2],xr[1],xi[1], ...,
				xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
				(note the exceptional second & last entries)

			13 = Complex signal in the time domain
			     xr[0],xi[0], ..., xr[N-1],xi[N-1]

			14 = Fourier transformed - amplitude/phase
			     a[0],p[0], ..., a[N-1],p[N-1]

			15 = Complex time signal - amplitude/phase
			     a[0],p[0], ..., a[N-1],p[N-1]

			16 = Real part of complex trace from 0 to Nyquist

			17 = Imag part of complex trace from 0 to Nyquist

			18 = Amplitude of complex trace from 0 to Nyquist

			19 = Phase of complex trace from 0 to Nyquist

			21 = Wavenumber time domain (k-t)

			22 = Wavenumber frequency (k-omega)

			23 = Envelope of the complex time trace

			24 = Phase of the complex time trace

			25 = Frequency of the complex time trace

			30 = Depth-Range (z-x) traces

			101 = Seismic data packed to bytes (by supack1)
			
			102 = Seismic data packed to 2 bytes (by supack2)
			*/

	short nvs;	/* number of vertically summed traces (see vscode
			   in bhed structure) */

	short nhs;	/* number of horizontally summed traces (see vscode
			   in bhed structure) */

	short duse;	/* data use:
				1 = production
				2 = test */

	int offset;	/* distance from source point to receiver
			   group (negative if opposite to direction
			   in which the line was shot) */

	int gelev;	/* receiver group elevation from sea level
			   (above sea level is positive) */

	int selev;	/* source elevation from sea level
			   (above sea level is positive) */

	int sdepth;	/* source depth (positive) */

	int gdel;	/* datum elevation at receiver group */

	int sdel;	/* datum elevation at source */

	int swdep;	/* water depth at source */

	int gwdep;	/* water depth at receiver group */

	short scalel;	/* scale factor for previous 7 entries
			   with value plus or minus 10 to the
			   power 0, 1, 2, 3, or 4 (if positive,
			   multiply, if negative divide) */

	short scalco;	/* scale factor for next 4 entries
			   with value plus or minus 10 to the
			   power 0, 1, 2, 3, or 4 (if positive,
			   multiply, if negative divide) */

	int  sx;	/* X source coordinate */

	int  sy;	/* Y source coordinate */

	int  gx;	/* X group coordinate */

	int  gy;	/* Y group coordinate */

	short counit;	/* coordinate units code:
				for previous four entries
				1 = length (meters or feet)
				2 = seconds of arc (in this case, the
				X values are longitude and the Y values
				are latitude, a positive value designates
				the number of seconds east of Greenwich
				or north of the equator */

	short wevel;	/* weathering velocity */

	short swevel;	/* subweathering velocity */

	short sut;	/* uphole time at source */

	short gut;	/* uphole time at receiver group */

	short sstat;	/* source static correction */

	short gstat;	/* group static correction */

	short tstat;	/* total static applied */

	short laga;	/* lag time A, time in ms between end of 240-
			   byte trace identification header and time
			   break, positive if time break occurs after
			   end of header, time break is defined as
			   the initiation pulse which maybe recorded
			   on an auxiliary trace or as otherwise
			   specified by the recording system */

	short lagb;	/* lag time B, time in ms between the time break
			   and the initiation time of the energy source,
			   may be positive or negative */

	short delrt;	/* delay recording time, time in ms between
			   initiation time of energy source and time
			   when recording of data samples begins
			   (for deep water work if recording does not
			   start at zero time) */

	short muts;	/* mute time--start */

	short mute;	/* mute time--end */

	unsigned short ns;	/* number of samples in this trace */

	unsigned short dt;	/* sample interval; in micro-seconds */

	short gain;	/* gain type of field instruments code:
				1 = fixed
				2 = binary
				3 = floating point
				4 ---- N = optional use */

	short igc;	/* instrument gain constant */

	short igi;	/* instrument early or initial gain */

	short corr;	/* correlated:
				1 = no
				2 = yes */

	short sfs;	/* sweep frequency at start */

	short sfe;	/* sweep frequency at end */

	short slen;	/* sweep length in ms */

	short styp;	/* sweep type code:
				1 = linear
				2 = cos-squared
				3 = other */

	short stas;	/* sweep trace length at start in ms */

	short stae;	/* sweep trace length at end in ms */

	short tatyp;	/* taper type: 1=linear, 2=cos^2, 3=other */

	short afilf;	/* alias filter frequency if used */

	short afils;	/* alias filter slope */

	short nofilf;	/* notch filter frequency if used */

	short nofils;	/* notch filter slope */

	short lcf;	/* low cut frequency if used */

	short hcf;	/* high cut frequncy if used */

	short lcs;	/* low cut slope */

	short hcs;	/* high cut slope */

	short year;	/* year data recorded */

	short day;	/* day of year */

	short hour;	/* hour of day (24 hour clock) */

	short minute;	/* minute of hour */

	short sec;	/* second of minute */

	short timbas;	/* time basis code:
				1 = local
				2 = GMT
				3 = other */

	short trwf;	/* trace weighting factor, defined as 1/2^N
			   volts for the least sigificant bit */

	short grnors;	/* geophone group number of roll switch
			   position one */

	short grnofr;	/* geophone group number of trace one within
			   original field record */

	short grnlof;	/* geophone group number of last trace within
			   original field record */

	short gaps;	/* gap size (total number of groups dropped) */

	short otrav;	/* overtravel taper code:
				1 = down (or behind)
				2 = up (or ahead) */

	/* local assignments */
	float d1;	/* sample spacing for non-seismic data */

	float f1;	/* first sample location for non-seismic data */

	float d2;	/* sample spacing between traces */

	float f2;	/* first trace location */

	float ungpow;	/* negative of power used for dynamic
			   range compression */

	float unscale;	/* reciprocal of scaling factor to normalize
			   range */

	int ntr; 	/* number of traces */

	short mark;	/* mark selected traces */

        short shortpad; /* alignment padding */


	short unass[14];	/* unassigned--NOTE: last entry causes 
			   a break in the word alignment, if we REALLY
			   want to maintain 240 bytes, the following
			   entry should be an odd number of short/UINT2
			   OR do the insertion above the "mark" keyword
			   entry */

	float  data[SU_NFLTS];

} segy;


typedef struct {	/* bhed - binary header */

	int jobid;	/* job identification number */

	int lino;	/* line number (only one line per reel) */

	int reno;	/* reel number */

	short ntrpr;	/* number of data traces per record */

        short nart;	/* number of auxiliary traces per record */

	unsigned short hdt; /* sample interval in micro secs for this reel */

	unsigned short dto; /* same for original field recording */

	unsigned short hns; /* number of samples per trace for this reel */

	unsigned short nso; /* same for original field recording */

	short format;	/* data sample format code:
				1 = floating point (4 bytes)
				2 = fixed point (4 bytes)
				3 = fixed point (2 bytes)
				4 = fixed point w/gain code (4 bytes) */

	short fold;	/* CDP fold expected per CDP ensemble */

	short tsort;	/* trace sorting code: 
				1 = as recorded (no sorting)
				2 = CDP ensemble
				3 = single fold continuous profile
				4 = horizontally stacked */

	short vscode;	/* vertical sum code:
				1 = no sum
				2 = two sum ...
				N = N sum (N = 32,767) */

	short hsfs;	/* sweep frequency at start */

	short hsfe;	/* sweep frequency at end */

	short hslen;	/* sweep length (ms) */

	short hstyp;	/* sweep type code:
				1 = linear
				2 = parabolic
				3 = exponential
				4 = other */

	short schn;	/* trace number of sweep channel */

	short hstas;	/* sweep trace taper length at start if
			   tapered (the taper starts at zero time
			   and is effective for this length) */

	short hstae;	/* sweep trace taper length at end (the ending
			   taper starts at sweep length minus the taper
			   length at end) */

	short htatyp;	/* sweep trace taper type code:
				1 = linear
				2 = cos-squared
				3 = other */

	short hcorr;	/* correlated data traces code:
				1 = no
				2 = yes */

	short bgrcv;	/* binary gain recovered code:
				1 = yes
				2 = no */

	short rcvm;	/* amplitude recovery method code:
				1 = none
				2 = spherical divergence
				3 = AGC
				4 = other */

	short mfeet;	/* measurement system code:
				1 = meters
				2 = feet */

	short polyt;	/* impulse signal polarity code:
				1 = increase in pressure or upward
				    geophone case movement gives
				    negative number on tape
				2 = increase in pressure or upward
				    geophone case movement gives
				    positive number on tape */

	short vpol;	/* vibratory polarity code:
				code	seismic signal lags pilot by
				1	337.5 to  22.5 degrees
				2	 22.5 to  67.5 degrees
				3	 67.5 to 112.5 degrees
				4	112.5 to 157.5 degrees
				5	157.5 to 202.5 degrees
				6	202.5 to 247.5 degrees
				7	247.5 to 292.5 degrees
				8	293.5 to 337.5 degrees */

	short hunass[170];	/* unassigned */

} bhed;

/* The following refer to the trid field in segy.h		*/
/* CHARPACK represents byte packed seismic data from supack1	*/
#define		CHARPACK	101
/* SHORTPACK represents 2 byte packed seismic data from supack2	*/
#define		SHORTPACK	102

/* TREAL represents real time traces 				*/
#define		TREAL		1
/* TDEAD represents dead time traces 				*/
#define		TDEAD		2
/* TDUMMY represents dummy time traces 				*/
#define		TDUMMY		3
/* TBREAK represents time break traces 				*/
#define		TBREAK		4
/* UPHOLE represents uphole traces 				*/
#define		UPHOLE		5
/* SWEEP represents sweep traces 				*/
#define		SWEEP		6
/* TIMING represents timing traces 				*/
#define		TIMING		7
/* WBREAK represents timing traces 				*/
#define		WBREAK		8

/* TCMPLX represents complex time traces 			*/
#define		TCMPLX		13
/* TAMPH represents time domain data in amplitude/phase form	*/
#define		TAMPH		15
/* FPACK represents packed frequency domain data 		*/
#define		FPACK		12
/* FUNPACKNYQ represents complex frequency domain data 		*/
#define		FUNPACKNYQ	11
/* FCMPLX represents complex frequency domain data 		*/
#define		FCMPLX		10
/* FAMPH represents freq domain data in amplitude/phase form	*/
#define		FAMPH		14
/* REALPART represents the real part of a trace to Nyquist	*/
#define		REALPART	16
/* IMAGPART represents the real part of a trace to Nyquist	*/
#define		IMAGPART	17
/* AMPLITUDE represents the amplitude of a trace to Nyquist	*/
#define		AMPLITUDE	18
/* PHASE represents the phase of a trace to Nyquist		*/
#define		PHASE		19
/* KT represents wavenumber-time domain data 			*/
#define		KT		21
/* KOMEGA represents wavenumber-frequency domain data		*/
#define		KOMEGA		22
/* ENVELOPE represents the envelope of the complex time trace	*/
#define		ENVELOPE	23
/* INSTPHASE represents the phase of the complex time trace	*/
#define		INSTPHASE	24
/* INSTFREQ represents the frequency of the complex time trace	*/
#define		INSTFREQ	25
/* DEPTH represents traces in depth-range (z-x)			*/
#define		TRID_DEPTH	30

#define ISSEISMIC(id) (( (id)==0 || (id)==TREAL || (id)==TDEAD || (id)==TDUMMY ) ? cwp_true : cwp_false )

/* FUNCTION PROTOTYPES */
#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif
int fgettr(FILE *fp, segy *tp);
int fvgettr(FILE *fp, segy *tp);
void fputtr(FILE *fp, segy *tp);
int fgettra(FILE *fp, segy *tp, int itr);

/* hdrpkge */
char *hdtype(const char *key);
char *getkey(const int index);
int getindex(const char *key);
void swaphval(segy *tp, int index);
void swapbhval(bhed *bhp, int index);
void printheader(const segy *tp);
/*int xdrbhdrsub(XDR *segyxdr, bhed *bhp);*/

void tabplot(segy *tp, int itmin, int itmax);

#ifdef __cplusplus /* if C++, end external linkage specification */
}
#endif

#endif
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern int xdrbhdrsub(XDR*, bhed*);
extern int xdrhdrsub(XDR *segyxdr, segy *tp);
extern int vtoi(register char* type, Value val);
long vtol(register char* type, Value val);
float vtof(register char* type, Value val);
double vtod(register char* type, Value val);
void atoval(char* type, char* keyval, Value *valp);
Value valtoabs(char* type, Value val);
int valcmp(register char* type, Value val1, Value val2);
void printfval(register char* type, Value val);
void fprintfval(FILE *stream, register char* type, Value val);
void scanfval(register char* type, Value *valp);
void *alloc1 (size_t n1, size_t size);
void *realloc1 (void *v, size_t n1, size_t size);
void **alloc2 (size_t n1, size_t n2, size_t size);
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size);
void ****alloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void *****alloc5 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size);
void ******alloc6 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6,
                   size_t size);

void free1 (void *p);
void free2 (void **p);
void free3 (void ***p);
void free4 (void ****p);
void free5 (void *****p);
void free6 (void ******p);
int *alloc1int (size_t n1);
int *realloc1int (int *v, size_t n1);
int **alloc2int (size_t n1, size_t n2);
int ***alloc3int (size_t n1, size_t n2, size_t n3);
float *alloc1float (size_t n1);
float *realloc1float (float *v, size_t n1);
float **alloc2float (size_t n1, size_t n2);
float ***alloc3float (size_t n1, size_t n2, size_t n3);

float ****alloc4float (size_t n1, size_t n2, size_t n3, size_t n4);
void free4float (float ****p);
float *****alloc5float (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
void free5float (float *****p);
float ******alloc6float (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6);
void free6float (float ******p);
int ****alloc4int (size_t n1, size_t n2, size_t n3, size_t n4);
void free4int (int ****p);
int *****alloc5int (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
void free5int (int *****p);
unsigned short ******alloc6ushort(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5, size_t n6);
unsigned char *****alloc5uchar(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5);
void free5uchar(unsigned char *****p);
unsigned short *****alloc5ushort(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5);
void free5ushort(unsigned short *****p);
unsigned char ******alloc6uchar(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5, size_t n6);
void free6uchar(unsigned char ******p);
unsigned short ******alloc6ushort(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5, size_t n6);
void free6ushort(unsigned short ******p);

double *alloc1double (size_t n1);
double *realloc1double (double *v, size_t n1);
double **alloc2double (size_t n1, size_t n2);
double ***alloc3double (size_t n1, size_t n2, size_t n3);
OLDcomplex *alloc1complex (size_t n1);
OLDcomplex *realloc1complex (OLDcomplex *v, size_t n1);
OLDcomplex **alloc2complex (size_t n1, size_t n2);
OLDcomplex ***alloc3complex (size_t n1, size_t n2, size_t n3);

dcomplex *alloc1dcomplex (size_t n1);
dcomplex *realloc1dcomplex (dcomplex *v, size_t n1);
dcomplex **alloc2dcomplex (size_t n1, size_t n2);
dcomplex ***alloc3dcomplex (size_t n1, size_t n2, size_t n3);

void free1int (int *p);
void free2int (int **p);
void free3int (int ***p);
void free1float (float *p);
void free2float (float **p);
void free3float (float ***p);

void free1double (double *p);
void free2double (double **p);
void free3double (double ***p);
void free1complex (OLDcomplex *p);
void free2complex (OLDcomplex **p);
void free3complex (OLDcomplex ***p);

void free1dcomplex (dcomplex *p);
void free2dcomplex (dcomplex **p);
void free3dcomplex (dcomplex ***p);
void *ealloc1 (size_t n1, size_t size);
void *erealloc1 (void *v, size_t n1, size_t size);
void **ealloc2 (size_t n1, size_t n2, size_t size);
void ***ealloc3 (size_t n1, size_t n2, size_t n3, size_t size);
void ****ealloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void ****ealloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void *****ealloc5 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size);
void ******ealloc6 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5,
                  size_t n6, size_t size);

int *ealloc1int(size_t n1);
int *erealloc1int(int *v, size_t n1);
int **ealloc2int(size_t n1, size_t n2);
int ***ealloc3int(size_t n1, size_t n2, size_t n3);
float *ealloc1float(size_t n1);
float *erealloc1float(float *v, size_t n1);
float **ealloc2float(size_t n1, size_t n2);
float ***ealloc3float(size_t n1, size_t n2, size_t n3);

int ****ealloc4int(size_t n1, size_t n2, size_t n3, size_t n4);
int *****ealloc5int(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
float ****ealloc4float(size_t n1, size_t n2, size_t n3, size_t n4);
float *****ealloc5float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
float ******ealloc6float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5,
       size_t n6);

unsigned short *****ealloc5ushort(size_t n1, size_t n2,
     size_t n3, size_t n4, size_t n5);
unsigned char *****ealloc5uchar(size_t n1, size_t n2,
    size_t n3, size_t n4, size_t n5);
unsigned short ******ealloc6ushort(size_t n1, size_t n2,
     size_t n3, size_t n4, size_t n5, size_t n6);

double *ealloc1double(size_t n1);
double *erealloc1double(double *v, size_t n1);
double **ealloc2double(size_t n1, size_t n2);
double ***ealloc3double(size_t n1, size_t n2, size_t n3);
OLDcomplex *ealloc1complex(size_t n1);
OLDcomplex *erealloc1complex(OLDcomplex *v, size_t n1);
OLDcomplex **ealloc2complex(size_t n1, size_t n2);
OLDcomplex ***ealloc3complex(size_t n1, size_t n2, size_t n3);
void intcub (int ideriv, int nin, float xin[], float ydin[][4],
  int nout, float xout[], float yout[]);
void cakima (int n, float x[], float y[], float yd[][4]);
void cmonot (int n, float x[], float y[], float yd[][4]);
void csplin (int n, float x[], float y[], float yd[][4]);
void xindex (int nx, float ax[], float x, int *index);
_XFUNCPROTOEND
#else
extern int xdrbhdrsub();
extern int xdrhdrsub();
int vtoi();
long vtol();
float vtof();
double vtod();
void atoval();
Value valtoabs();
int valcmp();
void printfval();
void fprintfval();
void scanfval();
#endif

#ifdef __cplusplus
}
#endif
#endif
