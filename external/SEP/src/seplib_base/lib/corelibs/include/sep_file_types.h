#ifndef SEP_SIZE
#define SEP_SIZE 1
/*this is to get around a problem on SUN 4 OS and LINUX who have no*/
/*concept of a 64 bit int, therefore placing a limitation of 2 GB on*/
/*file size (not totally correct on LINUX) sep_size corresponds to*/
/*a number that can handle integers of + 2GB.  sep_off_t is related*/
/*to off_t, basically a pointer into to a file*/
/*don't care anymore long long  is valid on any recent ARCH-RGC*/

/*typedef long sep_off_t;*/
typedef long long sep_off_t;
typedef long long sep_file_size_t;
	

#define MAX_INT_SIZE 2147483647

#endif
