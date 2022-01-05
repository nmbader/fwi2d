#include<cfortran.h>
#include<sepmath.h>

#define QUANTILE(k,n,bb,value) CCALLSFSUB4(QUANTILE,quantile,INT,INT,FLOATV,FLOAT,k,n,bb,value)
#define RANDM(randum) CCALLSFSUB1(RANDM,randm,FLOAT,randum)
#define RSEED(seedin) CCALLSFSUB1(RSEED,rseed,INT,seedin)
#define GRAND(gr) CCALLSFSUB1(GRAND,grand,FLOAT,gr)
/*PROTOCCALLSFFUN1(FLOAT,RAND01,rand01,INT)*/
#define RAND01(iseed) CCALLSFFUN1(RAND01,rand01,INT,iseed)
