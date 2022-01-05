/*$

=head1 NAME

evaluate_expression - evaluate a mathematical expression

=head1 SYNOPSIS

C<int evaluate_expression(exp,get_val,nvals,result)>

=head1 INPUT PARAMETERS

=over 4

=item exp  - char*  

      mathematical expression

=item get_val - func 

      (get_val(char*,double*)) given a name function should
      return values

=item nvals  - int   

      number of values to calculate

=back

=head1 OUTPUT PARAMETERS

=over 4

=item result  - double* 

      expression result

=back

=head1 DESCRIPTION

User provides a routine (get_val) that is able to return double values
when given a string, evaluate_expression will evaluate a complicated
expression.

Supported Functions (specified by @ at begining):

COS	SIN	TAN
ACOS	ASIN	ATAN
COSH	SINH	INT
EXP	LOG	SQRT
ABS SGN @SCI:num:

=head1 SEE ALSO

L<Math>, L<Headermath>

=head1 LIBRARY

B<sep>

=cut
>*/
/*
Author:Robert G. Clapp

12/12/95 Begun
*/
#define IS_FLOAT 1
#define IS_INT 0
#define MAX_KEYS 40
#define KEEP_KEY 0
#define DELETE_KEY 1
#define MAX_EQN_LEN 1024
#define MAX_STR_LEN 128
#define FLOAT_TYPE "scalar_float"
#define INT_TYPE "scalar_int"
#define NOT_SET -1


#include<sepConfig.h>
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include<sep_main_external.h>
#include<math.h>
#include<ctype.h>
int first_zero_warn=0;

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
/*typedef int (*GET_VALS)(char *,double*);*/
int primitive( double *result);
int is_in(char ch,char *s);
int isdelim(char c,char *delim_string );
int get_token(void);
int unary(char *o, double *r);
int arith( double* r ,  double* h);
int level6(double *result);
int level5( double *result);
int level4(double *result);
int level3(double *result);
int level2(double *result);
_XFUNCPROTOEND
#else
int primitive();
int is_in();
int isdelim();
int get_token();
int unary();
int arith();
int level6();
int level5();
int level4();
int level3();
int level2();
#endif

char *equation;
int this_pass;
char token[4096];
char tok_type;
char op2[4096];
/*page 555 The Complete C Reference*/
#define DELIMITER 1
#define VARIABLE 2
#define NUMBER 3

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int (*get_oper)(char *, double*);
_XFUNCPROTOEND
#else
int (*get_oper)();
#endif



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int evaluate_expression(char *exp, int get_val(char *, double *), int nvals,double *result)
_XFUNCPROTOEND
#else
int evaluate_expression(exp,get_val,nvals,result)
char *exp; 
int get_val();
int nvals;
double *result;
#endif
{
int ierr;

equation=(char *) malloc(sizeof(char)*(strlen(exp)+2));
strcpy(equation,exp);
this_pass=nvals;
get_oper=get_val;


get_token();
if(!*token){
	ierr= seperr("not a valid equation %s \n",equation);

free(equation);
	return 0;
}
level2(result);
ierr=ierr;
return 0;
}


/*add or subract */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int level2(double *result)
_XFUNCPROTOEND
#else
int level2(result)
double *result;
#endif
{
register char op;
double *hold;
int ierr;
level3(result);
while ((op= *token) == '+' || op == '-') {
	get_token();
	hold=(double *)malloc(sizeof(double) * this_pass);
	level3(hold);
	sprintf(op2,"%c",op);
	ierr=arith(result,hold);
	free(hold);
	}
ierr=ierr;
return 0;
}



/*multiply or divide two factors*/
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int level3(double *result)
_XFUNCPROTOEND
#else
int level3(result)
double *result;
#endif
{
char op;
double *hold;
int ierr;
level4(result);
while ((op= *token) == '*' || op == '/' || op == '%') {
	get_token();
	hold=(double *)malloc(sizeof(double) * this_pass);
	level4(hold);
	sprintf(op2,"%c",op);
	ierr=arith(result,hold);
	free(hold);
	}
ierr=ierr;
return 0;
}

/* process an exponent */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int level4(double *result)
_XFUNCPROTOEND
#else
int level4(result)
double *result;
#endif
{
char op;
double *hold;
int ierr;

level5(result);
while ((op= *token) == '^') {
	get_token();
	hold=(double *)malloc(sizeof(double) * this_pass);
	level4(hold);
	sprintf(op2,"%c",'^');
	ierr=arith(result,hold);
	free(hold);
	}
ierr=ierr;
op=op;
return 0;
}


/* unary + or - */
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int level5( double *result)
_XFUNCPROTOEND
#else
int level5(result)
double *result;
#endif
{
/*register char op;*/
char op[16];
strcpy(op,"RDFG");
if ((tok_type== DELIMITER) && *token=='+' || *token == '-'){
	sprintf(op,"%c",*token);
	get_token();
}
else if((tok_type== DELIMITER) && *token=='@'){
	get_token();
	strcpy(op,token);
	get_token();
}
	level6(result);
	if(0!=strcmp(op,"RDFG")) unary(op,result);
return 0;
}

/*parenthesized equation*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int level6(double *result)
_XFUNCPROTOEND
#else
int level6(result)
double *result;
#endif
{
if ((tok_type== DELIMITER) && (*token=='(')){
	get_token();
	level2(result);
	if(*token != ')') seperr("bad equation %s \n",equation);
	get_token();
   }
else {
primitive(result);
}
return 0;
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int arith( double* r ,  double* h)
_XFUNCPROTOEND
#else
int arith( r ,  h)
double *r,*h;
#endif
{
register double t,ex;
int i1;
switch (op2[0]){

case '-':
	for(i1=0; i1< this_pass; i1++) r[i1]=r[i1]-h[i1];
	break;
case '+':
	for(i1=0; i1< this_pass; i1++) r[i1]=r[i1]+h[i1];
	break;
case '*':
	for(i1=0; i1< this_pass; i1++) r[i1]=r[i1] *  h[i1];
	break;
case '/':
	for(i1=0; i1< this_pass; i1++){
		if (h[i1] == 0){
	   if(first_zero_warn==0){
				fprintf(stderr,"warning dividing by 0 \n");
				first_zero_warn=1;
			}
	   	r[i1]=0;
		}
		else r[i1]=r[i1] / h[i1] ;
	}
	break;
case '%':
	for(i1=0; i1< this_pass; i1++)
	for(i1=0; i1< this_pass; i1++){
		t=r[i1] / h[i1] ;
		r[i1] = r[i1] - (t * h[i1]);
	}
	break;
case '^':
/* old
	for(i1=0; i1< this_pass; i1++){
		ex= r[i1];
		if(h[i1]==0){
			r[i1]=1;
		}
		else{
			for(t=h[i1]-1;t>0;--t) r[i1]=r[i1]*ex;
		}
	}
*/
	for(i1=0; i1< this_pass; i1++)
       r[i1]=pow(r[i1],(double)h[i1]);

	break;
}
return 0;
}
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int unary(char *o, double *r)
_XFUNCPROTOEND
#else
int unary(o,r)
char *o;
double *r;
#endif
{
int tempi,i1;


if (*o=='-') {
	for(i1=0; i1< this_pass; i1++)
	r[i1] =  r[i1] * -1.0;
}
else if (*o=='+') {
	for(i1=0; i1< this_pass; i1++) r[i1] =  r[i1];
}
else if (0==strcmp(o,"COS")){
	for(i1=0; i1< this_pass; i1++) r[i1] = cos( r[i1]);
}
else if (0==strcmp(o,"SIN")){
	for(i1=0; i1< this_pass; i1++) r[i1] = sin( r[i1]);
}
else if (0==strcmp(o,"TAN")){
	for(i1=0; i1< this_pass; i1++) r[i1] = tan( r[i1]);
}
else if (0==strcmp(o,"ACOS")){
	for(i1=0; i1< this_pass; i1++) {
		if(r[i1]>-1 && r[i1]<1) r[i1] = acos( r[i1]);
		else{
	   if(first_zero_warn==0){
				fprintf(stderr,"warning acos out of range setting to 0 %f \n",r[i1]);
				first_zero_warn=1;
			}
	   r[i1]=0.000;
		}
	}
}
else if (0==strcmp(o,"ASIN")){
	for(i1=0; i1< this_pass; i1++){
		if(r[i1]>-1 && r[i1]<1) r[i1] = asin( r[i1]);
		else{
	   if(first_zero_warn==0){
				fprintf(stderr,"warning acos out of range setting to 0 %f \n",r[i1]);
				first_zero_warn=1;
			}
	   r[i1]=0.000;
		}
	}
}
else if (0==strcmp(o,"ATAN")){
	for(i1=0; i1< this_pass; i1++) r[i1] = atan( r[i1]);
}
else if (0==strcmp(o,"INT")){
	for(i1=0; i1< this_pass; i1++){
		tempi=(int) r[i1];
		r[i1] = (double) tempi ;
	}
}
else if (0==strcmp(o,"COSH")){
	for(i1=0; i1< this_pass; i1++) r[i1] = cosh( r[i1]);
}
else if (0==strcmp(o,"SGN")){
	for(i1=0; i1< this_pass; i1++) r[i1] = r[i1]/fabs( r[i1]);
}
else if (0==strcmp(o,"ABS")){
	for(i1=0; i1< this_pass; i1++) r[i1] = fabs( r[i1]);
}
else if (0==strcmp(o,"EXP")){
	for(i1=0; i1< this_pass; i1++) r[i1] = exp ( r[i1]);
}
else if (0==strcmp(o,"LOG")){
	for(i1=0; i1< this_pass; i1++) r[i1] = log( r[i1]);
}
else if (0==strcmp(o,"SQRT")){
	for(i1=0; i1< this_pass; i1++) r[i1] = sqrt(r[i1]);
}
else seperr("Unknown function %s \n",o);
return 0;
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int get_token(void)
_XFUNCPROTOEND
#else
int get_token()
#endif
{
register char *temp;
char base1[10],base2[10],*base;
char b1;
tok_type=0;
temp=token;

strcpy(base1,"@+-%*/^()");
strcpy(base2,"@%*/^()");

if (is_in(*equation,base1)){
	tok_type=DELIMITER;
	*temp++ =*equation++;
}
/* 	ADD STUFF HERE FOR ADDITION FUNCTIONS-	*/
else if(isalpha(*equation)) {
	while(!isdelim(*equation,base1)) {
	*temp++=*equation++;
 	}
	tok_type=VARIABLE;
}
else if(isdigit(*equation)){
       b1='3';
       base=base1;
       while(!isdelim(*equation,base)){
	 *temp=*equation;
         b1=*equation;
         base=base1; 
         if(b1=='e') base=base2;
	 temp++; equation++;
	}
	tok_type=NUMBER;
}
*temp='\0';
return 0;
}




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int isdelim(char c,char *choices)
_XFUNCPROTOEND
#else
int isdelim(c,choices)
char c,*choices;
#endif
{
if (is_in(c,choices) ||(c==9) || (c == '\0') || (c==0))
	return 1;
return 0;
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int is_in(char ch,char *s)
_XFUNCPROTOBEGIN
#else
int is_in(ch,s)
char ch,*s;
#endif
{
while (*s) if (*s++==ch) return 1;
return 0;
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int primitive( double *result)
_XFUNCPROTOEND
#else
int primitive(result)
double *result;
#endif
{
int ierr;
double num;
int i1;

switch(tok_type){

case VARIABLE:
  	ierr=get_oper(token,result);
	get_token();
	   return 0 ;
case NUMBER:
	num=atof(token);
	for(i1=0; i1 < this_pass; i1++) result[i1]=num;
	get_token();
	return 0;
default:
	ierr=seperr("bad equation %s \n",equation);
}
ierr=ierr;
return 0;
}
