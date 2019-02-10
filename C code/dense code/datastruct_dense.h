#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define CREATE(Size,Type) (Type *)calloc(Size,sizeof(Type))
#define TRUE 1
#define FALSE 0
#define MIN(a,b) ( a<b ? a : b )
#define MAX(a,b) ( a>=b ? a : b )
#define ABS(a) (a<0.0 ? -a : a)

typedef struct VECTOR {                                    /* full vector */   
               int dim;                            /* dimension of vector */
               double *ve;                                      /* values */
               } VECTOR;

typedef struct FMATRIX {                                   /* full matrix */
               int m,n;                        /* number of rows, columns */
               double *Value;                       /* values of nonzeros */
               } FMATRIX;
               
extern void readdense();
extern VECTOR *v_get();
extern FMATRIX **create_dense_blocks();
extern double norm2();
extern VECTOR *vones();
extern FMATRIX *createtildeA();
extern void putvaluesintildeA();

FILE *fp,*fp1;

FMATRIX **sub_ptr,*A,*AA,**Rptr,*B,*C,*CC,*tildeA,**tildeAptr;

VECTOR *rhs,*residual,*s,*sol,*x,*xexact;  /* vectors used in computations */

VECTOR **D,**Ad,*pvec,**pvecptr,*littles,*vvec,*zvec,**Apptr,**tmpApptr;

double r_zero_norm, r_norm, err;
double epsilon, omega;

int ngroup,numnz;
