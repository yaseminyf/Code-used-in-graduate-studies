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

typedef char * STRING;
typedef struct MATRIX {                                  /* sparse matrix */   
               unsigned int rowdefined;
               int m,n;                        /* number of rows, columns */
               int nonzeros;                        /* number of nonzeros */
               int *Iind,*Jind;            /* row indexes, column indexes */
               int *rowptr,*colind;
               double *Valuecol,*Valuerow;          /* values of nonzeros */
               } MATRIX;

typedef struct VECTOR {                                    /* full vector */   
               int dim;                            /* dimension of vector */
               double *ve;                                      /* values */
               } VECTOR;

typedef struct FMATRIX {                                   /* full matrix */
               int m,n;                        /* number of rows, columns */
               double *Value;                       /* values of nonzeros */
               } FMATRIX;

typedef struct GROUP {      /* to keep track of the column groups */ 
               int num_columns; /* number of columns in a group */ 
               int first,last;  /* the beginning and ending columns */
               struct GROUP *next;
               } GROUP;

/* structures borrowed from SPARSPAK to be used in C code */

typedef struct comtype4 {
              int MSGLVB,IERRB,MAXSB,MCOLS,MSEQNS,MDEQNS,MSCONS,MDCONS;
              }comtype4;
struct comtype4 spbusr_;
typedef struct comtype5 {
              int STAGE , IOUNIT, MXUSED, MXREQD, NCOLS;
              int NSEQNS, NDEQNS, NSCONS, NDCONS, NZEQNS;
              int NZCONS, NZMAX , NEDGES, METHOD, NOFNZ;
              int NOFSUB, ICPAD[28],RKDEF1,RKDEF2,NSHORT,NLONG,NSP1,NLP1;
              }comtype5;
struct comtype5 spbcon_;
typedef struct comtype6 {
              int PERM,INVP,RHS,DIAG,XRNZ;
              int RNZ,XNZSUB,NZSUB,DSEQNS,DSBEQN;
              int DSCONS,DSBCON,ROWMSK,XZROWS,IMPAD[35],LLFREE;
              }comtype6;
struct comtype6 spbmap_;
typedef struct comtype1 { 
              int IPRNTE,IPRNTS,MAXINT; 
              float RATIOS,RATIOL,MCHEPS,TIME; 
              }comtype1; 
struct comtype1 spksys_;

/* new struct definition for some global variables */

typedef struct comtype2 {
              int mdeg, ehead, tag, mdnode;           /* genmmd.c */
              }comtype2;
struct comtype2 spkglob_;
              
/* extern FORTRAN codes translated to C */

extern void inxywb();
extern void orcolb();
extern void lsqslvsn3();
extern VECTOR *solution();
extern VECTOR *gsslv();

extern void adaij5();
extern void build();
extern void copysi();
extern void fmadjy();
extern void genls3();
extern void genmmd();
extern void gsfct();
extern void ipendb();
extern void mmdelm();
extern void mmdint();
extern void mmdnum();
extern void mmdupd();
extern void rcopyl();
extern void rdeqns();
extern void rkchk2();
extern void rwprep();
extern void smbfct();
extern void sorts1();
extern void zerols();
extern void zerorv();

/* extern C codes */

extern double *multAtransr();
extern VECTOR *m_fulvectormlt();
extern VECTOR *multAtransr2();
extern void readmatrix();
extern void readbigmatrix();
extern MATRIX **create_over();
extern MATRIX **create_nover();
extern int *initialize();
extern MATRIX **create_nover25();
extern MATRIX **create_over_defined();
extern MATRIX **create_dependent();
extern MATRIX **create_over_defined_new();
extern MATRIX **create_nover_diff();
extern void readnewmatrix();
extern VECTOR *multAtransr3();
extern VECTOR *v_get();
extern VECTOR *vones();
extern void sym_mult();
extern void sym_mult2();
extern void sym_mult3();
extern VECTOR *permdv();
extern double norm2();
extern VECTOR *v_zero();
extern void permuteA();
extern void multAD();
extern void createC();
extern void createC1();
extern void createC2();
extern void createC22();
extern void createC222();
extern MATRIX **createtildeAptr();
extern MATRIX *createtildeA();
extern void putvaluesintildeA();
extern void putvaluesintildeA3();

/* global variables */

/* variables to be used in SPARSPAK routines */

double TOL;
int TYPTOL,OPTION;
double Rhs,WEIGHT;
int TYPE,ROWNUM,NSUBS;
int *SUBS,*T;
double *VALUES;
int iladj, limit;                                 /* ipendb.c */
int nudeg, nuladj, dhead, qsize, marker, llist;   /* orcolb.c */
int invp, perm, newadj;                           /* orcolb.c */
int q, mask, clqtst, nzsub, xrnz, xnzsub, nofsub; /* orcolb.c */
int flag__, maxlnz, maxsub;                       /* orcolb.c */
double *ddiag, *rnzv, *rhsv;                      /* lsqslvsn3.c */
int *nzsubv, *xrnzv, *xnzsubv, *invpv, *permv;    /* lsqslvsn3.c */
int *rowmsk;                                      /* lsqslvsn3.c */
double *qb;                                       /* lsqslvsn3.c */
double *lnzv;                                     /* gsfct.c */
int *llink, *first;                               /* gsfct.c */ 
double *temp;                                     /* gsfct.c */  
double b;                                         /* rkchk2.c */
int *xzrows;                                      /* lsqslvsn3.c */
int nsubs;                                        /* rkchk2.c */
int *nusubs;
double *nuvals;
int delta;                                        /* genmmd.c */
int adjncy, ladj, xadj, deg;                      /* ipendb.c */
int *isub,*jsub;                                  /* sym_mult.c */
double *values;                                   /* sym_mult.c */

/* other variables */

VECTOR *s,*x,*rhs,*xexact,*sol,*residual,*rvec2;

MATRIX *A;

FILE *fp,*fp1;

MATRIX **sub_ptr;

int ngroup;
int **zeroptr,**groupptr;

int **xrnzptr,**nzsubptr,**xnzsubptr,**permptr,**invpptr;
double **ddiagptr,**rnzptr;
double r_zero_norm, r_norm, err;
double epsilon, omega;
int *zeroptr_par;
