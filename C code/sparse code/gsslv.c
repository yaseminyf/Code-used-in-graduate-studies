/* last update: 03.09.1998 */
/* gsslv called by main program */

#include "data_struct.h"  

/* ******     gsslv ..... general sparse symmetric solve ******* */
/*     purpose - to perform solution of a factored system, where */
/*        the matrix is stored in                                */
/*        the compressed subscript data format                   */
/*     input parameters -                                        */
/*        neqns - number of equations                            */
/*        (xrnzv,lnzv) - structure of nonzeros in L.             */
/*        (xnzsubv,nzsubv) - compressed subscript structure      */
/*     updated parameters -                                      */
/*        s - on input it the rhs vector, and on                 */
/*              output the solution vector.                      */
/* ************************************************************* */

VECTOR *gsslv(neqns, xrnzv, lnzv, xnzsubv, nzsubv, ddiag, s)
int neqns, *xrnzv;
double *lnzv;
int *xnzsubv, *nzsubv;
double *ddiag;
VECTOR *s;
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int isub;
    static double rhsj;
    static int i__, j;
    static double ss;
    static int istop, istrt, ii, jj;

/* *************************************************************** */

/*        ------------------       */
/*        forward substitution ... */
/*        ------------------       */

    /* Function Body */
    
#ifdef DEBUG
   fprintf(fp1,"\n inside gsslv..");
#endif    

    i__1 = neqns;
    for (j = 1; j <= i__1; j++) {
    
#ifdef DEBUG
printf("\n inside loop..");
#endif
    
	rhsj = s->ve[j-1] / ddiag[j];
	s->ve[j-1] = rhsj;
	istrt = xrnzv[j];
	istop = xrnzv[j + 1] - 1;
	if (istop < istrt) {
	    goto L200;
	}
	i__ = xnzsubv[j];
	i__2 = istop;
	for (ii = istrt; ii <= i__2; ii++) {
	    isub = nzsubv[i__];
	    s->ve[isub-1] -= lnzv[ii] * rhsj;
	    i__++;
	}
L200:
	;
    }
/*        ------------------------------------------------- */
/*        backward substitution ...                         */
/*        ------------------------------------------------- */
    j = neqns;
    i__1 = neqns;
    for (jj = 1; jj <= i__1; jj++) {

#ifdef DEBUG
  printf("\n inside second loop..");
#endif
    
	ss = s->ve[j-1];
	istrt = xrnzv[j];
	istop = xrnzv[j + 1] - 1;
	if (istop < istrt) {
	    goto L400;
	}
	i__ = xnzsubv[j];
	i__2 = istop;
	for (ii = istrt; ii <= i__2; ii++) {
	    isub = nzsubv[i__];
	    ss -= lnzv[ii] * s->ve[isub-1];
	    i__++;
	}
L400:
	s->ve[j-1] = ss / ddiag[j];
	j--;
    }
     
    return(s);
} /* gsslv_ */

