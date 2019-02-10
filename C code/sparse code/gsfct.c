/* last updated: 01.09.1998 */
/* gsfct called by main program */

/* #include "data_struct.h" */

/* ******     gsfct ..... general sparse symmetric factorization ******* */
/*     purpose - this routine performs the symmetric                     */
/*        factorization for a general sparse system, stored in           */
/*        the compressed subscript data format                           */
/*     input parameters -                                                */
/*        neqns - number of equations                                    */
/*        xrnzv - index vector for lnzv. xrnzv(i) points to the          */
/*                start of nonzeros in column i of factor L.             */
/*        (xnzsubv,nzsubv) - the compressed subscript data               */
/*                         structure for factor L.                       */
/*     updated parameters -                                              */
/*        lnzv - on input contains nonzeros of A, and on return          */
/*              the nonzeros of L.                                       */
/*        ddiag - the diagonal of L overwrites that of A.                */
/*        flag__ - the error flag. It is set to 1 if a zero or           */
/*                negative square root occurs during the factorization.  */
/*     working parameters -                                              */
/*        llink - at step j, the list in                                 */
/*               llink(j), llink(llink(j)), ....                         */
/*               consists of those columns that will modify              */
/*               the column L(*,j).                                      */
/*        first - temporary vector to point to the first                 */
/*                nonzero in each column that will be used               */
/*                next for modification.                                 */
/*        temp - a temporary vector to accumulate modifications.         */
/* ***************************************************************       */

void gsfct(neqns,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,llink,first,temp,flag__)
int neqns, *xrnzv;
double *lnzv;
int *xnzsubv, *nzsubv;
double *ddiag;
int *llink, *first;
double *temp;
int flag__;
{
    /* System generated locals */
    int i__1, i__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static int isub, newk, i__, j, k;
    static double diagj;
    static int istop, istrt, ii, kfirst;
    static double ljk;

/* *************************************************************** */

/*        ------------------             */
/*        initialize working vectors ... */
/*        ------------------             */


    /* Function Body */

#ifdef DEBUG
   fprintf(fp1,"\n inside gsfct..");  
#endif
     
    i__1 = neqns;
    for (i__ = 1; i__ <= i__1; i__++) {
	llink[i__] = 0;
	temp[i__] = 0.;
    }
#ifdef DEBUG
   fprintf(fp1,"\n initialized llink and temp..");  
#endif

/*        ------------------------------------------------- */
/*        compute column L(*,j) for j=1,.., neqns           */
/*        ------------------------------------------------- */

    i__1 = neqns;
    for (j = 1; j <= i__1; j++) {
/*            ------------------------------             */
/*            for each column L(*,j) that affects L(*,j) */
/*            ------------------------------             */
	diagj = 0.;
	newk = llink[j];
L200:
	k = newk;
	if (k == 0) {
	    goto L400;
	}
	newk = llink[k];
/*               -----------------------                 */
/*               outer product modification of L(*,j) by */
/*               L(*,k) starting at first(k) of L(*,k)   */
/*               -----------------------                 */
	kfirst = first[k];
	ljk = lnzv[kfirst];
	diagj += ljk * ljk;
	istrt = kfirst + 1;
	istop = xrnzv[k + 1] - 1;
	if (istop < istrt) {
	    goto L200;
	}
/*                  -----------------------                    */
/*                  before modification, update vectors first, */
/*                  and link for future modification steps.    */
/*                  -----------------------                    */
	first[k] = istrt;
	i__ = xnzsubv[k] + (kfirst - xrnzv[k]) + 1;
	isub = nzsubv[i__];
	llink[k] = llink[isub];
	llink[isub] = k;
/*                  -----------------------                 */
/*                  the actual mod is saved in vector temp. */
/*                  -----------------------                 */
	i__2 = istop;
	for (ii = istrt; ii <= i__2; ii++) {
	    isub = nzsubv[i__];
	    temp[isub] += lnzv[ii] * ljk;
	    i__++;
	}
	goto L200;
/*            ------------------------------                 */
/*            apply the modifications accumulated in temp to */
/*            column L(*,j)                                  */
/*            ------------------------------                 */
L400:
	diagj = ddiag[j] - diagj;
	if (diagj <= 0.) {
	    goto L700;
	}
	diagj = sqrt(diagj);
	ddiag[j] = diagj;
	istrt = xrnzv[j];
	istop = xrnzv[j + 1] - 1;
	if (istop < istrt) {
	    goto L600;
	}
	first[j] = istrt;
	i__ = xnzsubv[j];
	isub = nzsubv[i__];
	llink[j] = llink[isub];
	llink[isub] = j;
	i__2 = istop;
	for (ii = istrt; ii <= i__2; ii++) {
	    isub = nzsubv[i__];
	    lnzv[ii] = (lnzv[ii] - temp[isub]) / diagj;
	    temp[isub] = 0.;
	    i__++;
	}
L600:
	;
    }
    return;
/*        ---------------------------------------------------------- */
/*        error - zero or negative square root in factorization      */
/*        ---------------------------------------------------------- */
L700:
    flag__ = 1;
    return;
} /* gsfct_ */

