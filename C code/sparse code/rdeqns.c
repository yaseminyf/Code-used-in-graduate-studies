/* last revised 14.08.1998 */
/* rdeqns.c called by genls3.c and rkchk2.c */

#include "data_struct.h" 

/* *********     rdeqns ..... reduce sparse equation     ********* */
/*     purpose - this routine eliminates a given row using a       */
/*        sequence of elementary transformations and givens        */
/*        transformations.  the transformations are not saved.     */
/*        elementary transformations are used only when the        */
/*        pivot row comes from a constraint.                       */
/*     input parameters -                                          */
/*        NSUBS  - number of nonzeros in the given row.            */
/*        SUBS   - column subscripts of nonzeros in given row.     */
/*        (xnzsub,nzsub) - column subscripts of nonzeros in        */
/*                 the final upper trapezoidal matrix.             */
/*        rowmsk - an integer mask vector.  a pivot row from       */
/*                 the upper triangular matrix comes from a        */
/*                 constraint if the corresponding entry in this   */
/*                 vector is equal to 2.                           */
/*     updated parameters -                                        */
/*        nuvals    - row to be transformed.                       */
/*        Rhs    - corresponding entry in right hand side.         */
/*        qb     - partially-tranformed right hand side.           */
/*        diag   - off-diagonal nonzeros in the final upper        */
/*                 trapezoidal matrix.                             */
/*        (xrnz,rnz) - nonzeros of the final upper trapezoidal     */
/*                 matrix.                                         */
/* *************************************************************** */

void rdeqns(NSUBS,SUBS,nuvals,Rhs,qb,ddiag,xnzsub,nzsub,xrnz,rnzv,rowmsk,T)
int NSUBS, *SUBS;
double *nuvals, Rhs, *qb, *ddiag;
int xnzsub, nzsub, xrnz;
double *rnzv;
int *rowmsk;
int *T;
{
    /* System generated locals */
    int i__1;
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    static double beta;
    static int icol;
    static double temp, mult;
    static int i__, j, k;
    static double gamma, alpha, delta, sigma;
    static int istop, istrt, switch__, nzsubj;
    static double eta;

/* *************************************************************** */

    /* Function Body */
    
#ifdef DEBUG        
    fprintf(fp1,"\n inside rdeqns...");
#endif
    
L100:
    if (NSUBS < 1) {
	return;
    }

/*            ---------------------------- */
/*            locate the first nonzero ... */
/*            ---------------------------- */
    i__1 = NSUBS;
    for (i__ = 1; i__ <= i__1; i__++) {
	icol = SUBS[i__];
	
#ifdef DEBUG    	
	fprintf(fp1,"\n icol: %d ",icol);
	fprintf(fp1,"\n nuvals[%d]: %lf ",icol, nuvals[icol]);
#endif
	
	if (nuvals[icol] != 0.) {
	    goto L300;
	}
    }
    return;

L300:
/*            ------------------------------- */
/*            determine type of pivot row ... */
/*            ------------------------------- */
    switch__ = rowmsk[icol];
    switch ((int)switch__) {
	case 1:  goto L400;
	case 2:  goto L500;
    }

L400:
/*                ----------------------------------------- */
/*                construct a givens rotation to annihilate */
/*                the first nonzero in this row ... */
/*                ----------------------------------------- */
/* Computing MAX */
    d__3 = (d__1 = ddiag[icol], ABS(d__1)), d__4 = (d__2 = nuvals[icol], ABS(d__2));

    eta = MAX(d__3,d__4);
    alpha = ddiag[icol] / eta;
    beta = nuvals[icol] / eta;
    delta = sqrt(alpha * alpha + beta * beta);
    gamma = alpha / delta;
    sigma = beta / delta;
    ddiag[icol] = eta * delta;
    goto L600;

L500:
/*                -------------------------------------------- */
/*                construct an elementary transformation to */
/*                annihilate the first nonzero in this row ... */
/*                -------------------------------------------- */
    mult = -nuvals[icol] / ddiag[icol];

L600:
    nuvals[icol] = 0.;

/*            ------------------------------------------ */
/*            apply transformation to remaining elements */
/*            in the row ... */
/*            ------------------------------------------ */
    k = 0;
    istrt = T[xrnz+icol];
    istop = T[xrnz+icol + 1] - 1;
    if (istop < istrt) {
	goto L1100;
    }
    j = T[xnzsub+icol];
    i__1 = istop;
    for (i__ = istrt; i__ <= i__1; i__++) {
	nzsubj = T[nzsub+j];
	switch ((int)switch__) {
	    case 1:  goto L700;
	    case 2:  goto L800;
	}
L700:
	temp = gamma * rnzv[i__] + sigma * nuvals[nzsubj];
	nuvals[nzsubj] = gamma * nuvals[nzsubj] - sigma * rnzv[i__];
	rnzv[i__] = temp;
	goto L900;
L800:
	nuvals[nzsubj] += mult * rnzv[i__];
L900:
	k++;
	SUBS[k] = nzsubj;
	j++;
    }
L1100:
    NSUBS = k;

/*            ----------------------------------------------- */
/*            apply transformation to the corresponding entry */
/*            in the right hand side ... */
/*            ----------------------------------------------- */
    switch ((int)switch__) {
	case 1:  goto L1200;
	case 2:  goto L1300;
    }
L1200:
    temp = gamma * qb[icol] + sigma * Rhs;
    Rhs = gamma * Rhs - sigma * qb[icol];
    qb[icol] = temp;
    goto L100;
L1300:
    Rhs += mult * qb[icol];
    goto L100;

} /* rdeqns */

