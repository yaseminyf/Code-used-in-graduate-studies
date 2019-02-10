/* last revised 29.06.1998 */
/* rkchk2.c called by genls3.c */

#include "data_struct.h" 

/* *******     rkchk2 ..... check for rank deficiency     ******** */
/*     purpose - this routine determines the rank deficiency       */
/*        of an upper triangular matrix by counting the number     */
/*        of small diagonal elements.  the upper triangular        */
/*        matrix is resulted from the reduction of sparse          */
/*        equations.                                               */
/*     input parameters -                                          */
/*        n      - order of the upper triangular matrix.           */
/*        (xnzsub,nzsub) - column subscripts of nonzeros of the    */
/*                 upper trapezoidal matrix.                       */
/*        TOL    - user specified tolerance.                       */
/*        TYPTOL - type of tolerance:                              */
/*                 0 - absolute tolerance.                         */
/*                 1 - relative tolerance.                         */
/*        rowmsk - an integer vector for identifying the           */
/*                 constraints.                                    */
/*     updated parameters -                                        */
/*        qb     - transformed right hand side.                    */
/*        diag   - diagonal elements of upper trapezoidal matrix.  */
/*        (xrnz,rnz) - upper trapezoidal matrix.                   */
/*     output parameters -                                         */
/*        xzrows - pointers to columns in the upper triangular     */
/*                 matrix whose diagonal elements are zero.        */
/*        rkdef2 - rank deficiency of the triangular matrix.       */
/*     working arrays -                                            */
/*        (SUBS,VALUES) - for processing rows with small diagonal  */
/*                 elements.                                       */
/*     program subroutine -                                        */
/*        rdeqns.                                                  */
/* *************************************************************** */
    
void rkchk2(n,SUBS,VALUES,qb,ddiag,xnzsub,nzsub,xrnz,rnzv,rowmsk,xzrows,TOL,TYPTOL,rkdef2,T)
int n, *SUBS;
double *VALUES, *qb, *ddiag;
int xnzsub, nzsub, xrnz;
double *rnzv;
int *rowmsk, *xzrows;
double TOL;
int TYPTOL, rkdef2,*T;
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;

    /* Local variables */
    static double d__;
    static int i__, j, k;
    static int kstop, kstrt;
    static double diamax;
    static int nzsubj;
    static double eps;
    
/* *************************************************************** */

    /* Function Body */

#ifdef DEBUG    
    fprintf(fp1,"\n inside rkchk2...");
#endif    
    
    eps = TOL;
    if (TYPTOL == 0) {
	goto L200;
    }
/*            ------------------------------------ */
/*            determine the diagonal element which */
/*            has the largest magnitude ...        */
/*            ------------------------------------ */
    diamax = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; i__++) {
	d__ = (d__1 = ddiag[i__], ABS(d__1));
	if (d__ > diamax) {
	    diamax = d__;
	}
    }
    eps = diamax * TOL;
    if (eps < TOL) {
	eps = TOL;
    }

L200:
    rkdef2 = 0;

/*        ------------------------------------------ */
/*        determine the rows whose diagonal elements */
/*        are less than eps ... */
/*        ------------------------------------------ */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; i__++) {
	if ((d__1 = ddiag[i__], ABS(d__1)) > eps) {
	    goto L500;
	}
	rowmsk[i__] = 0;
	rkdef2++;
	xzrows[rkdef2] = i__;
	ddiag[i__] = 1.;
	nsubs = 0;
	kstrt = T[xrnz+i__];
	kstop = T[xrnz+i__ + 1] - 1;
	if (kstop < kstrt) {
	    goto L400;
	}
	j = T[xnzsub+i__];
	i__2 = kstop;
	for (k = kstrt; k <= i__2; k++) {
	    nzsubj = T[nzsub+j];
	    nsubs++;
	    SUBS[nsubs] = nzsubj;
	    VALUES[nzsubj] = rnzv[k];
	    rnzv[k] = 0.;
	    j++;
	}
L400:
	b = qb[i__];
	qb[i__] = 0.;
/*            --------------------------------------- */
/*            annihilate all elements in this row ... */
/*            --------------------------------------- */
	if (nsubs <= 0) {
	    goto L500;
	}
	
#ifdef DEBUG    	
	fprintf(fp1,"\n before rdeqns...");
	fprintf(fp1,"\n ddiag... \n");
        for ( i__=1;i__<=n;i__++ ) {
          fprintf(fp1,"%lf ",ddiag[i__]);
          if ( (i__%5)==0 )
            fprintf(fp1,"\n");
        } 
        fprintf(fp1,"\n nsubs = %d",nsubs);
#endif
        
	rdeqns(nsubs,SUBS,VALUES,b,qb,ddiag,xnzsub,nzsub,xrnz,rnzv,rowmsk,T);
	
#ifdef DEBUG    	
	fprintf(fp1,"\n after rdeqns...");
	fprintf(fp1,"\n rnzv after rdeqns...");
        for ( i__=1;i__<=spbcon_.NOFNZ;i__++ ) {
          fprintf(fp1,"%lf ",rnzv[i__]);
          if ( ((i__)%5)==0 )
            fprintf(fp1,"\n");
        } 
        fprintf(fp1,"\n ddiag after rdeqns...");
        for ( i__=1;i__<=spbcon_.NCOLS;i__++ ) {
          fprintf(fp1,"%lf ",ddiag[i__]);
          if ( ((i__)%5)==0 )
            fprintf(fp1,"\n");
        } 
#endif
        
L500:
	;
    }
    return;

} /* rkchk2 */

