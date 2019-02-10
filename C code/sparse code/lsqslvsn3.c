/* last revision 14.09.1998 */
/* lsqslvsn3.c called by the main program */

#include "data_struct.h" 
   
/* **********     lsqslv ..... least squares solve     *********** */
/*     purpose - this is a driver for solving a sparse linear      */
/*        least squares problem using givens rotations.  dense     */
/*        equations are allowed.  the algorithm can also handle    */
/*        both sparse and dense constraints. the constraints need  */
/*        not be consistent, and neither the constraints nor the   */
/*        equations need be of full rank.                          */
/*     input parameters -                                          */
/*        TOL    - user specified tolerance for determining null   */
/*                 rows.                                           */
/*        TYPTOL - type of tolerance:                              */
/*                 0 - absolute tolerance.                         */
/*                 1 - relative tolerance (relative to largest     */
/*                     diagonal element).                          */
/*     updated parameter -                                         */
/*        t      - the storage array with size maxsb.              */
/*                 on output, the first ncols locations of t       */
/*                 contain the solution.                           */
/*     program subroutines -                                       */
/*        genls3, zerols, zerorv.                                  */
/* *************************************************************** */

void lsqslvsn3(TOL,TYPTOL,T, xnzsub, nzsub, xrnz, invp, perm) 
double TOL;
int TYPTOL;
int *T;
int xnzsub, nzsub, xrnz, invp, perm;
{
    /* System generated locals */
    int i__1, i__;

    /* Local variables */
    
    static int neqs;
    static int nhc, nhr, k;


/* *************************************************************** */

    /* Function Body */

/* the following are commented out by me */
/*    roffs = spksys_1.ratios - (float).01; */
/*    roffl = spksys_1.ratiol - (float).01; */

/*    len = (int) (((float) spbcon_1.nzmax + roffl) / spksys_1.ratiol); */
/*    subs = spbmap_1.xzrows + spbcon_1.nshort; */
/*    values = subs + len; */
/*    nusubs = values + spbcon_1.nzmax; */
/*    nuvals = nusubs + spbcon_1.nshort; */
/*    spbcon_1.mxreqd = nuvals + spbcon_1.ncols - 1; */

/*        ------------------ */
/*        initialization ... */
/*        ------------------ */

/* create the necessary storage for the coming subroutines */

#ifdef DEBUG    
    fprintf(fp1,"\n inside lsqslvsn3...");
    fprintf(fp1,"\n TOL: %16.14e, TYPTOL: %d ",TOL,TYPTOL);
#endif
    
    rhsv=CREATE(spbcon_.NSEQNS,double);
    --rhsv;

    
#ifdef DEBUG        
    fprintf(fp1,"\n before zerorv...");
#endif
    
    zerorv(spbcon_.NSEQNS, rhsv);

#ifdef DEBUG        
    fprintf(fp1,"\n after zerorv...");
    
    fprintf(fp1,"\n rhsv initial condition...\n");
    for ( i__=1;i__<=spbcon_.NSEQNS;i__++ ) {
      fprintf(fp1,"%16.14e ",rhsv[i__]);
      if ( ((i__)%5)==0 )
        fprintf(fp1,"\n");
    } 
    fprintf(fp1,"\n before zerorv...");
#endif
    
/* the following are commented out because there are no */
/* dense equations or constraints */    
/*    i__1 = spbcon_.ncols * spbcon_.ndeqns; */
/*    zerorv(i__1, T[spbmap_1.dseqns]); */
/*    zerorv(spbcon_1.ndeqns, T[spbmap_1.dsbeqn]); */
/*    i__1 = spbcon_1.ncols * spbcon_1.ndcons; */
/*    zerorv(i__1, T[spbmap_1.dscons]); */
/*    zerorv(spbcon_1.ndcons, T[spbmap_1.dsbcon]); */

    neqs = spbcon_.NSEQNS + spbcon_.NDEQNS + spbcon_.NSCONS + spbcon_.NDCONS;

#ifdef DEBUG    
    fprintf(fp1,"\n neqs = %d",neqs);
#endif
    
    /* create the necessary storage */
    nuvals=CREATE(spbcon_.NCOLS,double); 
    --nuvals;	  

#ifdef DEBUG           
    fprintf(fp1,"\n before zerorv...");
#endif    

    zerorv(spbcon_.NCOLS, nuvals); 
    
#ifdef DEBUG        
    fprintf(fp1,"\n after zerorv...");

    fprintf(fp1,"\n nuvals initial condition...\n");
    for ( i__=1;i__<=spbcon_.NCOLS;i__++ ) {
      fprintf(fp1,"%16.14e ",nuvals[i__]);
      if ( ((i__)%5)==0 )
        fprintf(fp1,"\n");
    } 
#endif
        
/* the following subroutine is not needed. the only needed part */
/* is copied below. */
/*    genls1(neqs, spbcon_.ncols, spbcon_.nscons, spbcon_.ndcons, ddiag, // */
/*           xnzsubv, nzsubv, xrnzv, rnzv, rhsv, T[spbmap_.dscons], // */
/*           T[spbmap_.dsbcon], rowmsk, invpv, T[spbmap_.xzrows], TOL, // */
/*	   TYPTOL, spbcon_.rkdef1, subs, values, nusubs, nuvals); */
/*        ------------------------------------------------- */
/*        no sparse constraints ... set rowmsk to all ones, */
/*        and also set rank deficiency. */
/*        ------------------------------------------------- */

    rowmsk=CREATE(spbcon_.NCOLS,int);
    --rowmsk;
    spbcon_.RKDEF1 = spbcon_.NCOLS;
    i__1 = spbcon_.NCOLS;
    for (k = 1; k <= i__1; k++) {
	rowmsk[k] = 1;
    }	

#ifdef DEBUG    
    fprintf(fp1,"\n rowmsk initial condition...\n");
    for ( i__=1;i__<=spbcon_.NCOLS;i__++ ) {
      fprintf(fp1,"%d ",rowmsk[i__]);
      if ( ((i__)%5)==0 )
        fprintf(fp1,"\n");
    } 
    fprintf(fp1,"\n spbcon_.RKDEF1 = %d",spbcon_.RKDEF1);
#endif
    
/* the following is commented out by me */	    
/*    spbdta_1.slvstr = (float) spbcon_1.mxreqd; */
    nhc = spbcon_.NDCONS;
    nhr = spbcon_.NCOLS - spbcon_.RKDEF1;
    if (nhc < nhr) {
	goto L300;
    }
    nhc = nhr;
    nhr = spbcon_.NDCONS;
L300:
/* the following is commented out by me */
/*    lhc = max(nhc,1); */
/*    h__ = spbmap_.XZROWS + spbcon_.NSHORT; */
/*    g = h__ + lhr * nhc; */
/*    u = g + spbcon_.NDCONS; */
/*    v = u + lhr * nhc; */
/*    rscon = v + lhc * nhc; */
/*    sv = rscon + (spbcon_.NCOLS - spbcon_.RKDEF1); */
/*    work = sv + nhc; ( */
/*    spbcon_.MXREQD = work + nhc - 1; */
/*    i__1 = spbcon_.NCOLS - spbcon_.RKDEF1; */

/* the following is commented out by me */
/*    zerorv(i__1, T[rscon]); */
/* the following is not needed */
/*    genls2(spbcon_1.ncols, spbcon_1.nscons, spbcon_1.ndcons, T[ */
/*	    spbmap_1.diag], T[spbmap_1.xnzsub], T[spbmap_1.nzsub], T[ */
/*	    spbmap_1.xrnz], T[spbmap_1.rnz], T[spbmap_1.rhs], T[ */
/*	    spbmap_1.dscons], T[spbmap_1.dsbcon], T[spbmap_1.rowmsk], T[ */
/*	    spbmap_1.xzrows], lhr, lhc, T[h__], T[g], T[u], T[v], T[ */
/*	    rscon], T[sv], T[work], TOL, TYPTOL, spbcon_1.rkdef1,  */
/*	    spbcon_1.rkg, iflag); */

/* the following is commented out by me */
/*    if (spbdta_1.slvstr < (real) spbcon_1.mxreqd) { */
/*	spbdta_1.slvstr = (real) spbcon_1.mxreqd; */
/*    } */

/* the following is commented out by me */
/*    len = (int) (((float) spbcon_1.nzmax + roffl) / spksys_1.ratiol); */
/*    subs = spbmap_1.xzrows + spbcon_1.nshort; */
/*    values = subs + len; */
/*    nusubs = values + spbcon_1.nzmax; */
/*    nuvals = nusubs + spbcon_1.nshort; */
/*    spbcon_1.mxreqd = nuvals + spbcon_1.ncols - 1; */

    neqs = spbcon_.NSEQNS + spbcon_.NDEQNS + spbcon_.NSCONS + spbcon_.NDCONS;
    
    xzrows=CREATE(spbcon_.NSHORT,int);
    --xzrows;
    
    nusubs=CREATE(spbcon_.NSHORT,int); 
    --nusubs; 
    	    
/* this statment is not needed since no change has been done */
/* in nuvals since it was created */	    
/*    zerorv(spbcon_.NCOLS, nuvals); */

#ifdef DEBUG    
    fprintf(fp1,"\n before genls3...");
#endif
    
    genls3(spbcon_.NCOLS,spbcon_.NSEQNS,spbcon_.NDEQNS,ddiag,xnzsub,nzsub,xrnz,rnzv,rhsv,T[spbmap_.DSEQNS-1],T[spbmap_.DSBEQN-1],rowmsk,invp,xzrows,TOL,TYPTOL,spbcon_.RKDEF1,spbcon_.RKDEF2,nusubs,nuvals,T);
    
#ifdef DEBUG        
    fprintf(fp1,"\n after genls3...");
#endif
    
/* the following is commented out by me */
/*    if (spbdta_1.slvstr < (float) spbcon_1.mxreqd) { */
/*	spbdta_1.slvstr = (float) spbcon_1.mxreqd; */
/*    } */
/* the following is not needed */
/*    genls4(spbcon_1.ncols, spbcon_1.nscons, spbcon_1.ndcons,  */
/*	    spbcon_1.nseqns, spbcon_1.ndeqns, spbcon_1.rkg, T[ */
/*	    spbmap_1.diag], T[spbmap_1.xnzsub], T[spbmap_1.nzsub], T[ */
/*	    spbmap_1.xrnz], T[spbmap_1.rnz], T[spbmap_1.rhs], T[ */
/*	    spbmap_1.dscons], T[spbmap_1.dsbcon], T[spbmap_1.dseqns], T[ */
/*	    spbmap_1.dsbeqn], T[spbmap_1.rowmsk], spbcon_1.rkt1,  */
/*	    spbcon_1.rkt2, TOL, TYPTOL); */
/*    i__1 = spbcon_1.ncols; */
/*    for (i__ = 1; i__ <= i__1; ++i__) { */
/*	T[i__] = T[spbmap_1.diag + i__ - 1]; */
/*	ddiag[i__] = T[spbmap_1.diag + i__ - 1]; */
/*    } */
/*    copyrest(spbcon_1.ncols, spbcon_1.nofnz, T[spbmap_1.xrnz], T[ */
/*	    spbmap_1.xnzsub], T[spbmap_1.perm], T[spbmap_1.invp], T[ */
/*	    spbmap_1.nzsub], xrnzv[1], xnzsubv[1], permv[1], invpv[1],  */
/*	    nzsubv[1]); */
/*    i__1 = spbcon_1.nofnz; */
/*    for (i__ = 1; i__ <= i__1; ++i__) { */
/*	rnzv[i__] = T[spbmap_1.rnz + i__ - 1]; */
/*    } */

    i__1 = spbcon_.NCOLS;
    for (i__ = 1; i__ <= i__1; i__++) {
	xrnzv[i__] = T[xrnz+i__];
	xnzsubv[i__] = T[xnzsub+i__];
	permv[i__] = T[perm+i__];
	invpv[i__] = T[invp+i__];
    }
    xrnzv[i__1+1]=T[xrnz+i__+1];
    xnzsubv[i__+1]=T[xnzsub+i__+1];
    i__1 = spbcon_.NOFNZ;
    for (i__ = 1; i__ <= i__1; i__++) {
	nzsubv[i__] = T[nzsub+i__];
    }

#ifdef DEBUG        
    fprintf(fp1,"\n xrnzv after lsqslvsn3...\n");
    i__1=spbcon_.NCOLS+1;
    for ( i__=1;i__<=i__1;i__++ ) {
      fprintf(fp1,"%d ",xrnzv[i__]);
      if ( ((i__)%5)==0 )
        fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n xnzsubv after lsqslvsn3...\n");
    i__1=spbcon_.NCOLS+1;
    for ( i__=1;i__<=i__1;i__++ ) {
      fprintf(fp1,"%d ",xnzsubv[i__]);
      if ( ((i__)%5)==0 )
        fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n permv after lsqslvsn3...\n");
    i__1=spbcon_.NCOLS;
    for ( i__=1;i__<=i__1;i__++ ) {
      fprintf(fp1,"%d ",permv[i__]);
      if ( ((i__)%5)==0 )
        fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n invpv after lsqslvsn3...\n");
    i__1=spbcon_.NCOLS;
    for ( i__=1;i__<=i__1;i__++ ) {
      fprintf(fp1,"%d ",invpv[i__]);
      if ( ((i__)%5)==0 )
        fprintf(fp1,"\n");
    }
    i__1 = spbcon_.NOFNZ;
    fprintf(fp1,"\n nzsubv after lsqslvsn3...\n");
    for ( i__=1;i__<=i__1;i__++ ) {
      fprintf(fp1,"%d ",nzsubv[i__]);
      if ( ((i__)%5)==0 )
        fprintf(fp1,"\n");
    }
#endif
    
/* garbage collection */

    ++nuvals;
    ++rhsv;
    ++xzrows;
    ++nusubs;
    ++rowmsk;
    free(nuvals);
    free(rhsv);
    free(xzrows);
    free(nusubs);
    free(rowmsk);
       
    return;
    
} /* lsqslvsn3 */
