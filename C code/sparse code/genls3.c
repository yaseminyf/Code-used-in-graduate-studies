/* last revised 14.09.1998 */
/* genls3.c called by lsqslvsn3.c */

#include "data_struct.h" 

/* extern declarations */

/*    extern int nsubs;  
    extern double Rhs;  
    extern int *SUBS, NSUBS;  
    extern double *VALUES,;  
    extern MATRIX *A;          */
        
/* ******     genls3 ..... general least squares solve     ******* */
/*     purpose - this routine performs the third step of the       */
/*        general least squares solve.                             */
/*     input parameters -                                          */
/*        input  - unit number for input file.                     */
/*        mrows  - total number of equations and constraints.      */
/*        ncols  - number of variables.                            */
/*        nseqns - number of sparse equations.                     */
/*        ndeqns - number of dense equations.                      */
/*        (xnzsub,nzsub,xrnz) - nonzero structure of the upper     */
/*                 triangular matrix.                              */
/*        invp   - inverse of column permutation.                  */
/*        TOL    - tolerance for testing zero.                     */
/*        TYPTOL - type of tolerance.                              */
/*        rkdef1 - rank deficiency of upper triangular matrix.     */
/*     output parameters -                                         */
/*        (diag,rnz,qb) - upper triangular matrix.                 */
/*        dseqns - dense equations.  (row dimension is ncols).     */
/*        dsbeqn - right hand side of dense equations.             */
/*        rowmsk - mask vector for the rows in the upper           */
/*                 triangular matrix.                              */
/*        xzrows - pointers to empty rows in the upper triangular  */
/*                 matrix.                                         */
/*        rkdef2 - rank deficiency of the upper triangular         */
/*                 matrix after the reduction of sparse            */
/*                 equations.                                      */
/*     work space -                                                */
/*        (subs,values) - space for reading rows.                  */
/*        (nusubs,nuvals) - space for row reduction.               */
/*     program subroutines -                                       */
/*        rdeqns, rkchk2, rwprep, storow.                          */
/* *************************************************************** */

void genls3(ncols,nseqns,ndeqns,ddiag,xnzsub,nzsub,xrnz,rnzv,qb,dseqns,dsbeqn,rowmsk,invp,xzrows,TOL,TYPTOL,rkdef1,rkdef2,nusubs,nuvals,T)
int ncols, nseqns, ndeqns;
double *ddiag;
int xnzsub, nzsub, xrnz;
double *rnzv, *qb, *dseqns, *dsbeqn;
int *rowmsk, invp, *xzrows;
double TOL;
int TYPTOL, rkdef1, rkdef2;
int *nusubs;
double *nuvals;
int *T;
{
    /* System generated locals */
    int dseqns_dim1, dseqns_offset, i__1, i__2;

    /* Local variables */
    static int type__, k, m,j,count;
    static double weight;
    static int kde;
    

/* *************************************************************** */

/*        ------------------ */
/*        initialization ... */
/*        ------------------ */

#ifdef DEBUG    
    fprintf(fp1,"\n inside genls3...");
#endif
    
    dseqns_dim1 = ncols;
    dseqns_offset = dseqns_dim1 + 1;
    dseqns -= dseqns_offset;

    /* Function Body */
    
    if (nseqns + ndeqns == 0) {
	goto L400;
    }

/*        ----------------------------------- */
/*        step 3: reduce sparse equations and */
/*                store dense equations ...   */
/*        ----------------------------------- */

    kde = 0;
    for ( j=1;j<=A->m;j++ ) {
/*            ----------------------------    */
/*            process sparse equations and    */
/*            store dense equations ...       */
/*            ----------------------------    */
       if ( A->rowptr[j] != 0 ) { 
          count=1;
          while ( A->rowptr[j+count] == 0 )
             count++;
          NSUBS=A->rowptr[j+count]-A->rowptr[j];
          SUBS=CREATE(NSUBS,int);
          VALUES=CREATE(NSUBS,double);
          for ( k=0;k<NSUBS;k++ ) {
             *(SUBS+k)=A->colind[A->rowptr[j]+k];
             VALUES[k]=A->Valuerow[A->rowptr[j]+k];
          }
          --SUBS;
          --VALUES;
          Rhs=1.0;
          weight=1.0;
          type__=1;
          
#ifdef DEBUG              
          fprintf(fp1,"\n SUBS and VALUES as read...\n");
          for (i__1=1;i__1<=NSUBS;i__1++) {
            fprintf(fp1,"\n %d %16.14e ",SUBS[i__1],VALUES[i__1]);
          }
#endif
          
/*            --------------------------------------------------- */
/*            skip the constraints, but process the equations ... */
/*            --------------------------------------------------- */
	   switch ((int)type__) {
	       case 1:  goto L100;
	       case 2:  goto L200;
	       case 3:  goto L300;
	       case 4:  goto L300;
	   }

L100:
/*                ----------------------- */
/*                apply weight to row ... */
/*                ----------------------- */
/* not needed since weight=1 */
/*	aplywt(nsubs, subs, values, Rhs, weight); */

/*                ------------------------------------ */
/*                current row is a sparse equation ... */
/*                ------------------------------------ */
        
#ifdef DEBUG    
        fprintf(fp1,"\n before rwprep...");
#endif
        
	rwprep(NSUBS, SUBS, VALUES, invp, nusubs, nuvals, T);
	
#ifdef DEBUG    	
	fprintf(fp1,"\n after rwprep...");
	fprintf(fp1,"\n nusubs after rwprep...\n");
        for (i__1=1; i__1<=NSUBS;i__1++) {
          fprintf(fp1,"%d ",nusubs[i__1]);
          if ((i__1%10)==0) fprintf(fp1,"\n");
        }
        fprintf(fp1,"\n nuvals after rwprep...\n");
        for (i__1=1; i__1<=ncols;i__1++) {
          fprintf(fp1,"%16.14e ",nuvals[i__1]);
          if ((i__1%5)==0) fprintf(fp1,"\n");
        } 
        fprintf(fp1,"\n Rhs: %16.14e ",Rhs);
	fprintf(fp1,"\n before rdeqns...");
#endif	
	
	rdeqns(NSUBS,nusubs,nuvals,Rhs,qb,ddiag,xnzsub,nzsub,xrnz,rnzv,rowmsk,T);
	
#ifdef DEBUG    	
	fprintf(fp1,"\n after rdeqns...");
	fprintf(fp1,"\n nuvals after rdeqns...\n");
        for (i__1=1; i__1<=ncols;i__1++) {
          fprintf(fp1,"%16.14e ",nuvals[i__1]);
          if ((i__1%5)==0) fprintf(fp1,"\n");
        }
	fprintf(fp1,"\n Rhs: %16.14e ",Rhs);
	fprintf(fp1,"\n rnzv after rdeqns...\n");
        for ( i__1=1;i__1<=spbcon_.NOFNZ;i__1++ ) {
          fprintf(fp1,"%16.14e ",rnzv[i__1]);
          if ( ((i__1)%5)==0 )
            fprintf(fp1,"\n");
        } 
        fprintf(fp1,"\n ddiag after rdeqns...\n");
        for ( i__1=1;i__1<=spbcon_.NCOLS;i__1++ ) {
          fprintf(fp1,"%16.14e ",ddiag[i__1]);
          if ( (i__1%5)==0 )
            fprintf(fp1,"\n");
        } 
#endif
        
        goto L300;

L200:
/*                ----------------------- */
/*                apply weight to row ... */
/*                ----------------------- */
/* not needed since weight=1 */
/*	aplywt(nsubs, subs, values, Rhs, weight); */

/*                ----------------------------------- */
/*                current row is a dense equation ... */
/*                ----------------------------------- */

/* this is commented out for the time-being since there are no */
/* dense equations */
/*	++kde; */
/*	storow(ncols, kde, dseqns[dseqns_offset], dsbeqn[1], nsubs,  */
/*		subs, values, Rhs, invpv); */

L300:
	;
	}
    }
/*        -------------------------------------------------------- */
/*        if there are sparse equations, determine the rank        */
/*        deficiency of the upper triangular matrix ...            */
/*        rowmsk marks the rows of the upper triangular matrix.    */
/*        rowmsk(i) is 0 if the diagonal element of row i is zero. */
/*        rowmsk(i) is 1 if row i comes from the constraints,      */
/*        rowmsk(i) is 2 if row i comes from the equations.        */
/*        -------------------------------------------------------- */
    if (nseqns == 0) {
	goto L400;
    }
    
#ifdef DEBUG        
    fprintf(fp1,"\n before rkchk2...");
#endif
    
    rkchk2(ncols,nusubs,nuvals,qb,ddiag,xnzsub,nzsub,xrnz,rnzv,rowmsk,xzrows,TOL,TYPTOL,rkdef2,T);
    
#ifdef DEBUG        
    fprintf(fp1,"\n after rkchk2...");
    fprintf(fp1,"\n rnzv after rkchk2...\n");
    for ( i__1=1;i__1<=spbcon_.NOFNZ;i__1++ ) {
      fprintf(fp1,"%16.14e ",rnzv[i__1]);
      if ( (i__1%5)==0 )
        fprintf(fp1,"\n");
    } 
    fprintf(fp1,"\n ddiag after rkchk2...\n");
    for ( i__1=1;i__1<=spbcon_.NCOLS;i__1++ ) {
      fprintf(fp1,"%16.14e ",ddiag[i__1]);
      if ( (i__1%5)==0 )
        fprintf(fp1,"\n");
    } 
#endif
    
    return;

L400:
/*        -------------------------------------------------------- */
/*        no sparse equations ... set the diagonal elements to one */
/*        if they are zero and set the corresponding elements of   */
/*        rowmsk to zero, and also set initialize rkdef2.          */
/*        -------------------------------------------------------- */
    rkdef2 = rkdef1;
    i__1 = ncols;
    for (k = 1; k <= i__1; k++) {
	if (rowmsk[k] == 2) {
	    goto L500;
	}
	ddiag[k] = 1.;
	rowmsk[k] = 0;
L500:
	;
    }
    return;

/*        -------------- */
/*        end of step 3. */
/*        -------------- */

} /* genls3 */
