/* last revision 28.06.1998 */
/* orcolb.c called by main program */

#include "data_struct.h" 

/* extern declarations */


/*    extern int xadj; 
    extern int delta; 
    extern void ipendb(), genmmd(); 
    extern int adjncy; 
    extern void smbfct(); 
    extern void copysi(); 
    extern int nudeg, nuladj, dhead, qsize; 
    extern int marker, llist, newadj; 
    extern int q, mask, clqtst, nzsub, xrnz, xnzsub; 
    extern int flag__; 
    extern int invp, perm; 
    extern int nofsub,maxlnz,maxsub; 
    extern int  *//* mdeg, ehead, tag, mdnode, */  /* limit; 
    extern int iladj, deg, ladj; */
        
/* **********     orcolb ..... find column ordering     ********** */
/*     purpose - this is the driver routine for finding a good     */
/*        column ordering.  the minimum degree algorithm is used.  */
/*        quotient graphs are used in the implementation.  the     */
/*        technique of multiple elimination and idea of external   */
/*        minimum degree are used to speed up execution.           */
/*     updated parameter -                                         */
/*        t      - the storage array with size maxsb.              */
/*     program subroutines                                         */
/*        copysi, genmmd, ipendb, smbfct                           */
/* *************************************************************** */

void orcolb(T)
int *T;
{
    /* System generated locals */
    int i__1;

    /* Local variables */

    static int ne,i__,nsub;
    
/* *************************************************************** */

/*            ------------------------------------------------- */
/*            first time orcolb is called after row input. */
/*            create adjaceny structure for (a-tranpose)(a) ... */
/*            ------------------------------------------------- */

    /* Function Body */
    
#ifdef DEBUG        
    fprintf(fp1,"\n inside orcolb.."); 
#endif
    
    spbusr_.IERRB = 0;
    
#ifdef DEBUG        
    fprintf(fp1,"\n right before ipendb.."); 
#endif

    ipendb(T);

#ifdef DEBUG        
    fprintf(fp1,"\n right after ipendb.."); 
#endif
    
    if (spbusr_.IERRB > 0) {
	return;
    }

   
/* the next statement is commented out by me */
/*    roffs = spksys_1.ratios - (float).01; */

/*        ------------------------------------ */
/*        establish pointers to data structure */
/*        vectors for the ordering. */
/*        ------------------------------------ */
    xadj = 1;
    adjncy = spbcon_.NLP1 + 1;

    spbcon_.METHOD = 54;
    
/* the next statement is commented out by me */    
/*    ne = (int) (((float) spbcon_1.nedges + roffs) / spksys_1.ratios); */

/* the next statement is introduced by me */
    ne = spbcon_.NEDGES;
    
    spbcon_.MXREQD = spbcon_.MXUSED + (spbcon_.NSHORT << 3) + (ne << 1); 
/*    spbdta_1.colstr = (float) spbcon_1.mxreqd; */

    spbmap_.PERM = spbcon_.MXUSED + 1; 
    spbmap_.INVP = spbmap_.PERM + spbcon_.NSHORT; 

#ifdef DEBUG    
    fprintf(fp1,"\n spbmap_.PERM: %d, spbmap_.INVP: %d, spbcon_.MXUSED: %d ",spbmap_.PERM,spbmap_.INVP,spbcon_.MXUSED);
#endif    
    
/*        ---------------------------------------------------- */
/*        establish pointers to work vectors for the ordering. */
/*        ---------------------------------------------------- */

    dhead = spbmap_.INVP + spbcon_.NSHORT; 
    qsize = dhead + spbcon_.NSHORT; 
    llist = qsize + spbcon_.NSHORT; 
    marker = llist + spbcon_.NSHORT; 
    newadj = marker + spbcon_.NSHORT; 

#ifdef DEBUG       
    fprintf(fp1,"\n dhead: %d, qsize: %d, llist: %d, marker: %d, newadj: %d ",dhead,qsize,llist,marker,newadj);
#endif    
          
/* i__1 = 2*nedges */

    i__1 = spbcon_.NEDGES << 1;
        
    adjncy--;
    newadj--;
    
#ifdef DEBUG        
    fprintf(fp1,"\n before copysi...");
#endif

    copysi(i__1, T, adjncy, newadj);
    
#ifdef DEBUG        
    fprintf(fp1,"\n after copysi...");
#endif    
        
    delta = 0;
/*        ------------------------------------------ */
/*        order the given matrix by the md ordering. */
/*        ------------------------------------------ */
    
    invp=spbmap_.INVP-1;
    perm=spbmap_.PERM-1;
    dhead--;
    qsize--;
    llist--;
    marker--;
    xadj--;

    spbcon_.NOFSUB=0;
    
#ifdef DEBUG    
    fprintf(fp1,"\n before genmmd...");
#endif

    genmmd(spbcon_.NCOLS,xadj,newadj,invp,perm,delta,dhead,qsize,llist,marker,spksys_.MAXINT,spbcon_.NOFSUB,T);

#ifdef DEBUG    
    fprintf(fp1,"\n after genmmd...");
#endif
    
    spbcon_.MXUSED = spbmap_.INVP + spbcon_.NSHORT - 1; 

#ifdef DEBUG    
    fprintf(fp1,"\n spbcon_.MXUSED: %d ",spbcon_.MXUSED);
#endif    

/*        ------------------------------------------------ */
/*        establish pointers to data structure vectors for */
/*        storage allocation.                              */
/*        ------------------------------------------------ */

/* the following statement is commented out by me */
/*    rnsub = (float) spbcon_.nofsub; */
    nsub = spbcon_.NOFSUB; 
    spbcon_.MXREQD=spbcon_.MXUSED+ 3*spbcon_.NSHORT +(spbcon_.NLP1 << 1)+nsub; 

    spbmap_.XRNZ = spbmap_.INVP + spbcon_.NSHORT; 
    spbmap_.XNZSUB = spbmap_.XRNZ + spbcon_.NLP1; 
    spbmap_.NZSUB = spbmap_.XNZSUB + spbcon_.NLP1; 

#ifdef DEBUG    
    fprintf(fp1,"\n nofsub: %d, spbcon_.MXREQD: %d ",nsub,spbcon_.MXREQD);
    fprintf(fp1,"\n spbmap_.XRNZ: %d, spbmap_.XNZSUB: %d, spbmap_.NZSUB: %d ",spbmap_.XRNZ,spbmap_.XNZSUB,spbmap_.NZSUB);
#endif    
    
/*        -------------------------------------------------- */
/*        establish pointers to work vectors for allocation. */
/*        -------------------------------------------------- */

    q = spbmap_.NZSUB + nsub; 
    mask = q + spbcon_.NSHORT; 
    clqtst = mask + spbcon_.NSHORT; 

#ifdef DEBUG    
    fprintf(fp1,"\n q: %d, mask: %d, clqtst: %d ",q,mask,clqtst);
#endif    

/*        ------------------------------------------------- */
/*        perform the storage allocation by calling smbfct. */
/*        ------------------------------------------------- */
    
    q--;
    mask--;
    clqtst--;
    
    nzsub=spbmap_.NZSUB-1;
    xnzsub=spbmap_.XNZSUB-1;
    xrnz=spbmap_.XRNZ-1;

/* the next statements are here just for the sake of testing! */

/*    for ( i__=1;i__<=spbcon_.NOFSUB;i__++ ) */
/*      T[nzsub+i__]=0; */
/*    for ( i__=1;i__<=spbcon_.NCOLS;i__++)     { */
/*      T[q+i__]=0; */
/*      T[xrnz+i__]=0; */
/*      T[xnzsub+i__]=0; */
/*    } */
    
#ifdef DEBUG    
    fprintf(fp1,"\n before smbfct...");
#endif

    smbfct(spbcon_.NCOLS,xadj,adjncy,perm,invp,xrnz,spbcon_.NOFNZ,xnzsub,nzsub,spbcon_.NOFSUB,q,mask,clqtst,flag__,T);

#ifdef DEBUG    
    fprintf(fp1,"\n after smbfct...");
    fprintf(fp1,"\n spbcon_.NOFSUB: %d, spbcon_.NOFNZ: %d, flag__: %d ",spbcon_.NOFSUB,spbcon_.NOFNZ,flag__);
    fprintf(fp1,"\n xrnz after smbfct...\n");
    i__1=spbcon_.NCOLS+1;
    for (i__=1;i__<=i__1;i__++) {
      fprintf(fp1,"%d ",T[xrnz+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n xnzsub after smbfct...\n");
    i__1=spbcon_.NCOLS+1;
    for (i__=1;i__<=i__1;i__++) {
      fprintf(fp1,"%d ",T[xnzsub+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n nzsub after smbfct...\n");
    i__1=spbcon_.NOFSUB;
    for (i__=1;i__<=i__1;i__++) {
      fprintf(fp1,"%d ",T[nzsub+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
#endif
    
/*        ------------------------------------------------- */
/*        if data structure allocation has been successful, */
/*        reorganize the storage vector.                    */
/*        ------------------------------------------------- */

/* the following statements are commented out by me */
/*    rnsub = (float) spbcon_1.nofsub; */
/*    nsub = (int) ((rnsub + roffs) / spksys_1.ratios); */ 
/*    spbdta_1.alostr = (float) (spbcon_1.mxused + spbcon_1.nshort * 3); */

    nsub=spbcon_.NOFSUB;
    spbcon_.MXUSED = spbmap_.NZSUB + nsub - 1;
    
#ifdef DEBUG    
    fprintf(fp1,"\n nsub: %d, spbcon_.MXUSED: %d ",nsub,spbcon_.MXUSED);
#endif    
    
/* the following subroutine is not needed anymore */
/*    reorgb(T[1]); */

/*        ----------------------------------------------------- */
/*        establish pointers to storage vectors for the matrix. */
/*        ----------------------------------------------------- */

/* the following statements are left as they are for the sake of */
/* referances to dense equations and constraints in the coming */
/* subroutines */

    spbmap_.RHS = 1; 
    spbmap_.DIAG = spbcon_.MXUSED + 1; 
    spbmap_.RNZ = spbmap_.DIAG + spbcon_.NCOLS; 
    spbmap_.DSEQNS = spbmap_.RNZ + spbcon_.NOFNZ; 
    spbmap_.DSBEQN = spbmap_.DSEQNS + spbcon_.NCOLS * spbcon_.NDEQNS; 
    spbmap_.DSCONS = spbmap_.DSBEQN + spbcon_.NDEQNS; 
    spbmap_.DSBCON = spbmap_.DSCONS + spbcon_.NCOLS * spbcon_.NDCONS; 
    spbmap_.ROWMSK = spbmap_.DSBCON + spbcon_.NDCONS; 
    spbmap_.XZROWS = spbmap_.ROWMSK + spbcon_.NSHORT; 

/* the following statements are commented out by me */    
/*    spbdta_1.slvstr = (float) (spbmap_1.xzrows + spbcon_1.nshort - 1); */
/*    spbdta_1.overhd = (float) (spbcon_1.ncols + spbcon_1.nshort * 5 + ( */
/*	    spbcon_1.nlp1 << 1) + nsub); */
/*    spbcon_1.mxreqd = (int) spbdta_1.slvstr; */

    return;

} /* orcolb */

