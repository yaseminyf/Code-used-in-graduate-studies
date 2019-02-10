/* last revised 29.06.1998 */
/* smbfct.c called by orcolb.c */

#include "data_struct.h" 

/* *********     smbfct ..... symbolic factorization     ********* */
/*     purpose - this routine performs symbolic factorization      */
/*        on a permuted linear system and it also sets up the      */
/*        compressed data structure for the system.                */
/*     input parameters -                                          */
/*        neqns  - number of equations.                            */
/*        (xadj,adjncy) - the adjacency structure.                 */
/*        (perm,invp) - the permutation vector and its inverse.    */
/*     updated parameter -                                         */
/*        maxsub - size of the subscript array nzsub.  on return,  */
/*                 it contains the number of subscripts used       */
/*     output parameters -                                         */
/*        xrnz   - index into the nonzero storage vector lnz.      */
/*        (xnzsub,nzsub) - the compressed subscript vectors.       */
/*        maxlnz - the number of nonzeros found.                   */
/*        flag   - error flag.  positive value indicates that      */
/*                 nzsub array is too small.                       */
/*     working parameters -                                        */
/*        mask - a vector of size neqns.  at the kth step,         */
/*                 mask(k), mask(mask(k)) , ..., is a list         */
/*                 containing all those columns l(*,j) with j      */
/*                 less than k, such that its first off-diagonal   */
/*                 nonzero is l(k,j).  thus, the nonzero structure */
/*                 of column l(*,k) can be found by merging that   */
/*                 of such columns l(*,j) with the structure of    */
/*                 a(*,k).                                         */
/*        q - a vector of size neqns.  it is used to               */
/*                 accumulate the structure of each column l(*,k). */
/*                 at the end of the kth step, q(k),               */
/*                 q(q(k)), ..., is the list of                    */
/*                 positions of nonzeros in column k of the        */
/*                 factor l.                                       */
/*        clqtst - an integer vector of length neqns.  it is used  */
/*                 to test if mass symbolic elimination can be     */
/*                 performed.  that is, it is used to check        */
/*                 whether the structure of the current column k   */
/*                 being processed is completely determined by     */
/*                 the single column mask(k).                      */
/* **************************************************************** */

void smbfct(neqns,xadj,adjncy,perm,invp,xrnz,maxlnz,xnzsub,nzsub,maxsub,q,mask,clqtst,flag__,T)
int neqns, xadj, adjncy, perm, invp, xrnz, maxlnz, xnzsub, nzsub;
int maxsub, q, mask, clqtst, flag__, *T;
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int node, rchm, mrgk, lmax, i__, j, k, m, nabor, nzbeg, nzend; 
    static int kxsub, jstop, jstrt, mrkflg, np1, inz, knz;

/* **************************************************************** */

    /* Function Body */
    
#ifdef DEBUG        
    fprintf(fp1,"\n inside smbfct...");
    fprintf(fp1,"\n maxsub: %d, maxlnz: %d, flag__: %d ",maxsub,maxlnz,flag__);
#endif    
     
    if (neqns <= 0) {
	return;
    }

/*        ------------------ */
/*        initialization ... */
/*        ------------------ */
    nzbeg = 1;
    nzend = 0;
    T[xrnz+1] = 1;
    i__1 = neqns;
    for (k = 1; k <= i__1; k++) {
	T[mask+k] = 0;
	T[clqtst+k] = 0;
    }

/*        ---------------------------------------------- */
/*        for each column ... knz counts the number */
/*        of nonzeros in column k accumulated in q.      */
/*        ---------------------------------------------- */
    np1 = neqns + 1;
    i__1 = neqns;
    for (k = 1; k <= i__1; k++) {
	knz = 0;
	mrgk = T[mask+k];
	mrkflg = 0;
	T[clqtst+k] = k;
	if (mrgk != 0) {
	    T[clqtst+k] = T[clqtst+mrgk];
	}
	T[xnzsub+k] = nzend;
	node = T[perm+k];
	jstrt = T[xadj+node];
	jstop = T[xadj+node + 1] - 1;
	if (jstrt > jstop) {
	    goto L1600;
	}
/*                ----------------------------------- */
/*                use q to link through the           */
/*                structure of a(*,k) below diagonal. */
/*                ----------------------------------- */
	T[q+k] = np1;
	i__2 = jstop;
	for (j = jstrt; j <= i__2; j++) {
	    nabor = T[adjncy+j];
	    nabor = T[invp+nabor];
	    if (nabor <= k) {
		goto L300;
	    }
	    rchm = k;
L200:
	    m = rchm;
	    rchm = T[q+m];
	    if (rchm <= nabor) {
		goto L200;
	    }
	    knz++;
	    T[q+m] = nabor;
	    T[q+nabor] = rchm;
	    if (T[clqtst+nabor] != T[clqtst+k]) {
		mrkflg = 1;
	    }
L300:
	    ;
	}
/*                -------------------------------------- */
/*                test for mass symbolic elimination ... */
/*                -------------------------------------- */
	lmax = 0;
	if ((mrkflg != 0) || (mrgk == 0)) {
	    goto L400;
	}
	if (T[mask+mrgk] != 0) {
	    goto L400;
	}
	T[xnzsub+k] = T[xnzsub+mrgk] + 1;
	knz = T[xrnz+mrgk + 1] - (T[xrnz+mrgk] + 1);
	goto L1500;
L400:
/*                ----------------------------------------------- */
/*                link through each column i that affects l(*,k). */
/*                ----------------------------------------------- */
	i__ = k;
L500:
	i__ = T[mask+i__];
	if (i__ == 0) {
	    goto L900;
	}
	inz = T[xrnz+i__ + 1] - (T[xrnz+i__] + 1);
	jstrt = T[xnzsub+i__] + 1;
	jstop = T[xnzsub+i__] + inz;
	if (inz <= lmax) {
	    goto L600;
	}
	lmax = inz;
	T[xnzsub+k] = jstrt;
L600:
/*                        ---------------------------------- */
/*                        merge structure of l(*,i) in nzsub */
/*                        into q.                            */
/*                        ---------------------------------- */
	rchm = k;
	i__2 = jstop;
	for (j = jstrt; j <= i__2; j++) {
	    nabor = T[nzsub+j];
L700:
	    m = rchm;
	    rchm = T[q+m];
	    if (rchm < nabor) {
		goto L700;
	    }
	    if (rchm == nabor) {
		goto L800;
	    }
	    knz++;
	    T[q+m] = nabor;
	    T[q+nabor] = rchm;
	    rchm = nabor;
L800:
	    ;
	}
	goto L500;
L900:
/*                ----------------------------------- */
/*                check if subscripts duplicate those */
/*                of another column.                  */
/*                ----------------------------------- */
	if (knz == lmax) {
	    goto L1500;
	}
/*               ----------------------------------------------- */
/*               or if tail of k-1st column matches head of kth. */
/*               ----------------------------------------------- */
	if (nzbeg > nzend) {
	    goto L1300;
	}
	i__ = T[q+k];
	i__2 = nzend;
	for (jstrt = nzbeg; jstrt <= i__2; jstrt++) {
	    if ((i__3 = T[nzsub+jstrt] - i__) < 0) {
		goto L1000;
	    } else if (i__3 == 0) {
		goto L1100;
	    } else {
		goto L1300;
	    }
L1000:
	    ;
	}
	goto L1300;
L1100:
	T[xnzsub+k] = jstrt;
	i__2 = nzend;
	for (j = jstrt; j <= i__2; j++) {
	    if (T[nzsub+j] != i__) {
		goto L1300;
	    }
	    i__ = T[q+i__];
	    if (i__ > neqns) {
		goto L1500;
	    }
	}
	nzend = jstrt - 1;
L1300:
/*                    ---------------------------------------- */
/*                    copy the structure of l(*,k) from q      */
/*                    to the data structure (xnzsub, nzsub).   */
/*                    ---------------------------------------- */
	nzbeg = nzend + 1;
	nzend += knz;
	i__ = k;
	i__2 = nzend;
	for (j = nzbeg; j <= i__2; j++) {
	    i__ = T[q+i__];
	    T[nzsub+j] = i__;
	    T[clqtst+i__] = k;
	}
	T[xnzsub+k] = nzbeg;
	T[clqtst+k] = k;
L1500:
/*                --------------------------------------------- */
/*                update the vector mask.  note column l(*,k)   */
/*                just found is required to determine column    */
/*                l(*,j), where l(j,k) is the first nonzero in  */
/*                l(*,k) below diagonal.                        */
/*                --------------------------------------------- */
	if (knz <= 1) {
	    goto L1600;
	}
	kxsub = T[xnzsub+k];
	i__ = T[nzsub+kxsub];
	T[mask+k] = T[mask+i__];
	T[mask+i__] = k;
L1600:
	T[xrnz+k + 1] = T[xrnz+k] + knz;
    }
    maxlnz = T[xrnz+neqns] - 1;
    maxsub = T[xnzsub+neqns];
    spbcon_.NOFNZ=maxlnz;
    spbcon_.NOFSUB=maxsub;
    T[xnzsub+neqns + 1] = T[xnzsub+neqns];
    flag__ = 0;
    
#ifdef DEBUG        
    fprintf(fp1,"\n maxlnz: %d ,maxsub: %d ",maxlnz,maxsub);
#endif
    
    return;

} /* smbfct */

