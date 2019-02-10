/* last revised 29.06.1998 */
/* genmmd.c called by orcolb.c */

#include "data_struct.h" 


/* ****     genmmd ..... multiple minimum external degree     **** */
/*     purpose - this routine implements the minimum degree        */
/*        algorithm.  it makes use of the implicit representation  */
/*        of elimination graphs by quotient graphs, and the        */
/*        notion of indistinguishable nodes.  it also implements   */
/*        the modifications by multiple elimination and minimum    */
/*        external degree.                                         */
/*        ---------------------------------------------            */
/*        caution - the adjacency vector newadj will be            */
/*        destroyed.                                               */
/*        ---------------------------------------------            */
/*     input parameters -                                          */
/*        neqns  - number of equations.                            */
/*        (xadj,newadj) - the adjacency structure.                 */
/*        delta  - tolerance value for multiple elimination.       */
/*        maxint - maximum machine representable (short) integer   */
/*                 (any smaller estimate will do) for marking      */
/*                 nodes.                                          */
/*     output parameters -                                         */
/*        perm   - the minimum degree ordering.                    */
/*        invp   - the inverse of perm.                            */
/*        nofsub - an upper bound on the number of nonzero         */
/*                 subscripts for the compressed storage scheme.   */
/*     working parameters -                                        */
/*        dhead  - vector for head of degree lists.                */
/*        invp   - used temporarily for degree forward link.       */
/*        perm   - used temporarily for degree backward link.      */
/*        qsize  - vector for size of supernodes.                  */
/*        llist  - vector for temporary linked lists.              */
/*        marker - a temporary marker vector.                      */
/*     program subroutines -                                       */
/*        mmdelm, mmdint, mmdnum, mmdupd.                          */
/* *************************************************************** */

void genmmd(neqns,xadj,newadj,invp,perm,delta,dhead,qsize,llist,marker,maxint,nofsub,T)
int neqns, xadj, newadj, invp, perm, delta, dhead, qsize, llist; 
int marker, maxint, nofsub, *T;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int  i__, mdlmt;
    static int nextmd, num;

/* *************************************************************** */

    /* Function Body */
    
#ifdef DEBUG        
    fprintf(fp1,"\n inside genmmd ...");
    fprintf(fp1,"\n neqns: %d, xadj: %d, newadj: %d ",neqns,xadj,newadj);
    fprintf(fp1,"\n invp: %d, perm: %d, delta: %d ",invp,perm,delta);
    fprintf(fp1,"\n dhead: %d, qsize: %d, llist: %d ",dhead,qsize,llist);
    fprintf(fp1,"\n marker: %d, maxint: %d, nofsub: %d ",marker,maxint,nofsub);
    fprintf(fp1,"\n perm vector at the start of genmmd...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[perm+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n invp vector at the start of genmmd...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[invp+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
#endif
    
    if (neqns <= 0) {
	return;
    }

/*        ------------------------------------------------ */
/*        initialization for the minimum degree algorithm. */
/*        ------------------------------------------------ */

    nofsub = 0;
    
#ifdef DEBUG        
    fprintf(fp1,"\n before mmdint...");
#endif
    
    mmdint(neqns,xadj,newadj,dhead,invp,perm,qsize,llist,marker,T);
    
#ifdef DEBUG        
    fprintf(fp1,"\n after mmdint...");
    fprintf(fp1,"\n invp vector after mmdint...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[invp+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n perm vector after mmdint...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[perm+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n dhead vector after mmdint...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[dhead+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n qsize vector after mmdint...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[qsize+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n llist vector after mmdint...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[llist+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n marker vector after mmdint...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[marker+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
#endif
    
/*        ---------------------------------------------- */
/*        num counts the number of ordered nodes plus 1. */
/*        ---------------------------------------------- */
    num = 1;

/*        ----------------------------- */
/*        eliminate all isolated nodes. */
/*        ----------------------------- */
    nextmd = T[dhead+1];
L100:
    if (nextmd <= 0) {
	goto L200;
    }
    spkglob_.mdnode = nextmd;
    nextmd = T[invp+spkglob_.mdnode];
    T[marker+spkglob_.mdnode] = maxint;
    T[invp+spkglob_.mdnode] = -num;
    num++;
    goto L100;

L200:
/*        ---------------------------------------- */
/*        search for node of the minimum degree.   */
/*        mdeg is the current minimum degree;      */
/*        tag is used to facilitate marking nodes. */
/*        ---------------------------------------- */

    if (num > neqns) {
	goto L1000;
    }
    spkglob_.tag = 1;
    T[dhead+1] = 0;
    spkglob_.mdeg = 2;
L300:
    if (T[dhead+spkglob_.mdeg] > 0) {
	goto L400;
    }
    spkglob_.mdeg++;  
    goto L300;
L400:
/*            ------------------------------------------------- */
/*            use value of delta to set up mdlmt, which governs */
/*            when a degree update is to be performed.          */
/*            ------------------------------------------------- */

#ifdef DEBUG    
    fprintf(fp1,"\n delta: %d ",delta);
#endif
    
    mdlmt = spkglob_.mdeg + delta;
    spkglob_.ehead = 0;

L500:
    spkglob_.mdnode = T[dhead+spkglob_.mdeg];
    if (spkglob_.mdnode > 0) {
	goto L600;
    }
    spkglob_.mdeg++;
    if (spkglob_.mdeg > mdlmt) {
	goto L900;
    }
    goto L500;
L600:
/*                ---------------------------------------- */
/*                remove mdnode from the degree structure. */
/*                ---------------------------------------- */
    nextmd = T[invp+spkglob_.mdnode];
    T[dhead+spkglob_.mdeg] = nextmd;
    if (nextmd > 0) {
	T[perm+nextmd] = -spkglob_.mdeg;
    }
    T[invp+spkglob_.mdnode] = -num;
    nofsub = nofsub + spkglob_.mdeg + T[qsize+spkglob_.mdnode] - 2;
    spbcon_.NOFSUB=nofsub; 
    
#ifdef DEBUG        
    fprintf(fp1,"\n spbcon_.NOFSUB: %d ",spbcon_.NOFSUB);
#endif
    
    if ((num + T[qsize+spkglob_.mdnode]) > neqns) {
	goto L1000;
    }
/*                ---------------------------------------------- */
/*                eliminate mdnode and perform quotient graph    */
/*                transformation.  reset tag value if necessary. */
/*                ---------------------------------------------- */
    spkglob_.tag++;
    if (spkglob_.tag < maxint) {
	goto L800;
    }
    spkglob_.tag = 1;
    i__1 = neqns;
    for (i__ = 1; i__ <= i__1; i__++) {
	if (T[marker+i__] < maxint) {
	    T[marker+i__] = 0;
	}
    }
L800:

#ifdef DEBUG    
    fprintf(fp1,"\n right before mmdelm..."); 
    fprintf(fp1,"\n spkglob_.tag: %d, spkglob_.mdnode: %d ",spkglob_.tag,spkglob_.mdnode);
#endif
    
    mmdelm(spkglob_.mdnode,xadj,newadj,dhead,invp,perm,qsize,llist,marker,maxint,spkglob_.tag,T);
    
#ifdef DEBUG        
    fprintf(fp1,"\n right after mmdelm...");
    fprintf(fp1,"\n invp vector after mmdelm...\n"); 
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[invp+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n perm vector after mmdelm...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[perm+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n dhead vector after mmdelm...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[dhead+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n qsize vector after mmdelm...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[qsize+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n llist vector after mmdelm...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1," %d ",T[llist+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n marker vector after mmdelm...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[marker+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n xadj vector after mmdelm...\n");
    i__1=neqns+1;
    for (i__=1;i__<=i__1;i__++) {
      fprintf(fp1," %d ",T[xadj+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n spbcon_.NEDGES: %d ",spbcon_.NEDGES);
    fprintf(fp1,"\n newadj vector after mmdelm...\n");
    i__1=2*spbcon_.NEDGES;
    for (i__=1;i__<=i__1;i__++) {
      fprintf(fp1,"%d ",T[newadj+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
#endif
    
    num += T[qsize+spkglob_.mdnode];
    T[llist+spkglob_.mdnode] = spkglob_.ehead;
    spkglob_.ehead = spkglob_.mdnode;
    if (delta >= 0) {
	goto L500;
    }
L900:
/*            ------------------------------------------- */
/*            update degrees of the nodes involved in the */
/*            minimum degree nodes elimination. */
/*            ------------------------------------------- */
    if (num > neqns) {
	goto L1000;
    }

#ifdef DEBUG        
    fprintf(fp1,"\n right before mmdupd...");
    fprintf(fp1,"\n  spkglob_.tag: %d, spkglob_.ehead: %d, spkglob_.mdeg: %d ",spkglob_.tag,spkglob_.ehead,spkglob_.mdeg);
#endif
    
    mmdupd(spkglob_.ehead,neqns,xadj,newadj,delta,spkglob_.mdeg,dhead,invp,perm,qsize,llist,marker,maxint,spkglob_.tag,T);
    
#ifdef DEBUG        
    fprintf(fp1,"\n right after mmdupd..."); 
    fprintf(fp1,"\n spkglob_.tag: %d ",spkglob_.tag);
    fprintf(fp1,"\n spkglob_.mdeg: %d ",spkglob_.mdeg);
    fprintf(fp1,"\n invp vector after mmdupd...\n"); 
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[invp+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n perm vector after mmdupd...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[perm+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n dhead vector after mmdupd...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[dhead+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n qsize vector after mmdupd...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[qsize+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n llist vector after mmdupd...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[llist+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n marker vector after mmdupd...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[marker+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
#endif
    
    goto L300;

L1000:

#ifdef DEBUG    
    fprintf(fp1,"\n right before mmdnum..."); 
#endif
    
    mmdnum(neqns,perm,invp,qsize,T);
    
#ifdef DEBUG        
    fprintf(fp1,"\n right after mmdnum..."); 
    fprintf(fp1,"\n invp vector after mmdnum...\n"); 
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[invp+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n perm vector after mmdnum...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[perm+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
#endif
    
    return;

} /* genmmd */

