/* last revised 29.06.1998 */
/* mmdnum.c called by genmmd.c */

/* #include "data_struct.h" */

/* *****     mmdnum ..... multi minimum degree numbering     ***** */
/*     purpose - this routine performs the final step in           */
/*        producing the permutation and inverse permutation        */
/*        vectors in the multiple elimination version of the       */
/*        minimum degree ordering algorithm.                       */
/*     input parameters -                                          */
/*        neqns  - number of equations.                            */
/*        qsize  - size of supernodes at elimination.              */
/*     updated parameters -                                        */
/*        invp   - inverse permutation vector.  on input,          */
/*                 if qsize(node)=0, then node has been merged     */
/*                 into the node -invp(node); otherwise,           */
/*                 -invp(node) is its inverse labelling.           */
/*     output parameters -                                         */
/*        perm   - the permutation vector.                         */
/* *************************************************************** */

void mmdnum(neqns, perm, invp, qsize, T)
int neqns, perm, invp, qsize, *T;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int node, root, nextf, father, nqsize, num;

/* *************************************************************** */

    /* Function Body */
    
#ifdef DEBUG        
    fprintf(fp1,"\n inside mmdnum...");
#endif
    
    i__1 = neqns;
    for (node = 1; node <= i__1; node++) {
	nqsize = T[qsize+node];
	if (nqsize <= 0) {
	    T[perm+node] = T[invp+node];
	}
	if (nqsize > 0) {
	    T[perm+node] = -T[invp+node];
	}
    }
/*        ------------------------------------------------------ */
/*        for each node which has been merged, do the following. */
/*        ------------------------------------------------------ */
    i__1 = neqns;
    for (node = 1; node <= i__1; node++) {
	if (T[perm+node] > 0) {
	    goto L500;
	}
/*                ----------------------------------------- */
/*                trace the merged tree until one which has */
/*                not been merged, call it root. */
/*                ----------------------------------------- */
	father = node;
L200:
	if (T[perm+father] > 0) {
	    goto L300;
	}
	father = -T[perm+father];
	goto L200;
L300:
/*                ----------------------- */
/*                number node after root. */
/*                ----------------------- */
	root = father;
	num = T[perm+root] + 1;
	T[invp+node] = -num;
	T[perm+root] = num;
/*                ------------------------ */
/*                shorten the merged tree. */
/*                ------------------------ */
	father = node;
L400:
	nextf = -T[perm+father];
	if (nextf <= 0) {
	    goto L500;
	}
	T[perm+father] = -root;
	father = nextf;
	goto L400;
L500:
	;
    }
/*        ---------------------- */
/*        ready to compute perm. */
/*        ---------------------- */
    i__1 = neqns;
    for (node = 1; node <= i__1; node++) {
	num = -T[invp+node];
	T[invp+node] = num;
	T[perm+num] = node;
    }
    return;

} /* mmdnum */

