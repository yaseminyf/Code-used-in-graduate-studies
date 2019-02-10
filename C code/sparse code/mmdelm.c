/* last revised 29.06.1998 */
/* mmdelm.c called by genmmd.c */

/* #include "data_struct.h" */

/* **     mmdelm ..... multiple minimum degree elimination     *** */
/*     purpose - this routine eliminates the node mdnode of        */
/*        minimum degree from the adjacency structure, which       */
/*        is stored in the quotient graph format.  it also         */
/*        transforms the quotient graph representation of the      */
/*        elimination graph.                                       */
/*     input parameters -                                          */
/*        mdnode - node of minimum degree.                         */
/*        maxint - estimate of maximum representable (short)       */
/*                 integer.                                        */
/*        tag    - tag value.                                      */
/*     updated parameters -                                        */
/*        (xadj,newadj) - updated adjacency structure.             */
/*        (dhead,invp,perm) - degree doubly linked structure.    */
/*        qsize  - size of supernode.                              */
/*        marker - marker vector.                                  */
/*        llist  - temporary linked list of eliminated nabors.     */
/* *************************************************************** */

void mmdelm(mdnode,xadj,newadj,dhead,invp,perm,qsize,llist,marker,maxint,tag,T)
int mdnode, xadj, newadj, dhead, invp, perm, qsize, llist;
int marker, maxint, tag, *T;
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int node, link, rloc, rlmt, i__, j, nabor, rnode, elmnt, xqnbr;
    static int istop, jstop, istrt, jstrt, nxnode, pvnode, nqnbrs, npv;

/* *************************************************************** */

/*        ----------------------------------------------- */
/*        find reachable set and place in data structure. */
/*        ----------------------------------------------- */

    /* Function Body */

#ifdef DEBUG    
    fprintf(fp1,"\n inside mmdelm...");
    fprintf(fp1,"\n mdnode: %d, tag: %d ",mdnode,tag);
#endif
    
    T[marker+mdnode] = tag;
    istrt = T[xadj+mdnode];
    istop = T[xadj+mdnode + 1] - 1;
/*        ------------------------------------------------------- */
/*        elmnt points to the beginning of the list of eliminated */
/*        nabors of mdnode, and rloc gives the storage location */
/*        for the next reachable node. */
/*        ------------------------------------------------------- */
    elmnt = 0;
    rloc = istrt;
    rlmt = istop;
    i__1 = istop;
    for (i__ = istrt; i__ <= i__1; i__++) {
	nabor = T[newadj+i__];
	if (nabor == 0) {
	    goto L300;
	}
	if (T[marker+nabor] >= tag) {
	    goto L200;
	}
	T[marker+nabor] = tag;
	if (T[invp+nabor] < 0) {
	    goto L100;
	}
	T[newadj+rloc] = nabor;
	rloc++;
	goto L200;
L100:
	T[llist+nabor] = elmnt;
	elmnt = nabor;
L200:
	;
    }
L300:
/*            ----------------------------------------------------- */
/*            merge with reachable nodes from generalized elements. */
/*            ----------------------------------------------------- */
    if (elmnt <= 0) {
	goto L1000;
    }
    T[newadj+rlmt] = -elmnt;
    link = elmnt;
L400:
    jstrt = T[xadj+link];
    jstop = T[xadj+link + 1] - 1;
    i__1 = jstop;
    for (j = jstrt; j <= i__1; j++) {
	node = T[newadj+j];
	link = -node;
	if (node < 0) {
	    goto L400;
	} else if (node == 0) {
	    goto L900;
	} else {
	    goto L500;
	}
L500:
	if ((T[marker+node] >= tag) || (T[invp+node] < 0)) {
	    goto L800;
	}
	T[marker+node] = tag;
/*                            --------------------------------- */
/*                            use storage from eliminated nodes */
/*                            if necessary. */
/*                            --------------------------------- */
L600:
	if (rloc < rlmt) {
	    goto L700;
	}
	link = -T[newadj+rlmt];
	rloc = T[xadj+link];
	rlmt = T[xadj+link + 1] - 1;
	goto L600;
L700:
	T[newadj+rloc] = node;
	rloc++;
L800:
	;
    }
L900:
    elmnt = T[llist+elmnt];
    goto L300;
L1000:
    if (rloc <= rlmt) {
	T[newadj+rloc] = 0;
    }
/*        -------------------------------------------------------- */
/*        for each node in the reachable set, do the following ... */
/*        -------------------------------------------------------- */
    link = mdnode;
L1100:
    istrt = T[xadj+link];
    istop = T[xadj+link + 1] - 1;
    i__1 = istop;
    for (i__ = istrt; i__ <= i__1; i__++) {
	rnode = T[newadj+i__];
	link = -rnode;
	if (rnode < 0) {
	    goto L1100;
	} else if (rnode == 0) {
	    goto L1800;
	} else {
	    goto L1200;
	}
L1200:
/*                -------------------------------------------- */
/*                if rnode is in the degree list structure ... */
/*                -------------------------------------------- */
	pvnode = T[perm+rnode];
	if ((pvnode == 0) || (pvnode == -(maxint))) {
	    goto L1300;
	}
/*                    ------------------------------------- */
/*                    then remove rnode from the structure. */
/*                    ------------------------------------- */
	nxnode = T[invp+rnode];
	if (nxnode > 0) {
	    T[perm+nxnode] = pvnode;
	}
	if (pvnode > 0) {
	    T[invp+pvnode] = nxnode;
	}
	npv = -pvnode;
	if (pvnode < 0) {
	    T[dhead+npv] = nxnode;
	}
L1300:
/*                ---------------------------------------- */
/*                purge inactive quotient nabors of rnode. */
/*                ---------------------------------------- */
	jstrt = T[xadj+rnode];
	jstop = T[xadj+rnode + 1] - 1;
	xqnbr = jstrt;
	i__2 = jstop;
	for (j = jstrt; j <= i__2; j++) {
	    nabor = T[newadj+j];
	    if (nabor == 0) {
		goto L1500;
	    }
	    if (T[marker+nabor] >= tag) {
		goto L1400;
	    }
	    T[newadj+xqnbr] = nabor;
	    xqnbr++;
L1400:
	    ;
	}
L1500:
/*                ---------------------------------------- */
/*                if no active nabor after the purging ... */
/*                ---------------------------------------- */
	nqnbrs = xqnbr - jstrt;
	if (nqnbrs > 0) {
	    goto L1600;
	}
/*                    ----------------------------- */
/*                    then merge rnode with mdnode. */
/*                    ----------------------------- */
	T[qsize+mdnode] += T[qsize+rnode];
	T[qsize+rnode] = 0;
	T[marker+rnode] = maxint;
	T[invp+rnode] = -mdnode;
	T[perm+rnode] = -maxint;
	goto L1700;
L1600:
/*                -------------------------------------- */
/*                else flag rnode for degree update, and */
/*                add mdnode as a nabor of rnode. */
/*                -------------------------------------- */
	T[invp+rnode] = nqnbrs + 1;
	T[perm+rnode] = 0;
	T[newadj+xqnbr] = mdnode;
	xqnbr++;
	if (xqnbr <= jstop) {
	    T[newadj+xqnbr] = 0;
	}

L1700:
	;
    }
L1800:
#ifdef DEBUG    
    fprintf(fp1,"\n tag: %d ",tag);
#endif    
    return;

} /* mmdelm */

