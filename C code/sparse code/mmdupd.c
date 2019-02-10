/* last revised 29.06.1998 */
/* mmdupd.c called by genmmd.c */

#include "data_struct.h" 

/* *****     mmdupd ..... multiple minimum degree update     ***** */
/*     purpose - this routine updates the degrees of nodes         */
/*        after a multiple elimination step.                       */
/*     input parameters -                                          */
/*        ehead  - the beginning of the list of eliminated         */
/*                 nodes (i.e., newly formed elements).            */
/*        neqns  - number of equations.                            */
/*        (xadj,newadj) - adjacency structure.                     */
/*        delta  - tolerance value for multiple elimination.       */
/*        maxint - maximum machine representable (short)           */
/*                 integer.                                        */
/*     updated parameters -                                        */
/*        mdeg   - new minimum degree after degree update.         */
/*        (dhead,invp,perm) - degree doubly linked structure.      */
/*        qsize  - size of supernode.                              */
/*        llist  - working linked list.                            */
/*        marker - marker vector for degree update.                */
/*        tag    - tag value.                                      */
/* *************************************************************** */

void mmdupd(ehead,neqns,xadj,newadj,delta,mdeg,dhead,invp,perm,qsize,llist,marker,maxint,tag,T)
int ehead, neqns, xadj, newadj, delta, mdeg, dhead, invp, perm;
int qsize, llist, marker, maxint, tag, *T;
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int node, mtag, link, mdeg0, i__, j, enode, fnode, nabor; 
    static int elmnt, istop, jstop, q2head, istrt, jstrt, qxhead;
    static int iq2, deg, deg0;

/* *************************************************************** */

    /* Function Body */
    
#ifdef DEBUG        
    fprintf(fp1,"\ninside mmdupd..");
#endif
    
    mdeg0 = mdeg + delta;
    elmnt = ehead;
L100:
/*            ------------------------------------------------------- */
/*            for each of the newly formed element, do the following. */
/*            (reset tag value if necessary.)                         */
/*            ------------------------------------------------------- */
    if (elmnt <= 0) {
        spkglob_.tag=tag;
        
#ifdef DEBUG            
        fprintf(fp1,"\n spkglob_.tag: %d ",spkglob_.tag);
        fprintf(fp1,"\n spkglob_.mdeg: %d ",spkglob_.mdeg);
#endif
        
	return;
    }
    mtag = tag + mdeg0;
    if (mtag < maxint) {
	goto L300;
    }
    tag=1;
    i__1 = neqns;
    for (i__ = 1; i__ <= i__1; i__++) {
	if (T[marker+i__] < maxint) {
	    T[marker+i__] = 0;
	}
    }
    mtag = tag + mdeg0;
L300:
/*            --------------------------------------------- */
/*            create two linked lists from nodes associated */
/*            with elmnt: one with two nabors (q2head) in */
/*            adjacency structure, and the other with more */
/*            than two nabors (qxhead).  also compute deg0, */
/*            number of nodes in this element. */
/*            --------------------------------------------- */
    q2head = 0;
    qxhead = 0;
    deg0 = 0;
    link = elmnt;
L400:
    istrt = T[xadj+link];
    istop = T[xadj+link + 1] - 1;
    i__1 = istop;
    for (i__ = istrt; i__ <= i__1; i__++) {
	enode = T[newadj+i__];
	link = -enode;
	if (enode < 0) {
	    goto L400;
	} else if (enode == 0) {
	    goto L800;
	} else {
	    goto L500;
	}

L500:
	if (T[qsize+enode] == 0) {
	    goto L700;
	}
	deg0 += T[qsize+enode];
	T[marker+enode] = mtag;
/*                        ---------------------------------- */
/*                        if enode requires a degree update, */
/*                        then do the following.             */
/*                        ---------------------------------- */
	if (T[perm+enode] != 0) {
	    goto L700;
	}
/*                    --------------------------------------- */
/*                    place either in qxhead or q2head lists. */
/*                    --------------------------------------- */
	if (T[invp+enode] == 2) {
	    goto L600;
	}
	T[llist+enode] = qxhead;
	qxhead = enode;
	goto L700;
L600:
	T[llist+enode] = q2head;
	q2head = enode;
L700:
	;
    }
L800:
/*            -------------------------------------------- */
/*            for each enode in q2 list, do the following. */
/*            -------------------------------------------- */
    enode = q2head;
    iq2 = 1;
L900:
    if (enode <= 0) {
	goto L1500;
    }
    if (T[perm+enode] != 0) {
	goto L2200;
    }
    tag++;
    deg = deg0;
/*                    ------------------------------------------ */
/*                    identify the other adjacent element nabor. */
/*                    ------------------------------------------ */
    istrt = T[xadj+enode];
    nabor = T[newadj+istrt];
    if (nabor == elmnt) {
	nabor = T[newadj+istrt + 1];
    }
/*                    ------------------------------------------------ */
/*                    if nabor is uneliminated, increase degree count. */
/*                    ------------------------------------------------ */
    link = nabor;
    if (T[invp+nabor] < 0) {
	goto L1000;
    }
    deg += T[qsize+nabor];
    goto L2100;
L1000:
/*                        -------------------------------------------- */
/*                        otherwise, for each node in the 2nd element, */
/*                        do the following.                            */
/*                        -------------------------------------------- */
    istrt = T[xadj+link];
    istop = T[xadj+link + 1] - 1;
    i__1 = istop;
    for (i__ = istrt; i__ <= i__1; i__++) {
	node = T[newadj+i__];
	link = -node;
	if (node == enode) {
	    goto L1400;
	}
	if (node < 0) {
	    goto L1000;
	} else if (node == 0) {
	    goto L2100;
	} else {
	    goto L1100;
	}

L1100:
	if (T[qsize+node] == 0) {
	    goto L1400;
	}
	if (T[marker+node] >= tag) {
	    goto L1200;
	}
/*                       ------------------------------------- */
/*                       case when node is not yet considered. */
/*                       ------------------------------------- */
	T[marker+node] = tag;
	deg += T[qsize+node];
	goto L1400;
L1200:
/*                       ---------------------------------------- */
/*                       case when node is indistinguishable from */
/*                       enode.  merge them into a new supernode. */
/*                       ---------------------------------------- */
	if (T[perm+node] != 0) {
	    goto L1400;
	}
	if (T[invp+node] != 2) {
	    goto L1300;
	}
	T[qsize+enode] += T[qsize+node];
	T[qsize+node] = 0;
	T[marker+node] = maxint;
	T[invp+node] = -enode;
	T[perm+node] = -maxint;
	goto L1400;
L1300:
/*                      -------------------------------------- */
/*                      case when node is outmatched by enode. */
/*                      -------------------------------------- */
	if (T[perm+node] == 0) {
	    T[perm+node] = -maxint;
	}
L1400:
	;
    }
    goto L2100;
L1500:
/*                ------------------------------------------------ */
/*                for each enode in the qx list, do the following. */
/*                ------------------------------------------------ */
    enode = qxhead;
    iq2 = 0;
L1600:
    if (enode <= 0) {
	goto L2300;
    }
    if (T[perm+enode] != 0) {
	goto L2200;
    }
    tag++;
    deg = deg0;
/*                        --------------------------------- */
/*                        for each unmarked nabor of enode, */
/*                        do the following.                 */
/*                        --------------------------------- */
    istrt = T[xadj+enode];
    istop = T[xadj+enode + 1] - 1;
    i__1 = istop;
    for (i__ = istrt; i__ <= i__1; i__++) {
	nabor = T[newadj+i__];
	if (nabor == 0) {
	    goto L2100;
	}
	if (T[marker+nabor] >= tag) {
	    goto L2000;
	}
	T[marker+nabor] = tag;
	link = nabor;
/*                           ------------------------------ */
/*                           if uneliminated, include it in */
/*                           deg count.                     */
/*                           ------------------------------ */
	if (T[invp+nabor] < 0) {
	    goto L1700;
	}
	deg += T[qsize+nabor];
	goto L2000;
L1700:
/*                           ------------------------------- */
/*                           if eliminated, include unmarked */
/*                           nodes in this element into the  */
/*                           degree count.                   */
/*                           ------------------------------- */
	jstrt = T[xadj+link];
	jstop = T[xadj+link + 1] - 1;
	i__2 = jstop;
	for (j = jstrt; j <= i__2; j++) {
	    node = T[newadj+j];
	    link = -node;
	    if (node < 0) {
		goto L1700;
	    } else if (node == 0) {
		goto L2000;
	    } else {
		goto L1800;
	    }

L1800:
	    if (T[marker+node] >= tag) {
		goto L1900;
	    }
	    T[marker+node] = tag;
	    deg += T[qsize+node];
L1900:
	    ;
	}
L2000:
	;
    }
L2100:
/*                    ------------------------------------------- */
/*                    update external degree of enode in degree */
/*                    structure, and mdeg (min deg) if necessary. */
/*                    ------------------------------------------- */
    deg = deg - T[qsize+enode] + 1;
    fnode = T[dhead+deg];
    T[invp+enode] = fnode;
    T[perm+enode] = -deg;
    if (fnode > 0) {
	T[perm+fnode] = enode;
    }
    T[dhead+deg] = enode;
    if (deg < mdeg) {
        mdeg = deg;
	spkglob_.mdeg = deg;
    }
L2200:
/*                    ---------------------------------- */
/*                    get next enode in current element. */
/*                    ---------------------------------- */
    enode = T[llist+enode];
    if (iq2 == 1) {
	goto L900;
    }
    goto L1600;
L2300:
/*            ----------------------------- */
/*            get next element in the list. */
/*            ----------------------------- */
    tag=mtag;
    elmnt = T[llist+elmnt];
    goto L100;

} /* mmdupd */

