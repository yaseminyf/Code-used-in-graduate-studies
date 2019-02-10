/* last revised 29.06.1998 */
/* mmdint.c called by genmmd.c */

/* #include "data_struct.h" */

/* ***     mmdint ..... mult minimum degree initialization     *** */
/*     purpose - this routine performs initialization for the      */
/*        multiple elimination version of the minimum degree       */
/*        algorithm.                                               */
/*     input parameters -                                          */
/*        neqns  - number of equations.                            */
/*        (xadj,newadj) - adjacency structure.                     */
/*     output parameters -                                         */
/*        (dhead,invp,perm) - degree doubly linked structure.    */
/*        qsize  - size of supernode (initialized to one).         */
/*        llist  - linked list.                                    */
/*        marker - marker vector.                                  */
/* *************************************************************** */

void mmdint(neqns,xadj,newadj,dhead,invp,perm,qsize,llist,marker,T)
int neqns, xadj, newadj, dhead, invp, perm, qsize, llist, marker, *T;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int ndeg, node, fnode;

/* *************************************************************** */

    /* Function Body */

#ifdef DEBUG        
    fprintf(fp1,"\n inside mmdint...");
#endif
    
    i__1 = neqns;
    
    for (node = 1; node <= i__1; node++) {
        T[dhead+node]=0;
	T[qsize+node] = 1;
	T[marker+node] = 0;
	T[llist+node] = 0;
    }
    
/*        ------------------------------------------ */
/*        initialize the degree doubly linked lists. */
/*        ------------------------------------------ */
    i__1 = neqns;
    for (node = 1; node <= i__1; node++) {
	ndeg = T[xadj+node + 1] - T[xadj+node] + 1;
        fnode = T[dhead+ndeg]; 
	T[invp+node] = fnode;
	T[dhead+ndeg] = node; 
	if (fnode > 0) {
	    T[perm+fnode] = node;
	}
	T[perm+node] = -ndeg;
/* L200: */
    }

    return;

} /* mmdint */
