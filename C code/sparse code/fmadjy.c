/* last revised 28.06.1998 */
/* fmadjy.c called by ipendb.c */

/* #include "data_struct.h" */

/* ********     fmadjy ..... form adjacency structure     ******** */
/*     purpose - to form the entire adjacency structure from       */
/*        the lower adjacency structure.                           */
/*     input parameters -                                          */
/*        neqns  - number of equations.                            */
/*        nuladj   - lower adjacency structure, where neighbors of */
/*                 a node are stored sequentially in nuladj.       */
/*        nudeg    - the degree of the variables.                  */
/*     updated parameter -                                         */
/*        xadj   - on input, contains the degree in the lower      */
/*                 adjacency structure.  on output, the index      */
/*                 to adjncy.                                      */
/*     output parameter -                                          */
/*        adjncy - the adjacency structure vector.                 */
/* *************************************************************** */

void fmadjy(neqns, T, xadj, adjncy, nuladj, nudeg)
int neqns, *T, xadj, adjncy, nuladj, nudeg;
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int next, i__, j, iadj, nabor, jstop, nuxadj, nbrsze;

/* *************************************************************** */

    /* Function Body */

#ifdef DEBUG        
    fprintf(fp1,"\n inside fmadjy...");
    fprintf(fp1,"\n neqns: %d, xadj: %d, adjncy: %d, nuladj: %d, nudeg: %d ",neqns,xadj,adjncy,nuladj,nudeg);
#endif
    
    if (neqns <= 0) {
	return;
    }

    next = 1;
    iadj = 0;
/*        ------------------------------------------------------ */
/*        move neighbors from the lower adjacency structure to   */
/*        the adjacency structure.  allow room for the neighbors */
/*        in the upper adjacency structure.                      */
/*        ------------------------------------------------------ */
    i__1 = neqns;
    for (i__ = 1; i__ <= i__1; i__++) {
	nbrsze = T[xadj+i__];
	jstop = next + nbrsze - 1;
	T[xadj+i__] = jstop;
	if (nbrsze <= 0) {
	    goto L200;
	}
	i__2 = jstop;
	for (j = next; j <= i__2; j++) {
	    iadj++;
	    nabor = T[nuladj+iadj];
	    T[adjncy+j] = nabor;
	    nuxadj = T[xadj+nabor] + 1;
	    T[xadj+nabor] = nuxadj;
	    T[adjncy+nuxadj] = i__;
	}
L200:
	next += T[nudeg+i__];
    }

    if (neqns == 1) {
	goto L500;
    }
    j = neqns + 1;
    i__1 = neqns;
    for (i__ = 1; i__ <= i__1; i__++) {
	T[xadj+j] = T[xadj+j - 1] + 1;
	j--;
    }
    T[xadj+1] = 1;
    return;

L500:
    T[xadj+1] = 1;
    T[xadj+2] = 1;
    
#ifdef DEBUG        
    fprintf(fp1,"\n xadj after fmadjy..\n");
    i__1=spbcon_.NCOLS+1;
    for (i__=1;i__<=i__1;i__++) {
      fprintf(fp1,"xadj[%d]: %d ",i__,T[xadj+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
#endif
    
    return;

} /* fmadjy */

