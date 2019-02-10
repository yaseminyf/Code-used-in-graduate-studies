/* last revised 29.06.1998 */
/* build.c called by ipendb.c */

/* #include "data_struct.h" */

/* *****     build ..... build lower adjacency structure     ***** */
/*     purpose - to build a lower adjacency structure from         */
/*        adjacency linked lists.  the transformation can be       */
/*        performed in place, but more space in the storage        */
/*        vector will speed up the task.                           */
/*     input parameters -                                          */
/*        neqns  - number of equations.                            */
/*        nedges - number of edges in the lower triangle.          */
/*        iladj  - the starting point for the lower adjacency      */
/*                 structure in the array t.                       */
/*        limit - last free storage before linked lists in t.      */
/*     updated parameter -                                         */
/*        t      - the storage vector.  the first neqns            */
/*                 locations, on input, contains the header for    */
/*                 each of the linked lists, and on output, the    */
/*                 degree of the nodes in the lower adjacency      */
/*                 structure.                                      */
/*     output parameter -                                          */
/*        deg    - the degree of each node.  (an array of size     */
/*                 neqns.)                                         */
/* *************************************************************** */

void build(neqns, nedges, T, deg, iladj, limit)
int neqns, nedges, *T, deg, iladj, limit;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int ideg, link, next, i__, nabor, hlink, ohlink, nulink, 
	    maxtop, top;

/* *************************************************************** */

    /* Function Body */

#ifdef DEBUG    
    fprintf(fp1,"\n inside build..");  
    fprintf(fp1,"\n neqns: %d, nedges: %d, iladj: %d, limit: %d, deg: %d ",neqns,nedges,iladj,limit,deg);
#endif
    
    if (neqns <= 0) {
	return;
    }

/*        --------------------------------------------------- */
/*        initialization ...                                  */
/*        next points to the position where a neighbor can be */
/*        placed.  top points to the position of the current  */
/*        top (link,nbr) pair in the vector s.                */
/*        --------------------------------------------------- */
    next = iladj;
    top = limit + 1;
    maxtop = iladj + nedges;

    T[1] = 0; 
    T[deg+1] = 0;

#ifdef DEBUG    
    fprintf(fp1,"\n initialized T[1] and deg[1].."); 
#endif
     
    if (neqns == 1) {
	return;
    }

/*        ---------------------------------------------- */
/*        create lower adjacency lists for each node ... */
/*        ---------------------------------------------- */
    i__1 = neqns;
    for (i__ = 2; i__ <= i__1; i__++) {
	ideg = 0;
	nulink = T[i__];

/*            ---------------------------------------------- */
/*            for each node in the linked list of node i ... */
/*            ---------------------------------------------- */
L100:
	link = nulink;
	if (link <= 0) {
	    goto L700;
	}
	nulink = T[link];
	nabor = T[link + 1];

	if (next >= top && link != top) {
	    goto L200;
	}
/*                    ---------------------------------------- */
/*                    free the space for (link,nabor) pair ... */
/*                    ---------------------------------------- */
	T[link] = 0;
	T[link + 1] = 0;

	if (link != top) {
	    goto L600;
	}
	goto L500;

L200:
/*                ------------------------------------------------ */
/*                collision occurs -- clear off some space for the */
/*                nabor by moving the (link,nbr) pair at the top   */
/*                downwards and adjust the associated link fields. */
/*                ------------------------------------------------ */
	T[i__] = link;
	hlink = top;
L300:
	hlink = T[hlink];
	if (hlink > 0) {
	    goto L300;
	}

	hlink = -hlink;
L400:
	ohlink = hlink;
	hlink = T[ohlink];
	if (hlink != top) {
	    goto L400;
	}

	T[link] = T[top];
	T[link + 1] = T[top + 1];
	if (ohlink != link) {
	    T[ohlink] = link;
	}
	if (nulink == top) {
	    nulink = link;
	}

L500:
	top += 2;
	if (top >= maxtop) {
	    goto L600;
	}
	if (T[top] != 0) {
	    goto L600;
	}
	goto L500;

L600:
/*                ----------------------------------------- */
/*                put nabor in the neighbor list of node i. */
/*                ----------------------------------------- */
	T[next] = nabor;
	next++;
	ideg++;
	T[deg+nabor]++;
	goto L100;

L700:
/*            ------------------- */
/*            set the degrees ... */
/*            ------------------- */
	T[deg+i__] = ideg;
	T[i__] = ideg;

/* L800: */
    }

#ifdef DEBUG    
    fprintf(fp1,"\n deg after build...\n");
    for (i__=1;i__<=neqns;i__++) {
      fprintf(fp1,"%d ",T[deg+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
#endif

    return;

} /* build */

