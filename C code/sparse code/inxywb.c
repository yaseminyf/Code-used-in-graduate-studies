/* last revision 18.06.1998 */
/* inxywb.c called by main program */

#include "data_struct.h" 


/* *****     inxywb ..... input equations and constraints     **** */
/*     purpose - the user supplies an equation or a constraint     */
/*        to the package by calling this routine which will        */
/*        write it out onto an external file.                      */
/*        (the part where the matrix is written to an external     */
/*        file is no longer present)                               */
/*        at the same                                              */
/*        time, insubs is called to create the linked              */
/*        structure representing the structure of                  */
/*        (a-transpose)(a).                                        */
/*     input parameters -                                          */
/*        ROWNUM - an integer associated with the equation or      */
/*                 constraint which determines an initial          */
/*                 ordering for the equations and constraints.     */
/*        type   - type of equation or constraint:                 */
/*                 1 - sparse equation.                            */
/*                 2 - dense equation.                             */
/*                 3 - sparse constraint.                          */
/*                 4 - dense constraint.                           */
/*        NSUBS  - number of nonzeros in the incoming equation     */
/*                 or constraint.                                  */
/*        SUBS   - column subscripts of nonzeros.                  */
/*        VALUES - numerical values of nonzeros.                   */
/*        Rhs    - right hand side component.                      */
/*        WEIGHT - the diagonal element of the weight matrix       */
/*                 corresponding to the incoming row.  if the      */
/*                 problem is not weighted, weight should be       */
/*                 set to one.                                     */
/*     updated parameter -                                         */
/*        t      - the storage array of size maxsb.                */
/* *************************************************************** */

void inxywb(ROWNUM, TYPE, NSUBS, SUBS, T)
int ROWNUM, TYPE, NSUBS, *SUBS;
int *T;
{
    /* System generated locals */
    
    int i__1, i__2;

    /* Local variables */
    
    static int link, imax, isub, jsub, i__, j;
    static int nabor, olink, isubm1;

/* *************************************************************** */
    
/*            ------------------------------------------ */
/*            first time inxywb is called, print heading */
/*            and initialize common blocks and storage */
/*            array. */
/*            ------------------------------------------ */

    /* Function Body */

/* this routine is not needed since we don't write in file */
/* other things done in this routine are done in the main program */
/*   initpb(t);  */

/* this routine is not needed since the entrance of the rows is done in */
/* ascending order */
/*    sortlr(NSUBS, SUBS, VALUES); */

    imax = SUBS[NSUBS];
    if (imax > spbcon_.NCOLS) {
	spbcon_.NCOLS = imax;
    }
    if (NSUBS > spbcon_.NZMAX) {
	spbcon_.NZMAX = NSUBS;
    }

    switch (TYPE) {
	case 1:  goto L200;
	case 2:  goto L1100;
	case 3:  goto L300;
	case 4:  goto L1200;
    }
L200:
/*        ------------------- */
/*        sparse equation ... */
/*        ------------------- */
    spbcon_.NSEQNS++;
    spbcon_.NZEQNS += NSUBS;
    goto L400;
L300:
/*        --------------------- */
/*        sparse constraint ... */
/*        --------------------- */
    spbcon_.NSCONS++;
    spbcon_.NZCONS += NSUBS;
L400:
/*        ------------------------------------------ */
/*        create linked structure if incoming row is */
/*        either an equation or a constraint ... */
/*        ------------------------------------------ */
    if (NSUBS == 1) {
	return;
    }
    i__1 = NSUBS;
    for (isub = 2; isub <= i__1; isub++) {
	i__ = SUBS[isub];
	if (spbmap_.LLFREE < spbcon_.NCOLS) {
	    T[i__] = -i__;
	}
	link = i__;
	isubm1 = isub - 1;
	i__2 = isubm1;
	for (jsub = 1; jsub <= i__2; jsub++) {
	    j = SUBS[jsub];
L500:
/*                    -------------------------------------- */
/*                    go through the linked list for node i. */
/*                    -------------------------------------- */
	    olink = link;
	    link = T[link];
	    if (link <= 0) {
		goto L600;
	    }
	    nabor = T[link + 1];
	    if (nabor == j) {
		goto L900;
	    }
	    if (nabor < j) {
		goto L500;
	    }
L600:
/*             -----------------------------------------------  */
/*              insert node j into the neighbor list of node i. */
/*              llfree points to the next free space for linked */
/*              lists.                                          */
/*             ----------------------------------------------- */
	    spbcon_.NEDGES++;
	    if (spbusr_.IERRB > 0) {
		goto L800;
	    }
	    if (spbmap_.LLFREE - 2 < spbcon_.NCOLS) {
		goto L700;
	    }
	    T[spbmap_.LLFREE] = j;
	    T[spbmap_.LLFREE - 1] = link;
	    T[olink] = spbmap_.LLFREE - 1;
	    spbmap_.LLFREE += -2;
	    link = olink;    
	    
	    goto L900;
L700:
	    spbusr_.IERRB = 227;
	    spbcon_.STAGE = 10;
	    if (spbusr_.MSGLVB >= 1) {
	       printf("Error 227: Not enough space allocated \n");
	       exit;
	    }
L800:
	    link = olink;
L900:
	    ;
	}
    }
    return;

L1100:
/*        ------------------ */
/*        dense equation ... */
/*        ------------------ */
    spbcon_.NDEQNS++;
    spbcon_.NZEQNS += NSUBS;
    return;

L1200:
/*        -------------------- */
/*        dense constraint ... */
/*        -------------------- */
    spbcon_.NDCONS++;
    spbcon_.NZCONS += NSUBS;
    return;
} /* inxywb */
