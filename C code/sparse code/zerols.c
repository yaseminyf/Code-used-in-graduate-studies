/* last revision 29.06.1998 */
/* zerols.c called by lsqslvsn3.c */

/* #include "data_struct.h" */

/* *********     zerols ..... zero a level structure     ********* */
/*     purpose - this routine is used to zero out a given          */
/*        level structure.                                         */
/*     input parameters -                                          */
/*        nlvl   - number of levels in the level structure.        */
/*        xls    - index vector to the level structure vector.     */
/*     output parameter -                                          */
/*        ls     - storage for level structure.                    */
/* *************************************************************** */

void zerols(nlvl, T, xls, ls)
int nlvl, *T, xls;
double *ls;
{
    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int i__, l, lstop, lstrt;

/* *************************************************************** */

    /* Function Body */

#ifdef DEBUG    
    fprintf(fp1,"\n inside zerols...");
#endif
    if (nlvl <= 0) {
	return;
    }
    i__1 = nlvl;
    for (l = 1; l <= i__1; l++) {
	lstrt = T[xls+l];
	lstop = T[xls+l + 1] - 1;
	if (lstrt > lstop) {
	    goto L200;
	}
	i__2 = lstop;
	for (i__ = lstrt; i__ <= i__2; i__++) {
	    ls[i__] = 0.;
	}
L200:
	;
    }
    return;

} /* zerols */

