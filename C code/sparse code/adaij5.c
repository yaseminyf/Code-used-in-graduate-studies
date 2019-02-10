/* last updated 03.09.1998 */
/* adaij5.c utility routine */

/* #include "data_struct.h" */

/* *********     ADAIJ5 ..... ADD ENTRY INTO MATRIX     ********** */
/* *************************************************************** */
/* *************************************************************** */
/*     PURPOSE - THIS ROUTINE ADDS A NUMBER INTO THE (I,J)-TH      */
/*        POSITION OF A MATRIX STORED IN A COMPRESSED SUBSCRIPT    */
/*        DATA STRUCTURE.                                          */
/*     INPUT PARAMETERS -                                          */
/*        (ISUB,JSUB) - SUBSCRIPTS OF THE NUMBER TO BE ADDED.      */
/*                 ASSUMPTION - ISUB .GE. JSUB.                    */
/*        VALUE  - VALUE OF THE NUMBER TO BE ADDED.                */
/*        INVP   - INVP(I) IS THE NEW POSITION OF THE VARIABLE     */
/*                 WHOSE ORIGINAL NUMBER IS I.                     */
/*     UPDATED PARAMETERS -                                        */
/*        DIAG   - ARRAY CONTAINING THE DIAGONAL ELEMENTS OF THE   */
/*                 COEFFICIENT MATRIX.                             */
/*        (XLNZ,LNZ) - LEVEL STRUCTURE CONTAINING THE NONZERO      */
/*                 ELEMENTS OF THE COLUMNS OF L BELOW THE          */
/*                 DIAGONAL.                                       */
/*        (XNZSUB,NZSUB) - COMPRESSED SUBSCRIPT PAIR FOR THE       */
/*                 COLUMNS OF L.                                   */
/* *************************************************************** */

void adaij5(isub, jsub, value, invpv, ddiag, xrnzv, lnzv, xnzsubv, nzsubv)
int isub, jsub;
double value;
int *invpv;
double *ddiag;
int *xrnzv;
double *lnzv;
int *xnzsubv, *nzsubv;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int ksub, i__, j, k, itemp, kstop, kstrt;

/* *************************************************************** */

    /* Function Body */
    i__ = invpv[isub];
    j = invpv[jsub];
    if (i__ == j) {
	goto L400;
    }
    if (i__ > j) {
	goto L100;
    }
    itemp = i__;
    i__ = j;
    j = itemp;
L100:
    ksub = xnzsubv[j];
    kstrt = xrnzv[j];
    kstop = xrnzv[j + 1] - 1;
    if (kstop < kstrt) {
	goto L500;
    }
/*                ------------------------------------------ */
/*                THE COMPONENT LIES IN THE LOWER TRIANGULAR */
/*                PORTION.                                   */
/*                ------------------------------------------ */
    i__1 = kstop;
    for (k = kstrt; k <= i__1; k++) {
	if (i__ == nzsubv[ksub]) {
	    goto L300;
	}
	ksub++;
    }
    goto L500;
L300:
    lnzv[k] = value;
    return;
L400:
/*        ------------------------------------------------- */
/*        THE COMPONENT LIES ON THE DIAGONAL OF THE MATRIX. */
/*        ------------------------------------------------- */
    ddiag[i__] = value;
    return;

L500:

    return;

} /* adaij5 */

