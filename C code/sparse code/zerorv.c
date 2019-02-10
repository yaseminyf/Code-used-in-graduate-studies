/* last revision 04.02.1998 */
/* zerorv.c called by lsqslvsn3.c */

/* #include "data_struct.h" */

/* ***********     zerorv ..... zero a real vector     *********** */
/*     purpose - this routine is used to zero out a real vector.   */
/*     input parameter -                                           */
/*        n      - the size of the real vector.                    */
/*     output parameter -                                          */
/*        rvectr - the real vector.                                */
/* *************************************************************** */

void zerorv(n, rvectr)
int n;
double *rvectr;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;

/* *************************************************************** */

    /* Function Body */
    if (n <= 0) {
	return;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rvectr[i__] = 0.;
    }
    return;

} /* zerorv */

