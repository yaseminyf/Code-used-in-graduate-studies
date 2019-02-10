/* last revised 29.06.1998 */
/* copysi.c called by orcolb.c */

/* #include "data_struct.h" */

/* **********     copysi ..... copy integer vector     *********** */
/*     purpose - this routine copies the n integer elements from   */
/*        the vector a to b.  (both a and b are arrays of short    */
/*        integers.)                                               */
/*     input parameters -                                          */
/*        n      - size of vector a.                               */
/*        a      - the integer vector.                             */
/*     output parameter -                                          */
/*        b      - the output integer vector.                      */
/* *************************************************************** */

void copysi(n, T, a, b)
int n, *T, a, b;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;

/* *************************************************************** */

    /* Function Body */
    
#ifdef DEBUG        
    fprintf(fp1,"\n inside copysi...");
    fprintf(fp1,"\n a: %d, b: %d, n: %d ",a,b,n);
#endif

    if (n <= 0) {
	return;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; i__++) {
	T[b+i__] = T[a+i__];
    }
    return;

} /* copysi */

