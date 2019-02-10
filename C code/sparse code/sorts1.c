/* last revised 29.06.1998 */
/* sorts1.c called by rwprep.c */

/* #include "data_struct.h" */

/* *********     sorts1 ..... linear insertion sort     ********** */
/*     purpose - sorts1 uses linear insertion to sort the          */
/*        given array of short integers into increasing order.     */
/*     input parameter -                                           */
/*        NSUBS     - the size of integer array.                      */
/*     updated parameter -                                         */
/*        array  - the integer vector, which on output will be     */
/*                 in increasing order.                            */
/* *************************************************************** */

void sorts1(NSUBS, array)
int NSUBS, *array;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int node, k, l;

/* *************************************************************** */

    /* Function Body */

#ifdef DEBUG    
    fprintf(fp1,"\n inside sorts1...");
#endif
    
    if (NSUBS <= 1) {
	return;
    }
    i__1 = NSUBS;
    for (k = 2; k <= i__1; k++) {
	node = array[k];
	l = k - 1;

L100:
	if (l < 1) {
	    goto L200;
	}
	if (array[l] <= node) {
	    goto L200;
	}
	array[l + 1] = array[l];
	l--;
	goto L100;

L200:
	array[l + 1] = node;
    }
    return;

} /* sorts1 */

