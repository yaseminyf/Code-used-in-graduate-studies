/* last revised 28.06.1998 */
/* rcopyl.c called by ipendb.c */

/* #include "data_struct.h" */

/* **********     rcopyl ..... copy from bottom up     *********** */
/*     purpose - rcopyl copies the n elements from vector a        */
/*        to vector b from n till 1.  this is sometimes useful     */
/*        when a and b are the same vector.  (a and b are arrays   */
/*        of long integers.)                                       */
/*     input parameters -                                          */
/*        n      - size of vector a.                               */
/*        a      - the integer vector.                             */
/*     output parameter -                                          */
/*        b      - the output integer vector.                      */
/* *************************************************************** */

void rcopyl(n, T, a, b)
int n, *T, a, b;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, ip;

/* *************************************************************** */

    /* Function Body */

#ifdef DEBUG    
    fprintf(fp1,"\n inside rcopyl...");
    fprintf(fp1,"\n a: %d, b: %d, n: %d ",a,b,n);
#endif    
    
    if (n <= 0) {
	return;
    }
    ip = n;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; i__++) {
	T[b+ip] = T[a+ip];
	ip--;
    }
    return;

} /* rcopyl */

