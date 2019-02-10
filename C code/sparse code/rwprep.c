/* last revised 29.06.1998 */
/* rwprep.c called by genls3.c */

/* #include "data_struct.h" */

/* *****     rwprep ..... preparation for row transform     ****** */
/*     purpose - this routine prepares the given row for           */
/*        transformation.                                          */
/*     input paramters -                                           */
/*        NSUBS  - number of nonzeros in the given row.            */
/*        SUBS   - column subscripts of nonzeros in given row.     */
/*        VALUES - numerical values of nonzeros in given row.      */
/*        invp   - inverse permutation.                            */
/*     updated parameters -                                        */
/*        nusubs - permuted column subscripts of nonzeros in       */
/*                 given row.                                      */
/*        nuvals - current row in full representation.             */
/* *************************************************************** */

void rwprep(NSUBS, SUBS, VALUES, invp, nusubs, nuvals, T)
int NSUBS, *SUBS;
double *VALUES;
int invp, *nusubs;
double *nuvals;
int *T;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, k, ip;
/*    extern void sorts1();  */

/* *************************************************************** */

    /* Function Body */
    
#ifdef DEBUG        
    fprintf(fp1,"\n inside rwprep...");
#endif
    
    if (NSUBS <= 0) {
	return;
    }

    i__1 = NSUBS;
    for (i__ = 1; i__ <= i__1; i__++) {
	k = SUBS[i__];

#ifdef DEBUG    	
	fprintf(fp1,"\n k: %d ",k);
#endif
	
/*            -------------------------------------------- */
/*            determine the column permutation and put the */
/*            nonzero in the right position ...            */
/*            -------------------------------------------- */
	ip = T[invp+k];
	
#ifdef DEBUG    	
	fprintf(fp1,"\n ip: %d ",ip);
#endif
	
	nusubs[i__] = ip;
	nuvals[ip] = VALUES[i__];
	
#ifdef DEBUG    	
	fprintf(fp1,"\n nuvals[%d]:%16.14e ",ip,nuvals[ip]);
#endif
	
    }
/*        --------------------------------------------- */
/*        sort the nonzeros in ascending order of their */
/*        column subscripts ...                         */
/*        --------------------------------------------- */

#ifdef DEBUG    
    fprintf(fp1,"\n before sorts1...");
#endif
    
    sorts1(NSUBS, nusubs);

#ifdef DEBUG    
    fprintf(fp1,"\n after sorts1...");
#endif
    
    return;

} /* rwprep */

