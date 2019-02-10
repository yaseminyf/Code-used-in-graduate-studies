/* last revision 28.06.1998 */
/* ipendb.c called by orcolb.c */

#include "data_struct.h" 


/* ************     ipendb ..... end of row input     ************ */
/*     purpose - this routine is used to build the adjacency       */
/*        structure from the linked lists formed by calls to       */
/*        inxywb.  the transformation is done in two steps --      */
/*            linked         lower adjacency          adjacency    */
/*            lists  ------>    structure    -------> structure.   */
/*        the storage requirement for the transformation is only   */
/*        2*ncols + 2*nedges.  however, a larger storage vector    */
/*        will improve the efficiency of the transformation,       */
/*        and it will be most efficient if the storage is not      */
/*        less than 2*ncols + 3*nedges.                            */
/*     updated variable -                                          */
/*        t      - the storage vector.  on input, it contains      */
/*                 the linked lists.  on output, it contains       */
/*                 the adjacency structure.                        */
/*     program subroutines -                                       */
/*        build , fmadjy, rcopyl.                                  */
/* *************************************************************** */

void ipendb(T)
int *T;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int nex2l, nex2s;
    static int i__;
    static int /* roffs,*/ rn;
    static int nelong;
    static int rne; 
    
/* *************************************************************** */

/*        make sure there is enough store for the build-up. */
/*        ------------------------------------------------- */

    /* Function Body */
    
#ifdef DEBUG        
    fprintf(fp1,"\n inside ipendb.."); 
#endif
    
    rn = spbcon_.NCOLS;
    rne = spbcon_.NEDGES;

#ifdef DEBUG        
    fprintf(fp1,"\n spbcon_.NCOLS: %d ",spbcon_.NCOLS);
    fprintf(fp1,"\n spbcon_.NEDGES: %d ",spbcon_.NEDGES);
#endif
    
/* the following statements are commented out by me */    
/*    roffs = spksys_.ratios - (float).01; */
/*    roffl = spksys_.ratiol - (float).01; */

/* the following statements are commented out by me and replaced */
/* with the ones that follow */
/*    spbcon_.nshort = (int) ((rn + roffs) / spksys_1.ratios); */
/*    spbcon_.nlong = (int) ((rn + roffl) / spksys_1.ratiol); */
/*    spbcon_.nsp1 = (int) ((rn + (float)1. + roffs) / spksys_1.ratios); */
/*    spbcon_.nlp1 = (int) ((rn + (float)1. + roffl) / spksys_1.ratiol); */
/*    nelong = (int) ((rne + roffl) / spksys_1.ratiol); */
/*    nex2l = (int) ((rne + rne + roffl) / spksys_1.ratiol); */
/*    nex2s = (int) ((rne + rne + roffs) / spksys_1.ratios); */

    spbcon_.NSHORT = rn; 
    spbcon_.NLONG = rn; 
    spbcon_.NSP1 = rn+1; 
    spbcon_.NLP1 = rn+1; 
    nelong = rne; 
    nex2l = rne + rne; 
    nex2s = rne + rne; 

#ifdef DEBUG        
    fprintf(fp1,"\n spbcon_.NSHORT: %d ",spbcon_.NSHORT);
    fprintf(fp1,"\n spbcon_.NLONG: %d ",spbcon_.NLONG);
    fprintf(fp1,"\n spbcon_.NSP1: %d ",spbcon_.NSP1);
    fprintf(fp1,"\n spbcon_.NLP1: %d ",spbcon_.NLP1);
#endif
    
/* Computing MAX */
    i__1 = spbcon_.NLONG + nelong;
    spbcon_.MXREQD = spbcon_.NLP1 + spbcon_.NLONG + nelong + MAX(i__1,
	    nex2s);
/* the following already exist in the main program */	    
/*    spbusr_1.mcols = spbcon_1.ncols; */
/*    spbusr_1.mseqns = spbcon_1.nseqns; */
/*    spbusr_1.mscons = spbcon_1.nscons; */
/*    spbusr_1.mdeqns = spbcon_1.ndeqns; */
/*    spbusr_1.mdcons = spbcon_1.ndcons; */

/*            ---------------------------------------------------- */
/*            build the lower adj structure from the linked lists. */
/*            ladj points to where the structure starts in t. */
/*            ---------------------------------------------------- */
    xadj = 1;
    deg = xadj + spbcon_.NLP1;
    ladj = deg + spbcon_.NLONG;
    iladj = ladj;
    limit = spbmap_.LLFREE; 
    
#ifdef DEBUG        
    fprintf(fp1,"\n xadj: %d, deg: %d, ladj: %d, iladj: %d, limit: %d ",xadj,deg,ladj,iladj,limit);
#endif
    
    deg--;
    
#ifdef DEBUG        
    fprintf(fp1,"\n right before build..");
#endif
    
    build(spbcon_.NCOLS, spbcon_.NEDGES, T, deg, iladj, limit);
    
#ifdef DEBUG        
    fprintf(fp1,"\n right after build.."); 
#endif
    
/*            ---------------------------------------------- */
/*            copy degree vector deg and lower adj structure */
/*            vector ladj to bottom of the storage vector t. */
/*            and then form the entire adjacency structure.  */
/*            ---------------------------------------------- */
 
    nudeg = spbusr_.MAXSB - spbcon_.NLONG + 1; 
    nuladj = nudeg - nelong; 

#ifdef DEBUG    
    fprintf(fp1,"\n nudeg: %d, nuladj: %d ",nudeg,nuladj);
#endif
    
    nudeg--;
    nuladj--;
    
#ifdef DEBUG        
    fprintf(fp1,"\n before rcopyl.."); 
#endif
    
    rcopyl(spbcon_.NCOLS, T, deg, nudeg);
    
#ifdef DEBUG        
    fprintf(fp1,"\n after rcopyl first.."); 
#endif
    
    ladj--;
    rcopyl(spbcon_.NEDGES, T, ladj, nuladj);
    adjncy = spbcon_.NLP1 + 1;

#ifdef DEBUG        
    fprintf(fp1,"\n adjncy: %d ",adjncy);
#endif
    
    xadj--;
    adjncy--;
    
#ifdef DEBUG        
    fprintf(fp1,"\n before fmadjy...");
    fprintf(fp1,"\n xadj before fmadjy...\n");
    i__1=spbcon_.NCOLS+1;
    for (i__=1;i__<=i__1;i__++) {
      fprintf(fp1,"%d ",T[xadj+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
#endif
    
    fmadjy(spbcon_.NCOLS, T, xadj, adjncy, nuladj, nudeg);
    
#ifdef DEBUG        
    fprintf(fp1,"\n after fmadjy.."); 
    i__1=nex2s;
    fprintf(fp1,"\n adjncy after fmadjy...\n");
    for (i__=1;i__<=i__1;i__++) {
      fprintf(fp1,"%d ",T[adjncy+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n xadj after fmadjy...\n");
    i__1=spbcon_.NCOLS+1;
    for (i__=1;i__<=i__1;i__++) {
      fprintf(fp1,"%d ",T[xadj+i__]);
      if ((i__%10)==0) fprintf(fp1,"\n");
    }
#endif
    
    spbcon_.MXUSED = spbcon_.NLP1 + nex2s; 

    return;

} /* ipendb */

