/* last revised 21.09.1998 */
/* utility routine */

#include "data_struct.h"

/******************************************************************************/
/* This routine takes in a whole matrix stored in the Harwell Boeing structure*/
/* and a column permutation p. Returns the permuted matrix structure.         */
/******************************************************************************/

void permuteA(A,p)
MATRIX *A;
int *p;
{ static int i,j,k;
  int count,current;
  MATRIX *a;
  
  a=CREATE(1,MATRIX);
  a->rowdefined=TRUE;
  a->m=A->m;
  a->n=A->n;
  a->nonzeros=A->nonzeros;
  a->Iind=CREATE(a->nonzeros+1,int);
  a->Jind=CREATE(a->n+2,int);
  a->rowptr=CREATE(a->m+2,int);
  a->colind=CREATE(a->nonzeros+1,int);
  a->Valuerow=CREATE(a->nonzeros+1,double);
  a->Valuecol=CREATE(a->nonzeros+1,double);

/* for each column do */
  
  count=1;
  for (i=1;i<=a->n;i++) {
    a->Jind[i]=count;
    k=p[i];
    
/* copy the kth column to ith column */
 
    for ( j=A->Jind[k];j<A->Jind[k+1];j++ ) {
      a->Iind[count]=A->Iind[j];
      a->Valuecol[count]=A->Valuecol[j];
      count++;
    }
  }
  a->Jind[a->n+1]=a->nonzeros+1;
  
/* form the row oriented structure */

  count=1;
  for ( j=1;j<=a->m;j++ ) {
    a->rowptr[j]=count;
    current=1;
    for ( k=1;k<=a->nonzeros;k++ ) {
      if ( a->Iind[k]==j ) {
        while ( k >= a->Jind[current+1] )
          current++;
        a->colind[count]=current;
        a->Valuerow[count]=a->Valuecol[k];
        count++;
      }
    }
  }
  a->rowptr[a->m+1]=a->nonzeros+1;  
  free(A->Iind);
  A->Iind=a->Iind;
  free(A->Jind);
  A->Jind=a->Jind;
  free(A->rowptr);
  A->rowptr=a->rowptr;
  free(A->colind);
  A->colind=a->colind;
  free(A->Valuerow);
  A->Valuerow=a->Valuerow;
  free(A->Valuecol);
  A->Valuecol=a->Valuecol;
}
