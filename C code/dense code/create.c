/* utility file containing "create" routines for dense matrices */
/* written: 25.08.1999 */
/* last revised: 14.09.1999 */

#include "datastruct_dense.h"

FMATRIX **create_dense_blocks(A,sub_ptr,nc,ngroup)
FMATRIX *A,**sub_ptr;
int nc,ngroup;
{ static int i,j,k;
  FMATRIX *AA;
  int nclst,numnz,sind,eind;
  
  nclst=A->n-((ngroup-1)*nc);
  for (i=0;i<ngroup;i++) {
    AA=CREATE(1,FMATRIX);
    AA->m=A->m;
    if (i==ngroup-1)
      AA->n=nclst;
    else
      AA->n=nc;
    numnz=AA->n*AA->m;
    AA->Value=CREATE(numnz,double);
    sind=i*nc*AA->m;
    eind=sind+numnz;
    k=0;
    for (j=sind;j<eind;j++) {
      AA->Value[k]=A->Value[j];
      k++;
    }
    sub_ptr[i]=AA;
  }
  return(sub_ptr);
}      
/******************************************************************************/
FMATRIX *createtildeA(tildeA,A,Apptr,ngroup,i)
FMATRIX *tildeA,*A;
VECTOR **Apptr;
int ngroup,i;
{ static int j,k,l;
  int count,current;
  
  tildeA->m=A->m; 
  tildeA->n=A->n+ngroup-1;
  count=tildeA->m*tildeA->n;
  tildeA->Value=CREATE(count,double);
  k=0;
  for (j=0;j<ngroup;j++) {
    if (j!=i) {
      for (l=0;l<tildeA->m;l++) {
        tildeA->Value[k]=Apptr[j]->ve[l];
        k++;
      }
    }
    else {
      current=A->m*A->n;
      for (l=0;l<current;l++) {
        tildeA->Value[k]=A->Value[l];
        k++;
      }
    }
  }
  return(tildeA);
}
/******************************************************************************/
void putvaluesintildeA(tildeA,Apptr,AA,ngroup,i)
FMATRIX *tildeA;
VECTOR **Apptr;
FMATRIX *AA;
int ngroup,i;
{ static int j,k,l;
  int count;
  
  count=tildeA->m;
  k=0;
  for (j=0;j<ngroup;j++) {
    if (j!=i) {
      for (l=0;l<count;l++) {
        tildeA->Value[k]=Apptr[j]->ve[l];
        k++;
      }
    }
    else {
      l=AA->m*AA->n;
      k+=l;
    }
  }
}
