/* last revised 08.03.1999 */
/* multiply.c called by the main program */

#include "data_struct.h" 

/**************************************************************************
* This function is for multiplying the transpose of a packed matrix with a*
* full vector rrvec and storing the result in the vector rvec which is    *
* initialized to zeros before the beginning of the function.              *
***************************************************************************/
double *multAtransr(A,rrvec,rvec)
MATRIX *A;
double *rrvec,*rvec;
{ static int i,j;
  for ( i=0;i<=A->n;i++ )
    rvec[i]=0.0;
  for ( i=1;i<=A->n;i++ ) {
    for ( j=A->Jind[i];j<=A->Jind[i+1]-1;j++ ) {
      if ( (A->Valuecol[j]!=0.0)&(rrvec[A->Iind[j]]!=0.0)) {
        rvec[i]=rvec[i]+A->Valuecol[j]*rrvec[A->Iind[j]];
      }
    }
  }
  return(rvec);
}
/**************************************************************************
* This function is for multiplying the transpose of a packed matrix with a*
* meschach vector residual and storing the result in the meschach vector  *
* rvec2.                                                                  * 
***************************************************************************/
VECTOR *multAtransr2(A,residual,rvec2)
MATRIX *A;
VECTOR *residual,*rvec2;
{ static int i,j;
  int sind, eind;
  for ( i=1;i<=A->n;i++ ) {
    sind=A->Jind[i];
    eind=A->Jind[i+1]-1;
    for ( j=sind;j<=eind;j++ ) {
/*      if ( (A->Valuecol[j]!=0.0)&(residual->ve[A->Iind[j]-1]!=0.0)) { */
        rvec2->ve[i-1]+=A->Valuecol[j]*residual->ve[A->Iind[j]-1];
/*      } */
    }
  }
  return(rvec2);
}
/**************************************************************************
* This function is for multiplying the transpose of a packed matrix with  *
* the negative of a                                                       *
* meschach vector residual and storing the result in the meschach vector  *
* s.                                                                      * 
***************************************************************************/
VECTOR *multAtransr3(A,residual,s)
MATRIX *A;
VECTOR *residual,*s;
{ static int i,j;
  int sind, eind;
  
  for ( i=1;i<=A->n;i++ ) {
    sind=A->Jind[i];
    eind=A->Jind[i+1]-1;
    for ( j=sind;j<=eind;j++ ) {
      s->ve[i-1]-=A->Valuecol[j]*residual->ve[A->Iind[j]-1];
    }
  } 
  return(s);
}
/**************************************************************************
* This function is for multiplying a packed matrix with a full vector s   *
* and storing the result in the vector sol which is initialized to zeros  *
* before the beginning of the function and it has the size of A->m.       *
***************************************************************************/
VECTOR *m_fulvectormlt(A,s,sol)
MATRIX *A;
VECTOR *s,*sol;  
{ static int j,l;
/*  int numnonzero,count; */
  int sind, eind;
  
/*  for ( i=1;i<=A->m;i++ ) { */
/*    if ( A->rowptr[i]!=0 ) { */
/*      count=1; */
/*      while ( A->rowptr[i+count] == 0 ) */
/*        count++; */
/*      numnonzero=A->rowptr[i+count]-A->rowptr[i]; */
/*      for ( j=A->rowptr[i];j<numnonzero+A->rowptr[i];j++ ) { */
/*        sol->ve[i-1]=sol->ve[i-1]+A->Valuerow[j]*s->ve[A->colind[j]-1]; */
/*      } */
/*    } */
/*  } */

/* new column oriented approach */

  for ( j=1;j<=A->n;j++ ) {
    sind=A->Jind[j];
    eind=A->Jind[j+1]-1;
    for ( l=sind;l<=eind;l++ ) {
      sol->ve[A->Iind[l]-1]+=A->Valuecol[l]*s->ve[j-1];
    }
  }
  return(sol);
}
/**************************************************************************
* This function is for multiplying the transpose of a packed matrix with  *
* itself. It stores the nonzeros of the symmetric triangular part in lnz, *
* the diagonal entries in diag. Used for submatrices.                     * 
***************************************************************************/
void sym_mult(A,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values) 
MATRIX *A; 
double *lnzv,*ddiag;
int *invpv,*xnzsubv,*nzsubv,*xrnzv; 
int *isub, *jsub;
double *values;

{ register int i,j,k,l; 
  double temp; 
  int count; 
  int check1,check2,check3;
  int sind,eind;
  int i__,j__;
  int kk, ksub, kstop, kstrt;
  
  count=0;
  for (j=1; j<=A->n; j++) { 
    for (i=j; i<=A->n; i++) { 
      temp=0.0; 
      for (k=1; k<=A->m; k++) {      
        if (A->rowptr[k] != 0 ) {
        
/* multiply A[k][j] with A[k][i] if they are both nonzero */ 

/*        check1=0; */
/*        sind=A->Jind[j]; */
/*        eind=A->Jind[j+1]-1; */
/*        for (l=sind;l<=eind;l++) { */
/*          if ( A->Iind[l]==k ) { */
/*            check1=l; */
/*          } */
/*        } */
/*        check2=0; */
/*        sind=A->Jind[i]; */
/*        eind=A->Jind[i+1]-1; */
/*        for (l=sind;l<=eind;l++) { */
/*          if ( A->Iind[l]==k ) { */
/*            check2=l; */
/*          } */
/*        } */
          check1=0;
          check2=0;
          sind=A->rowptr[k];
          check3=1;
          while (A->rowptr[k+check3] == 0 )
            check3++;
          eind=A->rowptr[k+check3]-1;
          for (l=sind;l<=eind;l++) {
            if (A->colind[l]==j)
              check1=l;
            if (A->colind[l]==i)
              check2=l;
          }
        
          
          if ((check1!=0)&&(check2!=0)) {
/*          temp+=A->Valuecol[check1]*A->Valuecol[check2]; */
            temp+=A->Valuerow[check1]*A->Valuerow[check2];
             
          }
        }
      }

      if (temp !=0.0) {         

/* put the indices and the value in arrays for further permutation */

        count++;
        isub[count]=i;
        jsub[count]=j;
        values[count]=temp;
      }                      
    }
  }
  
/* put this matrix into the compressed data structure arrays */

  for (k=1;k<=count;k++) {
    adaij5(isub[k],jsub[k],values[k],invpv,ddiag,xrnzv,lnzv,xnzsubv,nzsubv);
  }
    
  return;
}
/**************************************************************************
* This function is for multiplying the transpose of a packed matrix with  *
* itself. It stores the nonzeros of the symmetric triangular part in lnz, *
* the diagonal entries in diag. Used for hatA matrix.                     * 
***************************************************************************/
void sym_mult2(A,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values) 
MATRIX *A; 
double *lnzv,*ddiag;
int *invpv,*xnzsubv,*nzsubv,*xrnzv; 
int *isub, *jsub;
double *values;

{ register int i,j,k,l; 
  double temp; 
  int count; 
  int check1,check2,check3;
  int sind,eind;
  int i__,j__;
  int kk, ksub, kstop, kstrt;
  
  count=0;
  for (j=1; j<=A->n; j++) { 
    for (i=j; i<=A->n; i++) { 
      temp=0.0; 
      for (k=1; k<=A->m; k++) {      
        
/* multiply A[k][j] with A[k][i] if they are both nonzero */ 

        check1=0;
        check2=0;
        sind=A->rowptr[k];
        eind=A->rowptr[k+1]-1;
        for (l=sind;l<=eind;l++) {
          if (A->colind[l]==j)
            check1=l;
          if (A->colind[l]==i)
            check2=l;
        }
          
        if ((check1!=0)&&(check2!=0)) {
          temp+=A->Valuerow[check1]*A->Valuerow[check2];
        }
      }

      if (temp !=0.0) {          

/* put the indices and the value in arrays for further permutation */

        count++;
        isub[count]=i;
        jsub[count]=j;
        values[count]=temp;
      }                       
    }
  }
  
/* put this matrix into the compressed data structure arrays */

  for (k=1;k<=count;k++) {
    adaij5(isub[k],jsub[k],values[k],invpv,ddiag,xrnzv,lnzv,xnzsubv,nzsubv);
  }
    
  return;
}
/**************************************************************************
* This function is for multiplying the transpose of a packed matrix with  *
* itself. It stores the nonzeros of the symmetric triangular part in lnz, *
* the diagonal entries in diag. Used for hatA matrix.                     * 
***************************************************************************/
void sym_mult3(A,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values) 
MATRIX *A; 
double *lnzv,*ddiag;
int *invpv,*xnzsubv,*nzsubv,*xrnzv; 
int *isub, *jsub;
double *values;

{ register int i,j,ii,jj,k; 
  double temp; 
  int count; 
  int sindi,eindi,sindj,eindj,indexi,indexj;

  count=0;
  for (i=1; i<=A->n; i++) { 
    sindi=A->Jind[i];
    eindi=A->Jind[i+1]-1;
    for (j=i; j<=A->n; j++) { 
      temp=0.0; 
      sindj=A->Jind[j];
      eindj=A->Jind[j+1]-1;
      indexi=sindi;
      indexj=sindj;
      while ((indexi<=eindi)&&(indexj<=eindj)) {
        ii=A->Iind[indexi];
        jj=A->Iind[indexj];

/* multiply A[k][j] with A[k][i] if they are both nonzero */ 

        if (ii==jj) {
          temp+=A->Valuecol[indexi]*A->Valuecol[indexj];
          indexi++;
          indexj++;
        }
        else if (ii>jj) {
          indexj++;
        }
        else {
          indexi++;
        }
      }

      if (temp !=0.0) {          

/* put the indices and the value in arrays for further permutation */

        count++;
        isub[count]=i;
        jsub[count]=j;
        values[count]=temp;
      }                       
    }
  }
  
/* put this matrix into the compressed data structure arrays */

  for (k=1;k<=count;k++) {
    adaij5(isub[k],jsub[k],values[k],invpv,ddiag,xrnzv,lnzv,xnzsubv,nzsubv);
  }
}
/**************************************************************************
* This function is for multiplying a packed matrix with a full matrix.    *
* The result is a packed matrix created from scratch and stored           *
* appropriately. A:m*n, D:n*ng, C:m*ng. C is created before sending in.   * 
***************************************************************************/
void multAD(A,D,C,ngroup)
MATRIX *A;
VECTOR **D;
MATRIX *C;
int ngroup;
{ static int i,j,k;
  int count,current;
  double *curcol;
  double **curcolptr;
  curcolptr=CREATE(ngroup,double *);
  --curcolptr;
  count=1;
  for ( j=1; j<=ngroup; j++ ) {
    curcol=CREATE(A->m,double);
    --curcol;
    for ( i=1; i<=A->m; i++ )
      curcol[i]=0.0;
    current=1;
    for ( k=1; k<=A->n; k++ ) {
      for ( i=1; i<=A->m; i++ ) {
        while ( (A->Iind[current]<i) && (current<A->Jind[k+1]))
          current++;
        if ((current<A->Jind[k+1]) && (A->Iind[current]==i)) {
          curcol[i]+=A->Valuecol[current]*D[j]->ve[k-1];
        }
      }
    }
    curcolptr[j]=curcol;     
  }
  for ( j=1; j<=ngroup; j++ ) {
    for ( i=1; i<=A->m; i++ ) {
      if (curcol[i]!=0.0) 
        count++;
    }
  }
  C->nonzeros=count;
  C->Jind=CREATE(ngroup+2,int);
  C->Iind=CREATE(C->nonzeros+1,int);
  C->Valuecol=CREATE(C->nonzeros+1,double);
  count=1;
  for ( j=1; j<=ngroup; j++ ) {
    C->Jind[j]=count;
    for ( i=1; i<=A->m; i++ ) {
      curcol=curcolptr[j];
      if (curcol[i]!=0.0) {
        C->Iind[count]=i;
        C->Valuecol[count]=curcol[i];
        count++;
      }
    }
  }
  C->Jind[C->n+1]=count;

/* create the row oriented data structure */

  C->rowptr=CREATE(C->m+2,int);
  C->colind=CREATE(C->nonzeros+1,int);
  C->Valuerow=CREATE(C->nonzeros+1,double);
  count=1;
  for ( j=1;j<=C->m;j++ ) {
    C->rowptr[j]=count;
    current=1;
    for ( k=1;k<=C->nonzeros;k++ ) {
      if ( C->Iind[k]==j ) {
        while ( k >= C->Jind[current+1] )
          current++;
        C->colind[count]=current;
        C->Valuerow[count]=C->Valuecol[k];
        count++;
      }
    }
  }
  C->rowptr[C->m+1]=C->nonzeros+1;  
}
