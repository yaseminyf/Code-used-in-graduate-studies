/* last revised 05.07.1999 */
/* utility file containing subroutines called by the main program */

#include "data_struct.h" 

/**************************************************************************
 * This function is for forming the non-overlapping sub matrices with     *
 * different number of groups, where each group has n/ngroup columns      *
 **************************************************************************/
MATRIX **create_nover_diff(A,sub_ptr,ngroup,zeroptr,nc)
MATRIX *A;
MATRIX **sub_ptr;
int ngroup;
int **zeroptr,nc;
{ static int i,j,k;
  int current,found,count,counter,nclst;
  MATRIX *a;
  int *zerorows;
  double number,num,*ip;
  
  nclst=A->n-((ngroup-1)*nc);
  for ( i=0;i<ngroup;i++ ) {
    a=CREATE(1,MATRIX);
    a->rowdefined=TRUE;
    a->m=A->m;
    if ( i==ngroup-1 ) {
      a->n=nclst;
      a->nonzeros=A->Jind[A->n+1]-A->Jind[i*nc+1];
    }
    else {
      a->n=nc;
      a->nonzeros=A->Jind[(i+1)*nc+1]-A->Jind[i*nc+1];
    }
    
    a->Iind=CREATE(a->nonzeros+1,int);
    a->Jind=CREATE(a->n+2,int);
    a->rowptr=CREATE(a->m+2,int);
    a->colind=CREATE(a->nonzeros+1,int);
    a->Valuerow=CREATE(a->nonzeros+1,double);
    a->Valuecol=CREATE(a->nonzeros+1,double);
    zerorows=CREATE(A->m+1,int);

/* form the column oriented data structure first */

    for ( k=0;k<a->nonzeros;k++ ) {
      a->Iind[k+1]=A->Iind[A->Jind[i*nc+1]+k];
      a->Valuecol[k+1]=A->Valuecol[A->Jind[i*nc+1]+k];
    }     
    a->Jind[1]=1;
    counter=1;
    if ( i == ngroup-1 ) {
      for ( j=2;j<=nclst;j++ ) {
        counter=counter+A->Jind[i*nc+j]-A->Jind[i*nc+j-1];
        a->Jind[j]=counter;
      }
    }
    else {
      for ( j=2;j<=nc;j++ ) {
        counter=counter+A->Jind[i*nc+j]-A->Jind[i*nc+j-1];
        a->Jind[j]=counter;
      }
    }
    
    a->Jind[a->n+1]=a->nonzeros+1;
    
/* form the row oriented data structure */

    count=1;
    for ( j=1;j<=a->m;j++ ) {
      found=0;
      a->rowptr[j]=count;
      current=1;
      for ( k=1;k<=a->nonzeros;k++ ) {
        if ( a->Iind[k]==j ) {
          found=1;
          while ( k >= a->Jind[current+1] )
            current++;
          a->colind[count]=current;
          a->Valuerow[count]=a->Valuecol[k];
          count++;
        }
      }

/* if this is a zero row then assign zero to the rowptr to identify it */
/* and assign 1 to the zerorowptr to identify the zero row */

      if ( found==0 ) {
        a->rowptr[j]=0;
        zerorows[j]=1;
        zerorows[0]+=1;
      }
    }
    a->rowptr[a->m+1]=a->nonzeros+1;
    sub_ptr[i]=a;
    zeroptr[i]=zerorows;
  }
   return(sub_ptr); 
}
/* This function is for creating the hatA matrix called C                     */
/******************************************************************************/
void createC1(hatA,ngroup,zeroptr,indicesptr) 
MATRIX *hatA;
int ngroup,**zeroptr,**indicesptr;
{ static int i,j,k;
  int count,current,size;
  
  size=hatA->m;
  count=0;
  for (i=0;i<ngroup;i++) {           /* the number of nonzeros */
    current=size-zeroptr[i][0];
    count+=current;
    indicesptr[i]=CREATE(current,int);
  }
  hatA->nonzeros=count; 
  hatA->Jind=CREATE(ngroup+2,int);
  hatA->Iind=CREATE(hatA->nonzeros+1,int);
  count=1;
  for (i=0;i<ngroup;i++) { 
    current=0; 
    hatA->Jind[i+1]=count;
    for (j=1;j<=size;j++) {
      if (zeroptr[i][j]==0) {
        hatA->Iind[count]=j;
        indicesptr[i][current]=j;
        current++;
        count++;
      }
    }
  }
  hatA->Jind[hatA->n+1]=count;

/* create the row oriented data structure */

  hatA->rowptr=CREATE(hatA->m+2,int);  
  hatA->colind=CREATE(hatA->nonzeros+1,int);  
  count=1;  
  for ( j=1;j<=hatA->m;j++ ) {  
    hatA->rowptr[j]=count;  
    current=1;  
    for ( k=1;k<=hatA->nonzeros;k++ ) {  
      if ( hatA->Iind[k]==j ) {  
        while ( k >= hatA->Jind[current+1] )  
          current++;  
        hatA->colind[count]=current;  
        count++;  
      }  
    }  
  }  
  hatA->rowptr[hatA->m+1]=hatA->nonzeros+1;      
  hatA->Valuecol=CREATE(hatA->nonzeros+1,double);
}
/* This function is for putting the nonzero values into the hatA matrix C    */
/*****************************************************************************/
void createC222(C,DS,ngroup)
MATRIX *C;
VECTOR **DS;
int ngroup;
{ int count,sind,eind,index;
  static int i,j,k;
  
  count=1;
  for ( j=1; j<=ngroup; j++ ) {
    eind=DS[j]->dim;
    for (i=0;i<eind;i++) {
      C->Valuecol[count]=DS[j]->ve[i];
      count++;
    }
  }
/*  C->Valuerow=CREATE(C->nonzeros+1,double); */
/*  count=1; */
/*   */
/*  for ( j=1;j<=C->m;j++ ) {  */
/*    sind=C->rowptr[j];  */
/*    eind=C->rowptr[j+1]-1;  */
/*    for ( k=sind;k<=eind;k++ ) {  */
/*      index=C->Jind[C->colind[k]];  */
/*      while ( C->Iind[index] != j )  */
/*        index++;  */
/*      C->Valuerow[count]=C->Valuecol[index];  */
/*      count++;  */
/*    }  */
/*  }    */
}
/* This function is for the sequential programs. It is for creating a pointer */
/* to the tildeA_i matrices.                                                  */
/******************************************************************************/
MATRIX **createtildeAptr(tildeAptr,sub_ptr,Apptr,ngroup)
MATRIX **tildeAptr,**sub_ptr;
VECTOR **Apptr;
int ngroup;
{ MATRIX *a;
  static int i,j,k,l;
  int count,current,colcounter,*countptr;
  int sind, eind;
  
  l=sub_ptr[0]->m;
  countptr=CREATE(ngroup,int);
  for (i=0;i<ngroup;i++ ) {
    for (k=0;k<l;k++) {
      if (Apptr[i]->ve[k]!=0.0)
        countptr[i]++;
    }
  }
  for ( i=0;i<ngroup;i++ ) {
    a=CREATE(1,MATRIX);
    a->m=sub_ptr[i]->m; 
    a->n=sub_ptr[i]->n+ngroup-1;
    count=0;
    for (j=0;j<ngroup;j++) {
      if (j!=i) 
        count+=countptr[j];
    }
    a->nonzeros=sub_ptr[i]->nonzeros+count;
    a->Jind=CREATE(a->n+2,int);
    a->Iind=CREATE(a->nonzeros+1,int);
    a->Valuecol=CREATE(a->nonzeros+1,double);
    count=1;
    colcounter=1;
    for ( j=1; j<=ngroup; j++ ) {
      if (j!=(i+1)) {
        a->Jind[colcounter]=count;
        for ( k=1; k<=a->m; k++ ) {
          if (Apptr[j-1]->ve[k-1]!=0.0) {
            a->Iind[count]=k;
            a->Valuecol[count]=Apptr[j-1]->ve[k-1];
            count++;
          }
        }
        colcounter++;
      }
      else {
        for (k=1;k<=sub_ptr[i]->n;k++) { /* for each column of Ai do */
          a->Jind[colcounter]=count;
          sind=sub_ptr[i]->Jind[k];
          eind=sub_ptr[i]->Jind[k+1];
          for (l=sind;l<eind;l++) {
            a->Iind[count]=sub_ptr[i]->Iind[l];
            a->Valuecol[count]=sub_ptr[i]->Valuecol[l];
            count++;
          }
          colcounter++;
        }
      }
    }
    a->Jind[a->n+1]=count;

/* create the row oriented data structure */

    a->rowptr=CREATE(a->m+2,int);
    a->colind=CREATE(a->nonzeros+1,int);
    a->Valuerow=CREATE(a->nonzeros+1,double);
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
    tildeAptr[i]=a;
  }
  return(tildeAptr);
}
/* This function is for putting the nonzero values of the tildeA matrix into  */
/* the matrix structure.                                                      */
/******************************************************************************/
void putvaluesintildeA(tildeAptr,sub_ptr,pvecptr,ngroup)
MATRIX **tildeAptr,**sub_ptr;
VECTOR **pvecptr;
int ngroup;
{ VECTOR **Apptr,*sol;
  static int i,j,k,l;
  int count,current,colcounter;
  MATRIX *AA;
  
  Apptr=CREATE(ngroup,VECTOR *);
  for ( i=0; i<ngroup; i++ ) {
     AA=sub_ptr[i];
     sol=v_get(AA->m);
     sol=m_fulvectormlt(AA,pvecptr[i],sol); 
     Apptr[i]=sol;
  }
  for ( i=0;i<ngroup;i++ ) {
    AA=tildeAptr[i];
    count=1;
    colcounter=1;
    for ( j=1; j<=ngroup; j++ ) {
      if (j!=(i+1)) {
        for ( k=1; k<=AA->m; k++ ) {
          if (Apptr[j-1]->ve[k-1]!=0.0) {
            AA->Valuecol[count]=Apptr[j-1]->ve[k-1];
            count++;
          }
        }
      }            
      else {
        count+=sub_ptr[i]->nonzeros;
/*        for (k=1;k<=sub_ptr[i]->n;k++) { */  /* for each column of Ai do */
/*          for (l=sub_ptr[i]->Jind[k];l<sub_ptr[i]->Jind[k+1];l++) { */
/*            AA->Valuecol[count]=sub_ptr[i]->Valuecol[l]; */
/*            count++; */
/*          } */
/*        } */
      } 
    }
    count=1;
    for ( j=1;j<=AA->m;j++ ) {
      current=1;
      for ( k=1;k<=AA->nonzeros;k++ ) {
        if ( AA->Iind[k]==j ) {
          while ( k >= AA->Jind[current+1] )
            current++;
          AA->Valuerow[count]=AA->Valuecol[k];
          count++;
        }
      }
    }
  }
  for ( i=0; i<ngroup; i++ ) {
    free((double *)(Apptr[i]->ve));
    free((VECTOR *)(Apptr[i]));
  }
  free((VECTOR **)(Apptr)); 
}
/******************************************************************************/
void putvaluesintildeA2(tildeAptr,Apptr,sub_ptr,ngroup)
MATRIX **tildeAptr;
VECTOR **Apptr;
MATRIX **sub_ptr;
int ngroup;
{ static int i,j,k,l;
  int count,current,colcounter;
  MATRIX *AA;
  
  for ( i=0;i<ngroup;i++ ) {
    AA=tildeAptr[i];
    count=1;
    colcounter=1;
    for ( j=1; j<=ngroup; j++ ) {
      if (j!=(i+1)) {
        for ( k=1; k<=AA->m; k++ ) {
          if (Apptr[j-1]->ve[k-1]!=0.0) {
            AA->Valuecol[count]=Apptr[j-1]->ve[k-1];
            count++;
          }
        }
      }            
      else {
        count+=sub_ptr[i]->nonzeros;
/*        for (k=1;k<=sub_ptr[i]->n;k++) {  */ /* for each column of Ai do */
/*          for (l=sub_ptr[i]->Jind[k];l<sub_ptr[i]->Jind[k+1];l++) { */
/*            AA->Valuecol[count]=sub_ptr[i]->Valuecol[l]; */
/*            count++; */
/*          } */
/*        } */
      } 
    }
    count=1;
    for ( j=1;j<=AA->m;j++ ) {
      current=1;
      for ( k=1;k<=AA->nonzeros;k++ ) {
        if ( AA->Iind[k]==j ) {
          while ( k >= AA->Jind[current+1] )
            current++;
          AA->Valuerow[count]=AA->Valuecol[k];
          count++;
        }
      }
    }
  } 
}    
/******************************************************************************/
void putvaluesintildeA3(tildeA,Apptr,A,ngroup,myid)
MATRIX *tildeA;
VECTOR **Apptr;
MATRIX *A;
int ngroup,myid;
{ int j,k;
  int count;
  
  count=1;
  for ( j=1; j<=ngroup; j++ ) {
    if (j!=myid) {
      for ( k=0; k<Apptr[j-1]->dim; k++ ) {
        tildeA->Valuecol[count]=Apptr[j-1]->ve[k];
        count++;
      }
    }            
    else {
      count+=A->nonzeros;
    } 
  }
} 
/* This function is for creating the matrix structure of tildaA_i matrix.     */   
/******************************************************************************/
MATRIX *createtildeA(tildeA,A,Apptr,ngroup,myid,indicesptr)
MATRIX *tildeA,*A;
VECTOR **Apptr;
int ngroup,myid,**indicesptr;
{ static int i,j,k,l;
  int count,current,colcounter;

  tildeA->m=A->m; 
  tildeA->n=A->n+ngroup-1;
  count=0;
  for (j=0;j<ngroup;j++) {
    if (j!=(myid-1))
      count+=Apptr[j]->dim;
  }
  tildeA->nonzeros=count+A->nonzeros;
  tildeA->Jind=CREATE(tildeA->n+2,int);
  tildeA->Iind=CREATE(tildeA->nonzeros+1,int);
  tildeA->Valuecol=CREATE(tildeA->nonzeros+1,double);
  count=1;
  colcounter=1;
  for ( j=1; j<=ngroup; j++ ) {
    if (j!=myid) {
      tildeA->Jind[colcounter]=count;
      for ( k=0; k<Apptr[j-1]->dim; k++ ) {
        tildeA->Iind[count]=indicesptr[j-1][k];
        tildeA->Valuecol[count]=Apptr[j-1]->ve[k];
        count++;
      }
      colcounter++;
    }
    else {
      for (k=1;k<=A->n;k++) { /* for each column of Ai do */
        tildeA->Jind[colcounter]=count;
        for (l=A->Jind[k];l<A->Jind[k+1];l++) {
          tildeA->Iind[count]=A->Iind[l];
          tildeA->Valuecol[count]=A->Valuecol[l];
          count++;
        }
        colcounter++;
      }
    }
  }
  tildeA->Jind[tildeA->n+1]=count;

/* create the row oriented data structure */

  tildeA->rowptr=CREATE(tildeA->m+2,int);
  tildeA->colind=CREATE(tildeA->nonzeros+1,int);
  
  count=1;
  for ( j=1;j<=tildeA->m;j++ ) {
    tildeA->rowptr[j]=count;
    current=1;
    for ( k=1;k<=tildeA->nonzeros;k++ ) {
      if ( tildeA->Iind[k]==j ) {
        while ( k >= tildeA->Jind[current+1] )
          current++;
        tildeA->colind[count]=current;
        count++;
      }
    }
  }
  tildeA->rowptr[tildeA->m+1]=tildeA->nonzeros+1;  
  
  return(tildeA);
}
/******************************************************************************/
/* This function is for fulling in the hatA matrix in the sequential programs */
/******************************************************************************/
void createC22(C,DS,ngroup)
MATRIX *C;
VECTOR **DS;
int ngroup;
{ int count,sind,eind,index;
  static int i,j,k;
  
  count=1;
  for ( j=1; j<=ngroup; j++ ) {
    sind=C->Jind[j];
    eind=C->Jind[j+1]-1;
    for (i=sind;i<=eind;i++) {
      C->Valuecol[count]=DS[j]->ve[C->Iind[i]-1];
      count++;
    }
  }
}