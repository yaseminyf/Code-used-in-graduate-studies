/* last revised 13.04.1999 */
/* readbigmatrix.c called by main program */

#include "data_struct.h" 

/**************************************************************************
 * For reading in the matrix A from the given data file. This matrix is   *
 * created by matlab and saved in column order. First are the column      *
 * indices, then row indices, then the values of nonzeros in column       *
 * oriented structure.                                                    *
 **************************************************************************/

void readbigmatrix(fp)
FILE *fp;
{
  int *Colptr, *Rowind, *Rowptr, *Colind;
  static int i, j, k;
  int m,n,numnz;  /* number of rows and columns and nonzeros */
  double verdi,*Values,*Valuesrow;
  int count,current;
  
  fscanf(fp,"%d\n",&m);
  fscanf(fp,"%d\n",&n);
  fscanf(fp,"%d\n",&numnz);

/* open the A matrix */

  A=CREATE(1,MATRIX);
  if ( !A ) {
    printf("Error in creating matrix.\n");
    exit;
  }

  A->rowdefined=TRUE;
  A->m=m;
  A->n=n;
  A->nonzeros=numnz;
  A->Iind=CREATE(numnz+1,int);
  Rowind=A->Iind;
  A->Jind=CREATE(n+2,int);
  Colptr=A->Jind;
  A->rowptr=CREATE(m+2,int);
  Rowptr=A->rowptr;
  A->colind=CREATE(numnz+1,int);
  Colind=A->colind;
  A->Valuecol=CREATE(numnz+1,double);
  Values=A->Valuecol;
  A->Valuerow=CREATE(numnz+1,double);
  Valuesrow=A->Valuerow;
  
/* read in the column pointers */

  count=n+1;
  for (i=1;i<=count;i++) {
    fscanf(fp,"%d\n",&k);
    Colptr[i]=k;
  }
  
/* read in the row indices */

  for (i=1;i<=numnz;i++) {
    fscanf(fp,"%d\n",&k);
    Rowind[i]=k;
  }   
  
/* read in the nonzero values in column oriented structure */

  for (i=1;i<=numnz;i++) {
    fscanf(fp,"%le\n",&verdi);
    Values[i]=verdi;
  }
  
/* form the row order of the matrix */

  count=1;
  for ( j=1;j<=A->m;j++ ) {
    Rowptr[j]=count;
    current=1;
    for ( i=1;i<=A->nonzeros;i++) {
      if ( A->Iind[i]==j ) {
        while ( i>=Colptr[current+1] )
          current++;
        Colind[count]=current;
        Valuesrow[count]=Values[i];
        count++;
      }
    }
  }
  Rowptr[A->m+1]=A->nonzeros+1;
  fclose(fp);
}
