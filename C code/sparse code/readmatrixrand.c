/* last revised 17.09.1998 */
/* readmatrixrand.c called by main program */

#include "data_struct.h" 

/**************************************************************************
 * For reading in the matrix A from the given data file.                  *
 **************************************************************************/

MATRIX *readmatrixrand(fp)
FILE *fp;
{
  int Totcrd, Ptrcrd, Indcrd, Valcrd, Nrow, Ncol, Nnzero;
  int *Colptr, *Rowind, *Rowptr, *Colind;
  int PtrNumLines, PtrFieldLen, IndNumLines, IndFieldLen, countElm;
  int LENGTH;
  int i, j, k;
  char strInput[81];
  char Title[72];
  char Key[8];
  char c1;
  char *tempstr1;
  double *Values, *Valuesrow;
  MATRIX *Aa;
  int counter,current;
  double verdi;
  FILE *fp2;
  
/* take the title and key */

  if ((fread (Title,72,1,fp))!=1) { 
    printf(" Error (readmatrix)        :  read Title\n"); 
    exit; 
  } 
  
  if ((fread (Key,8,1,fp))!=1) {
    printf(" Error (readmatrix)        :  read Key\n");
    exit;
  }

/* begin taking the details of the matrix */

  fscanf(fp,"%d%d%d%d%d",&Totcrd,&Nrow,&Ncol,&Nnzero,&Ptrcrd);

/* read until opening of a parantheses */

  c1=')';
  while ( c1 != '(') {
    if ((fread (&c1,1,1,fp))!=1) {
      printf(" Error (readmatrix)        :  read c1\n");
      exit;
    }
  }

  fscanf(fp,"%d",&PtrNumLines);

/* read I */

  if ((fread (&c1,1,1,fp))!=1) {
    printf(" Error (readmatrix)        :  read c1\n");
    exit;
  }

  fscanf(fp,"%d",&PtrFieldLen);

/* read closing of parantheses */

  if ((fread (&c1,1,1,fp))!=1) {
    printf(" Error (readmatrix)        :  read c1\n");
    exit;
  }

  fscanf(fp,"%d",&Indcrd);

/* read until opening of a parantheses */

  while ( c1 != '(') {
    if ((fread (&c1,1,1,fp))!=1) {
      printf(" Error (readmatrix)        :  read c1\n");
      exit;
    }
  }

  fscanf(fp,"%d",&IndNumLines);

/* read I */

  if ((fread (&c1,1,1,fp))!=1) {
    printf(" Error (readmatrix)        :  read c1\n");
    exit;
  }

  fscanf(fp,"%d",&IndFieldLen);

/* read closing of parantheses */

  if ((fread (&c1,1,1,fp))!=1) {
    printf(" Error (readmatrix)        :  read c1\n");
    exit;
  }

  fscanf(fp,"%d",&Valcrd);

/* read until the end of line */

  while ( c1 != '\n' ) {
    if ((fread (&c1,1,1,fp))!=1) {
      printf(" Error (readmatrix)        :  read c1\n");
      exit;
    }
  }

/* begin creating the data structure for the sparse matrix */

/* open the A matrix */

/* next statement is changed on 24.10.1997*/
/* Aa=NEW(MATRIX); */

  Aa=CREATE(1,MATRIX);
  if ( !Aa ) {
    printf("Error in creating matrix.\n");
    exit;
  }

  Aa->rowdefined=TRUE;
  Aa->m=Nrow;
  Aa->n=Ncol;
  Aa->nonzeros=Nnzero;
  Aa->Iind=CREATE(Nnzero+1,int);
  Rowind=Aa->Iind;
  Aa->Jind=CREATE(Ncol+2,int);
  Colptr=Aa->Jind;
  Aa->rowptr=CREATE(Nrow+2,int);
  Rowptr=Aa->rowptr;
  Aa->colind=CREATE(Nnzero+1,int);
  Colind=Aa->colind;
  Aa->Valuecol=CREATE(Nnzero+1,double);
  Values=Aa->Valuecol;
  Aa->Valuerow=CREATE(Nnzero+1,double);
  Valuesrow=Aa->Valuerow;

/* read in the column pointers */

  LENGTH=81;
  tempstr1=CREATE(PtrFieldLen,char);
  countElm=1;
  while (( countElm <= Ncol+2 )&&( LENGTH == 81 )) {
     if ((Ncol+2-countElm)<PtrNumLines) LENGTH=(Ncol+2-countElm)*PtrFieldLen+1;
     if ( fread (strInput,LENGTH,1,fp) != 1) {
        printf(" Error (readmatrix)        :  read Colptr\n");
     }

     i=0;
     while (( i<PtrNumLines )&&( countElm <= Ncol+1 )) {
        k=0;
        for ( j=0;j<PtrFieldLen;j++) {
           if ( strInput[i*PtrFieldLen+j] != ' ' ) {
              tempstr1[k]=strInput[i*PtrFieldLen+j];
              k++;
           }
        }

        Colptr[countElm]=atoi(tempstr1);

        for ( k=0;k<PtrFieldLen;k++) {
           tempstr1[k]=' ';
        }
        i++;
        countElm++;
     }
  }
  Colptr[Ncol+1]=Nnzero+1;
  free(tempstr1);

/* read in the row indexes */

  LENGTH=81;
  tempstr1=CREATE(IndFieldLen,char);
  countElm=1;
  while (( countElm <= Nnzero+1 )&&( LENGTH == 81 )) {
    if (( Nnzero+1-countElm ) < IndNumLines ) {                                  
       LENGTH=(Nnzero+1-countElm)*IndFieldLen+1;
    }
    if ( fread (strInput,LENGTH,1,fp) != 1) {
      printf(" Error (readmatrix)        :  read Rowind\n");
    }

    i=0;
    while (( i<IndNumLines )&&( countElm <= Nnzero )) {
      k=0;
      for ( j=0;j<IndFieldLen;j++) {
        if ( strInput[i*IndFieldLen+j] != ' ' ) {
          tempstr1[k]=strInput[i*IndFieldLen+j];
          k++;
          }
      }

      Rowind[countElm]=atoi(tempstr1);

      i++;
      countElm++;
      for ( k=0;k<IndFieldLen;k++) {
        tempstr1[k]=' ';
      }
    }
  }
  free(tempstr1);

/* initialize the values of nonzeros in each column */
fp2=fopen("aadata","w");
  for ( i=1;i<=Aa->nonzeros;i++ ) {
    if ( i%2 == 0 )
      Aa->Valuecol[i]=mrand();
    if ( i%2 == 1 )
      Aa->Valuecol[i]=-mrand();
    fprintf(fp2,"\n %lf",Aa->Valuecol[i]);
  }
fclose(fp2);
/* form the row order of the matrix */

  counter=1;
  for ( j=1;j<=Aa->m;j++ ) {
    Rowptr[j]=counter;
    current=1;
    for ( i=1;i<=Aa->nonzeros;i++) {
      if ( Aa->Iind[i]==j ) {
        while ( i>Colptr[current+1] )
          current++;
        if ( i==Colptr[current+1] ) current++;
        Colind[counter]=current;
        Valuesrow[counter]=Values[i];
        counter++;
      }
    }
  }
  Rowptr[Aa->m+1]=Aa->nonzeros+1;
  fclose(fp);
  return Aa;
}

