/* last revised 07.09.1998 */
/* readnewmatrix.c called by main program to read the bigger */
/* Harwell matrices */

#include "data_struct.h" 

/**************************************************************************
 * For reading in the matrix A from the given Harwell-Boeing file.        *
 **************************************************************************/

void readnewmatrix(fp,fp1)
FILE *fp,*fp1;
{
  int Totcrd,  /* total number of lines excluding header */ 
      Ptrcrd,  /* number of lines for pointers */
      Indcrd,  /* number of lines for row indices */
      Valcrd,  /* number of lines for right-hand sides */
      Rhscrd,  /* number of lines for right-hand sides */
      Nrow,    /* number of rows */
      Ncol,    /* number of columns */
      Nnzero,  /* number of nonzero elements */
      Neltvl;  /* number of elemental matrix entries(zero if assembled matrix)*/
  int PtrNumLines, /* number of column pointers per line */
      PtrFieldLen, /* field length of pointers */
      IndNumLines, /* number of row indices per line */
      IndFieldLen, /* field |ength of indices */
      ValNumLines, /* number of values of nonzeros per line */
      ValFieldLen, /* field length of values */
      RhsNumLines, /* number of right-hand side values per line */
      RhsFieldLen; /* field length of right-hand side values */ 
  int Nrhs,    /* number of right-hand sides */
      Nrhsix;  /* number of row indices */
  int LENGTH;
  int i,j,k,line;
  int countElm;
  int counter,current;
  char strInput[81];
  char Title[72]; /* title of the matrix */
  char Key[8],    /* key for the matrix */
       Mxtype[3], /* matrix type */
       Ptrfmt[16], /* pointer format */
       Indfmt[16],    /* index format */
       Valfmt[20],    /* value format */
       Rhsfmt[20];    /* right-hand side format */
  char Rhstyp,      /* right-hand side type */
       Guess,       /* G if a starting vector is supplied */
       Xsol;        /* X if an exact solution is supplied */     
  char c1,*c2;
  char *tempstr1;
  int *Rowind, *Colptr, *Rowptr, *Colind;
  double *Values, *Valuesrow;
  
/* read the first line */

  if ((fread (&Title,72,1,fp))!=1) {
    printf(" Error (read matrix)        :  read Title\n");
    exit;
  }
  if ((fread (&Key,8,1,fp))!=1) {
    printf(" Error (read matrix)        :  read Key\n");
    exit;
  }

#ifdef DEBUG
  fprintf(fp1," Debug (read matrix) : Title=[%72s] \n",Title);
  fprintf(fp1," Debug (read matrix) : Key=[%8s] \n",Key);
#endif

/* read the second line */

  fscanf(fp,"%d%d%d%d%d \n",&Totcrd,&Ptrcrd,&Indcrd,&Valcrd,&Rhscrd,&c1);
#ifdef DEBUG
  fprintf(fp1,"%d %d %d %d %d \n",Totcrd,Ptrcrd,Indcrd,Valcrd,Rhscrd);
#endif
   
/* read the third line */
 
  if ((fread (&strInput,81,1,fp))!=1) {
    printf(" Error (read matrix)        :  read third line\n");
    exit;
  } 
  fscanf(fp,"\n",&c1);
  k=0;
  for (i=0;i<3;i++) {
    if (strInput[i] != ' ') {
      Mxtype[k]=strInput[i];
      k++;
    }
  }
#ifdef DEBUG  
  printf(" Debug (read matrix) : Mxtype=[%3s] \n",Mxtype);
#endif
  
  tempstr1=CREATE(13,char);
  k=0;
  for (i=14;i<28;i++) { 
    if (strInput[i] != ' ') {
      tempstr1[k]=strInput[i]; 
      k++;
    }
  }
  Nrow=atoi(tempstr1);
  free(tempstr1);
  tempstr1=CREATE(13,char);
  k=0;
  for (i=28;i<42;i++) {
    if (strInput[i] != ' ') {
      tempstr1[k]=strInput[i];
      k++;
    }
  }
  Ncol=atoi(tempstr1);
  free(tempstr1);
  tempstr1=CREATE(13,char);
  k=0;
  for (i=42;i<56;i++) {
    if (strInput[i] != ' ') {
      tempstr1[k]=strInput[i];
      k++;
    }
  }
  Nnzero=atoi(tempstr1);
  free(tempstr1);
  tempstr1=CREATE(13,char);
  k=0;
  for (i=56;i<70;i++) {
    if (strInput[i] != ' ') {
      tempstr1[k]=strInput[i];
      k++;
    }
  }
  Neltvl=atoi(tempstr1);
  free(tempstr1);

#ifdef DEBUG  
  fprintf(fp1,"%d %d %d %d",Nrow,Ncol,Nnzero,Neltvl);
#endif

/* read the fourth line */

  if ((fread (&strInput,81,1,fp))!=1) {
    printf(" Error (read matrix)        :  read third line\n");
    exit;
  }
  fscanf(fp,"\n",&c1);
  for (i=0;i<16;i++)
    Ptrfmt[i]=strInput[i];
  tempstr1=CREATE(2,char);
  for(i=0;i<2;i++)
    tempstr1[i]=Ptrfmt[i+1];
  PtrNumLines=atoi(tempstr1);

#ifdef DEBUG  
  fprintf(fp1,"%d \n",PtrNumLines);
#endif
  
  free(tempstr1);
  tempstr1=CREATE(2,char);
  tempstr1[0]=Ptrfmt[3];
  tempstr1[1]=Ptrfmt[4];
  sscanf(tempstr1,"%c%d",&c1,&PtrFieldLen);

#ifdef DEBUG  
  fprintf(fp1,"%d \n",PtrFieldLen);
#endif
  
  free(tempstr1);
  for (i=16;i<32;i++)
    Indfmt[i-16]=strInput[i];
  tempstr1=CREATE(2,char);
  for(i=0;i<2;i++)
    tempstr1[i]=Indfmt[i+1];
  IndNumLines=atoi(tempstr1);
  
#ifdef DEBUG  
  fprintf(fp1,"%d \n",IndNumLines);
#endif
  
  free(tempstr1);
  tempstr1=CREATE(2,char);
  tempstr1[0]=Indfmt[3];
  tempstr1[1]=Indfmt[4];
  sscanf(tempstr1,"%c%d",&c1,&IndFieldLen);

#ifdef DEBUG  
  fprintf(fp1,"%d \n",IndFieldLen);
#endif
  
  free(tempstr1);
  for (i=32;i<52;i++)
    Valfmt[i-32]=strInput[i];
  tempstr1=CREATE(3,char);
  tempstr1[0]=Valfmt[3];
  tempstr1[1]=Valfmt[4];
  sscanf(tempstr1,"%c%d",&c1,&ValNumLines);
  free(tempstr1);
#ifdef DEBUG  
  fprintf(fp1,"%d \n",ValNumLines);
#endif
  
  tempstr1=CREATE(2,char);
  for(i=0;i<2;i++)
    tempstr1[i]=Valfmt[i+6];
  ValFieldLen=atoi(tempstr1);
  
#ifdef DEBUG  
  fprintf(fp1,"%d \n",ValFieldLen);
#endif
  
  free(tempstr1);
  for (i=52;i<72;i++)
    Rhsfmt[i-52]=strInput[i];
  tempstr1=CREATE(2,char);
  tempstr1[0]=Rhsfmt[3];
  tempstr1[1]=Rhsfmt[4];
  sscanf(tempstr1,"%c%d",&c1,&RhsNumLines); 
  
#ifdef DEBUG  
  fprintf(fp1,"%d \n",RhsNumLines);
#endif
  
  free(tempstr1);
  tempstr1=CREATE(2,char);
  for(i=0;i<2;i++)
    tempstr1[i]=Rhsfmt[i+6];
  RhsFieldLen=atoi(tempstr1);
  
#ifdef DEBUG  
  fprintf(fp1,"%d \n",RhsFieldLen);
#endif
  
  free(tempstr1);
 
/* read the fifth line */

  if ((fread (&strInput,81,1,fp))!=1) {
    printf(" Error (read matrix)        :  read third line\n");
    exit;
  }

  Rhstyp=strInput[0];
  Guess=strInput[1];
  Xsol=strInput[2];
  tempstr1=CREATE(13,char);
  k=0;
  for (i=14;i<28;i++) {
    if (strInput[i] != ' ') {
      tempstr1[k]=strInput[i];
      k++;
    }
  }
  Nrhs=atoi(tempstr1);
  free(tempstr1);
  tempstr1=CREATE(13,char);
  k=0;
  for (i=28;i<42;i++) {
    if (strInput[i] != ' ') {
      tempstr1[k]=strInput[i];
      k++;
    }
  }
  Nrhsix=atoi(tempstr1);
  free(tempstr1);

#ifdef DEBUG  
  fprintf(fp1,"\n %d %d \n",Nrhs,Nrhsix); 
#endif  

/* begin creating the data structure for the sparse matrix */

#ifdef DEBUG
  fprintf(fp1," Debug (read matrix) : A is created\n");
  fprintf(fp1," Debug (read matrix) : Nrow=[%d]\n",Nrow);
  fprintf(fp1," Debug (read matrix) : Ncol=[%d]\n",Ncol);
#endif

  A=CREATE(1,MATRIX);
  A->rowdefined=TRUE;
  A->m=Nrow;
  A->n=Ncol;
  A->nonzeros=Nnzero;
  A->Iind=CREATE(Nnzero+1,int);
  Rowind=A->Iind;
  A->Jind=CREATE(Ncol+2,int);
  Colptr=A->Jind;
  A->rowptr=CREATE(Nrow+2,int);
  Rowptr=A->rowptr;
  A->colind=CREATE(Nnzero+1,int);
  Colind=A->colind;
  A->Valuecol=CREATE(Nnzero+1,double);
  Values=A->Valuecol;
  A->Valuerow=CREATE(Nnzero+1,double);
  Valuesrow=A->Valuerow;

/* read in the column pointers */

  countElm=1;
  while ( countElm <= Ncol+1 ) {
    fscanf(fp,"%d",&Colptr[countElm]);
    if (( countElm % PtrNumLines) == 0 )
      fscanf(fp,"\n",&c1);

#ifdef DEBUG
  fprintf(fp1," Colptr[%d] :[%d] \n ",countElm,Colptr[countElm]);
#endif
    countElm++;
  }

  fscanf(fp,"\n",&c1);

/* read in the row indices */

  countElm=1;
  while ( countElm <= Nnzero ) {
    fscanf(fp,"%d",&Rowind[countElm]);

#ifdef DEBUG
  fprintf(fp1,"  Rowind[%d] :[%d] \n ",countElm,Rowind[countElm]);
#endif

    if (( countElm % IndNumLines ) == 0 )
      fscanf(fp,"\n",&c1);
    countElm++;
  }
  c1='a';
/* read until the end of line */

  while ( c1 != '\n' ) {
    if ((fread (&c1,1,1,fp))!=1) {
      printf(" Error (read matrix)        :  read c1\n");
      exit;
    }
  }

/*  fscanf(fp,"\n",&c1); */

/* read in the nonzero values */

  countElm=1; 
  while ( countElm <= Nnzero ) {  
    fscanf(fp,"%lf",&Values[countElm]);
    if ( (countElm % ValNumLines) == 0 )
      fscanf(fp,"\n",&c1);

#ifdef DEBUG
  fprintf(fp1," Values[%d] :[%14.10g] \n ",countElm,Values[countElm]);
#endif

    countElm++; 
  }

/* read until the end of line */

  c1='a';
  while ( c1 != '\n' ) {
    if ((fread (&c1,1,1,fp))!=1) {
      printf(" Error (read matrix)        :  read c1\n");
      exit;
    }
  }

  rhs=v_get(A->m);

/* read in the right-hand side vector */

  countElm=0;
  while ( countElm < Nrow ) {
    fscanf(fp,"%lf",&rhs->ve[countElm]);
    countElm++;
    if ( (countElm % RhsNumLines) == 0 )
      fscanf(fp,"\n",&c1);
  }
   
/* form the row order of the matrix */

  counter=1;
  for ( j=1;j<=A->m;j++ ) {
    Rowptr[j]=counter;
    current=1;
    for ( i=1;i<=A->nonzeros;i++) {
      if ( A->Iind[i]==j ) {
        while ( i>Colptr[current+1] )
          current++;
        if ( i==Colptr[current+1] ) current++;
        Colind[counter]=current;
        Valuesrow[counter]=Values[i];
        counter++;
      }
    }
  }
  Rowptr[A->m+1]=A->nonzeros+1;
  fclose(fp);
}

