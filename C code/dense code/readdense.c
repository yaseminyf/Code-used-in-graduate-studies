/* readdense.c */
/* written 25.08.1999 */
/* last revised 25.08.1999 */
/* this function is for reading in a full matrix to be used in dense programs */

#include "datastruct_dense.h"

void readdense(fp)
FILE *fp;
{
  static int i;
  int m,n,numnz;   /* number of rows and columns and nonzeros */
  double verdi;
  
  fscanf(fp,"%d\n",&m);
  fscanf(fp,"%d\n",&n);
  numnz=m*n;
  
/* open the A matrix */

  A=CREATE(1,FMATRIX);
  if ( !A ) {
    printf("Error in creating matrix.\n");
    exit(0);
  }  
  A->m=m;
  A->n=n;
  A->Value=CREATE(numnz,double);
  
/* read in the nonzero values in column oriented order */

  for (i=0;i<numnz;i++) {
    fscanf(fp,"%le\n",&verdi);
    A->Value[i]=verdi;
  }
  fclose(fp);
}  
