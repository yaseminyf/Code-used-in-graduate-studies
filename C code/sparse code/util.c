/* last revised 27.10.1998 */
/* utility file containing subroutines called by the main program */

#include "data_struct.h" 

/**************************************************************************
 * Dot product function.                                                  *
 **************************************************************************/
double dotproduct(sol,residual,check)
VECTOR *sol,*residual;
double check;
{ 
  static int i;
  for ( i=0;i<residual->dim;i++ ) 
    check+=sol->ve[i]*residual->ve[i];
  return check;
}

/**************************************************************************
 * For initializing an array to integer zeros                             *
 **************************************************************************/
int *initialize(gottennodes,ngroup)
int *gottennodes,ngroup;
{
  static int i;
  for (i=0;i<ngroup;i++)
    gottennodes[i]=0;
  return(gottennodes);
}

/**************************************************************************
 * For creating a double vector of dimension Size and initialize to zero  *
 **************************************************************************/
VECTOR *v_get(Size)
int Size;
{
  VECTOR *v;
  double *d;
  static int i;
  
  d=CREATE(Size,double);
  v=CREATE(1,VECTOR);
  v->dim=Size;
  v->ve=d;
  return(v);
}  
/**************************************************************************
 * For initializing a double vector to ones                               *
 **************************************************************************/
VECTOR *vones(r)
VECTOR *r;
{
  static int i;
  for ( i=0;i<r->dim;i++ )
    r->ve[i]=1.0;
  return(r);
}  
   
/**************************************************************************
 * For permuting a double vector of dimension Size                        *
 **************************************************************************/
VECTOR *permdv(Size,s,permv)
int Size;
VECTOR *s;
int *permv;
{
  static int i, next, now;
  double stemp, temp;
  
  for ( i=1; i<=Size; i++ ) {
    if ( permv[i]>=0 ) {
      next=permv[i];
      stemp=s->ve[i-1];
L100:
      if ( permv[next] >= 0 ) {
        temp=stemp;
        stemp=s->ve[next-1];
        s->ve[next-1]=temp;
        now=next;
        next=permv[now];
        permv[now]=-next;
        goto L100;
      }
    }
  }
  for ( i=1; i<=Size; i++ ) 
    permv[i]=-permv[i];
  return(s);
}  
/**************************************************************************
 * For finding the 2-norm of a VECTOR                                     *
 **************************************************************************/
double norm2(r)
VECTOR *r;
{ static int i;
  double temp, result;
  
  temp=0.0;
  for ( i=0;i<r->dim;i++ ) {
    temp+=r->ve[i]*r->ve[i];
  }
  result=sqrt(temp);
  return(result);
}  

