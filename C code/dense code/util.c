#include "datastruct_dense.h"

/**************************************************************************
 * For creating a double vector of dimension Size initialized to zero     *
 **************************************************************************/
VECTOR *v_get(Size)
int Size;
{
  VECTOR *v;
  double *d;
  
  d=CREATE(Size,double);
  v=CREATE(1,VECTOR);
  v->dim=Size;
  v->ve=d;
  return(v);
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
