/* last revised 29.06.1998 */
/* solution.c called by the main program */

#include "data_struct.h" 

/**************************************************************************
 * For backward and forward solution of a system.                         *
 **************************************************************************/
VECTOR *solution(nncols,s,rnzv,ddiag,xnzsubv,nzsubv,xrnzv,permv,invpv)
int nncols;
VECTOR *s;
double *rnzv,*ddiag;
int *xnzsubv,*nzsubv,*xrnzv,*permv,*invpv;
{ int i,j,ii,jj,istrt,istop,icol,isub,now,next;
  double rhsj,ss,save,temp;

#ifdef DEBUG    
  fprintf(fp1,"\n inside solution...");
#endif  
  
  for ( i=1;i<=nncols;i++ ) {
    if ( invpv[i]>=0 ) {
      next=invpv[i];
      save=s->ve[i-1];
      while ( invpv[next]>=0 ) {
        temp=save;
        save=s->ve[next-1];
        s->ve[next-1]=temp;
        now=next;
        next=invpv[now];
        invpv[now]=-1*next;

      }
    }
  }
  for ( i=1;i<=nncols;i++ ) {
    invpv[i]=-1*invpv[i];
  
  }
  for ( j=1;j<=nncols;j++ ) {
    rhsj=s->ve[j-1]/ddiag[j];

    s->ve[j-1]=rhsj;
    istrt=xrnzv[j];
    istop=xrnzv[j+1]-1;
    if ( istop>=istrt ) {
      i=xnzsubv[j];
      for ( ii=istrt;ii<=istop;ii++ ) {
        icol=nzsubv[i];
        s->ve[icol-1]=s->ve[icol-1]-rnzv[ii]*rhsj;
        i++;
      }
    }
  }
  j=nncols;
  for ( jj=1;jj<=nncols;jj++ ) {
    ss=s->ve[j-1];
    istrt=xrnzv[j];
    istop=xrnzv[j+1]-1;
    if ( istop>=istrt ) {
      i=xnzsubv[j];
      for ( ii=istrt;ii<=istop;ii++ ) {
        isub=nzsubv[i];
        ss=ss-rnzv[ii]*s->ve[isub-1];
        i++;
      }
    }
    s->ve[j-1]=ss/ddiag[j];
    j--;
  }
  for ( i=1;i<=nncols;i++ ) {
    if ( permv[i]>=0 ) {
      next=permv[i];
      save=s->ve[i-1];
      while ( permv[next]>=0 ) {
        temp=save;
        save=s->ve[next-1];
        s->ve[next-1]=temp;
        now=next;
        next=permv[now];
        permv[now]=-1*next;

      }
    }
  }
  for ( i=1;i<=nncols;i++ ) {
    permv[i]=-1*permv[i];
    
  }
  return(s);
}
