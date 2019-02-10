/* last revised 25.03.1999                                             */
/* jacobi_seq.c = routine for solving the LSQ problem iteratively      */
/* and using jacobi updates                                            */
/* subspace corrected approach is used                                 */
/***********************************************************************/
/* #include "mpi.h" */
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include "data_struct.h"

#define CHOL  

#define PONTUS 

/* #define PERM    */
/* #define MMDP    */
/* #define COLP      */
/* #define NEWPERM    */

  VECTOR **D,**DS;
  MATRIX *C,*AA;

/* #define XERROR  */

#ifdef PERM
  int *p,*ivp;
  double *tmp;
#endif
 
#ifdef XERROR
  VECTOR *xnormv;
#endif

#define DIFFGROUPSIZE    
                                     
/**************************************************************************/
int main(argc,argv)
int argc;
char *argv[];
{  

/* local variables */
   
   static int i,j,k,l;
   int rnum, count, counter;
   char c1;
   double verdi;
   double number,num,*ip;
   int nc,converged,sind,eind,choice1;
   clock_t time1,time2;
   double dtime;
   int *nonzeroptr;
   
#ifdef DIFFGROUPSIZE   
   int choice2;
#endif

   sscanf(argv[1],"%d",&choice1);

#ifdef DIFFGROUPSIZE   
   sscanf(argv[2],"%d",&choice2); /* get the number of groups */
#endif       
            
   switch(choice1) {
      case 1:
        if ((fp=fopen("matrix1","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
#ifdef PONTUS
          if ((fp1=fopen("ash219.d","r"))==NULL) {  
#else
        if ((fp1=fopen("matrix1dat","r"))==NULL) {
#endif        
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case 2:
        if ((fp=fopen("matrix2","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
#ifdef PONTUS
          if ((fp1=fopen("ash958.d","r"))==NULL) { 
#else
          if ((fp1=fopen("matrix2dat","r"))==NULL) { 
#endif         
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case 3:
        if ((fp=fopen("matrix3","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
#ifdef PONTUS
          if ((fp1=fopen("ash331.d","r"))==NULL) {  
#else
          if ((fp1=fopen("matrix3dat","r"))==NULL) { 
#endif        
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case 4:
        if ((fp=fopen("matrix4","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
#ifdef PONTUS
          if ((fp1=fopen("ash608.d","r"))==NULL) {  
#else
          if ((fp1=fopen("matrix4dat","r"))==NULL) { 
#endif        
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case 5:
        if ((fp=fopen("matrix5","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
#ifdef PONTUS
          if ((fp1=fopen("abb313.d","r"))==NULL) { 
#else
          if ((fp1=fopen("matrix5dat","r"))==NULL) { 
#endif        
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case 6:
          if ((fp=fopen("wl1033dat","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          break;
      case 7:
        if ((fp=fopen("wl1252dat","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case 8:
        if ((fp=fopen("wl1364dat","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case 9:
        if ((fp=fopen("wl1641dat","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case 10:
        if ((fp=fopen("wl1850dat","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case 11:
        if ((fp=fopen("wl1991dat","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case 12:
        if ((fp=fopen("wl2808dat","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;   
   }
   if ( (choice1 >= 1) && (choice1<= 5) ){ 
      epsilon=1.0e-5;   
      readmatrix(fp,fp1); 
      xexact=v_get(A->n);   
      
#ifdef PONTUS 
       switch(choice1) {
         case 1:
           fp=fopen("ash219x","r");
           break;
         case 2:
           fp=fopen("ash958x","r");
           break;
         case 3:
           fp=fopen("ash331x","r");
           break;
         case 4:
           fp=fopen("ash608x","r");
           break;
         case 5:
           fp=fopen("abb313x","r");
           break;
       }
#else  
      switch(choice1) {
        case 1:
          fp=fopen("matrix1x","r");
          break;
        case 2:
          fp=fopen("matrix2x","r");
          break;
        case 3:
          fp=fopen("matrix3x","r");
          break;
        case 4:
          fp=fopen("matrix4x","r");
          break;
        case 5:
          fp=fopen("matrix5x","r");
          break;
      }
#endif      
      for (i=0;i<xexact->dim;i++) {
        fscanf(fp,"%le\n",&verdi);
        xexact->ve[i]=verdi;
      }
      fclose(fp);

#ifdef PERM

/* permute the matrix A */

      switch(choice1) {
        case 1:
#ifdef COLP
          fp=fopen("m1colperm","r");
#endif
#ifdef MMDP          
          fp=fopen("m1colmmd","r");
#endif          
#ifdef NEWPERM
          fp=fopen("newperm1","r");
#endif          
          break;
        case 2:
#ifdef COLP
          fp=fopen("m2colperm","r");
#endif
#ifdef MMDP          
          fp=fopen("m2colmmd","r");
#endif           
#ifdef NEWPERM
          fp=fopen("newperm2","r");
#endif          
          break;
        case 3:
#ifdef COLP
          fp=fopen("m3colperm","r");
#endif
#ifdef MMDP          
          fp=fopen("m3colmmd","r");
#endif           
#ifdef NEWPERM
          fp=fopen("newperm3","r");
#endif          
          break;
        case 4:
#ifdef COLP
          fp=fopen("m4colperm","r");
#endif
#ifdef MMDP          
          fp=fopen("m4colmmd","r");
#endif           
#ifdef NEWPERM
          fp=fopen("newperm4","r");
#endif   
          break;
        case 5:
#ifdef COLP
          fp=fopen("m5colperm","r");
#endif
#ifdef MMDP          
          fp=fopen("m5colmmd","r");
#endif          
#ifdef NEWPERM
          fp=fopen("newperm5","r");
#endif 
          break;
      }
      p=CREATE(A->n,int);
      --p;
      for (i=1;i<=A->n;i++) {
        fscanf(fp,"%d\n",&j);
        p[i]=j;
      }
      fclose(fp);
      ivp=CREATE(A->n,int);
      --ivp;
      for (i=1;i<=A->n;i++) {
        j=p[i];
        ivp[i]=j;
      }
      permuteA(A,p);
#endif
      
      rhs=v_get(A->m);

#ifdef PERM

/* permute xexact */

      tmp=CREATE(A->n,double);
      for (i=1;i<=A->n;i++) {
        j=p[i];
        tmp[i-1]=xexact->ve[j-1];
      }
      free((double *)(xexact->ve));
      xexact->ve=tmp;
#endif
      
      rhs=m_fulvectormlt(A,xexact,rhs);   
   }
   else if ( (choice1 >= 6) && (choice1<= 12) ) {
      epsilon=1.0e-3; 
      readbigmatrix(fp);
      xexact=v_get(A->n);  
      switch(choice1) {
         case 6:
           fp=fopen("wl1033x","r");
           break;
         case 7:
           fp=fopen("wl1252x","r");
           break;
         case 8:
           fp=fopen("wl1364x","r");
           break;
         case 9:
           fp=fopen("wl1641x","r");
           break;
         case 10:
           fp=fopen("wl1850x","r");
           break;
         case 11:
           fp=fopen("wl1991x","r");
           break;
         case 12:
           fp=fopen("wl2808x","r");
           break;
      } 
      for (i=0;i<xexact->dim;i++) {
        fscanf(fp,"%le\n",&verdi);
        xexact->ve[i]=verdi;
      }
      fclose(fp); 
      rhs=v_get(A->m);
      rhs=m_fulvectormlt(A,xexact,rhs);  
   }  

/* garbage collection */

#ifndef XERROR
    free((double *)(xexact->ve));
    free((VECTOR *)(xexact));      
#endif   
       
   residual=v_get(rhs->dim);
   x=v_get(A->n);              /* initial x vector, all zeros */
   for ( i=0;i<rhs->dim;i++ )
     residual->ve[i]=-1.0*rhs->ve[i];

/* garbage collection */

   free((double *)(rhs->ve));
   free((VECTOR *)(rhs)); 
   
/* form the non-overlapping column groups of A */

#ifdef DIFFGROUPSIZE
   ngroup=choice2;
#else   
   ngroup=floor((A->n)/10)+1;
#endif      
    
   sub_ptr=CREATE(ngroup,MATRIX *);
   zeroptr=CREATE(ngroup,int *);

#ifdef DIFFGROUPSIZE
   ip=CREATE(1,double);
   number=(double)(A->n)/(double)(ngroup);
   num=modf(number,ip);
   if (num < 0.5)
     nc=floor(number);
   else if (num >= 0.5) 
     nc=ceil(number); 
   free((double *)(ip));
   sub_ptr=create_nover_diff(A,sub_ptr,ngroup,zeroptr,nc);   
#else      
   sub_ptr=create_nover(A,sub_ptr,ngroup,zeroptr);
   nc=10;
#endif

/* garbage collection */

   free((int *)(A->Jind));
   free((int *)(A->Iind));
   free((int *)(A->colind));
   free((int *)(A->rowptr));
   free((double *)(A->Valuerow));
   free((double *)(A->Valuecol));
     
   r_zero_norm=norm2(residual);
   printf("\n initial norm: %17.15e",r_zero_norm);

/* allocate the necessary storage to keep the vectors for each group */

   rnzptr=CREATE(ngroup,double *);
   ddiagptr=CREATE(ngroup,double *);
   xrnzptr=CREATE(ngroup,int *);
   nzsubptr=CREATE(ngroup,int *);
   xnzsubptr=CREATE(ngroup,int *);
   permptr=CREATE(ngroup,int *);
   invpptr=CREATE(ngroup,int *);
   nonzeroptr=CREATE(ngroup,int);
   --rnzptr;
   --ddiagptr;
   --xrnzptr;
   --nzsubptr;
   --xnzsubptr;
   --permptr;
   --invpptr;
   --nonzeroptr;
   
/* initialize for SPARSPAK */  

   spbusr_.MAXSB=7000;
   TOL=0.0000000001;
   TYPTOL=1;
   spksys_.MAXINT=SHRT_MAX;
   spbusr_.IERRB=0;
   spbusr_.MSGLVB=1;
    
/* for each group do the factorization and store the necessary vectors */

   for ( l=0; l<ngroup; l++ ) {
     AA=sub_ptr[l];

     spbusr_.MDEQNS=0;
     spbusr_.MSCONS=0;
     spbusr_.MDCONS=0;
     spbcon_.NSEQNS=0;
     spbcon_.NDEQNS=0;
     spbcon_.NSCONS=0;
     spbcon_.NDCONS=0;
     spbcon_.NOFNZ=0;
     spbcon_.NOFSUB=0;
    
     spbmap_.LLFREE=spbusr_.MAXSB;

     spbcon_.NEDGES=0;
     spbcon_.NZMAX=0;
     spbcon_.NZEQNS=0;
     spbcon_.NZCONS=0;

     spbusr_.MCOLS=AA->n;
     spbusr_.MSEQNS=AA->m-zeroptr[l][0];
    
     spbcon_.NCOLS=spbusr_.MCOLS;  
        
     T=CREATE(spbusr_.MAXSB,int);
     for (i=0;i<spbusr_.MAXSB;i++) {
       T[i]=-(i+1);
     }
     --T;

/* read in the matrix into SPARSPAK */

     TYPE=1;   
     rnum=1;
    
     for ( j=1;j<=AA->m;j++ ) {
       if ( AA->rowptr[j] != 0 ) { 
         count=1;
         while ( AA->rowptr[j+count] == 0 )
           count++;
         NSUBS=AA->rowptr[j+count]-AA->rowptr[j];
         spbcon_.NOFNZ+=NSUBS;
         ROWNUM=rnum;
         rnum++;
         SUBS=CREATE(NSUBS,int);
         for ( k=0;k<NSUBS;k++ ) {
           *(SUBS+k)=AA->colind[AA->rowptr[j]+k];
         }
         --SUBS;
         inxywb(ROWNUM,TYPE,NSUBS,SUBS,T);
         SUBS++;
         free((int *)(SUBS)); 
       }
     }
  
     orcolb(T);
 
/* compute A^T*A and initialize lnz and diag */

     xrnzv=CREATE(spbcon_.NCOLS+1,int);
     xnzsubv=CREATE(spbcon_.NCOLS,int);
     nzsubv=CREATE(spbcon_.NOFNZ,int);
     permv=CREATE(spbcon_.NCOLS,int);
     invpv=CREATE(spbcon_.NCOLS,int);
     --permv;
     --invpv;
     --xrnzv;
     --xnzsubv;
     --nzsubv;
     j = spbcon_.NCOLS+1;
     for (i = 1; i <= j; i++) {
       xnzsubv[i] = T[xnzsub+i];
       xrnzv[i] = T[xrnz+i];       
     }
     
     for (i = 1; i <= spbcon_.NOFNZ; i++) {
	nzsubv[i] = T[nzsub+i];
     }
     j = spbcon_.NCOLS;
     for (i = 1; i <= j; i++) {
       permv[i] = T[perm+i];
       invpv[i] = T[invp+i];
     }
 
/* garbage collection */

     ++T;
     free((int *)(T));
     
     xrnzptr[l+1]=xrnzv;
     nzsubptr[l+1]=nzsubv;
     xnzsubptr[l+1]=xnzsubv;
     permptr[l+1]=permv;
     invpptr[l+1]=invpv;
     nonzeroptr[l+1]=spbcon_.NOFNZ;
     
/* get the necessary storage for A^TA */

     ddiag=CREATE(AA->n,double);
     lnzv=CREATE(nonzeroptr[l+1],double);
     --ddiag;
     --lnzv;
     
     j=nonzeroptr[l+1]+AA->n;
     isub=CREATE(j,int);
     jsub=CREATE(j,int);
     values=CREATE(j,double);
     --isub;
     --jsub;
     --values;
     
     sym_mult(AA,lnzv,ddiag,invpptr[l+1],xnzsubptr[l+1],nzsubptr[l+1],xrnzptr[l+1],isub,jsub,values);

/* garbage collection */

     isub++;
     jsub++;
     values++;
     free((int *)(isub));
     free((int *)(jsub));
     free((double *)(values));     

/* the next three vectors are initialized in gsfct */
   
     temp=CREATE(AA->n,double);     
     first=CREATE(AA->n,int);
     llink=CREATE(AA->n,int);

     --llink;
     --temp;
     --first;
   
/* do Cholesky factorization */ 

     gsfct(AA->n,xrnzptr[l+1],lnzv,xnzsubptr[l+1],nzsubptr[l+1],ddiag,llink,first,temp,flag__); 

/* garbage collection */

     ++llink;
     ++temp;
     ++first;
     free((int *)(llink));
     free((double *)(temp));
     free((int *)(first));
   
     ddiagptr[l+1]=ddiag;
     rnzptr[l+1]=lnzv;
     
   } 
   
/* end of for loop for all groups */  

/* do the symmetric factorization of the C matrix */

   C=CREATE(1,MATRIX);
   C->m=residual->dim;     
   C->n=ngroup;
   createC1(sub_ptr,C,ngroup,nc,zeroptr);

   spbusr_.MDEQNS=0;
   spbusr_.MSCONS=0;
   spbusr_.MDCONS=0;
   spbcon_.NSEQNS=0;
   spbcon_.NDEQNS=0;
   spbcon_.NSCONS=0;
   spbcon_.NDCONS=0;
   spbcon_.NOFNZ=0;
   spbcon_.NOFSUB=0;
    
   spbmap_.LLFREE=spbusr_.MAXSB;

   spbcon_.NEDGES=0;
   spbcon_.NZMAX=0;
   spbcon_.NZEQNS=0;
   spbcon_.NZCONS=0;
   spbusr_.MCOLS=C->n;
   spbusr_.MSEQNS=C->m;
    
   spbcon_.NCOLS=spbusr_.MCOLS;  
   spbusr_.MAXSB=6000;
       
   T=CREATE(spbusr_.MAXSB,int);
   for (i=0;i<spbusr_.MAXSB;i++) {
     T[i]=-(i+1);
   }
   --T;

/* read in the matrix into SPARSPAK */

   TYPE=1;   
   rnum=1;
   for ( j=1;j<=C->m;j++ ) {
     if ( C->rowptr[j] != 0 ) { 
       count=1;
       while ( C->rowptr[j+count] == 0 )
         count++;
       NSUBS=C->rowptr[j+count]-C->rowptr[j];
       spbcon_.NOFNZ+=NSUBS;
       ROWNUM=rnum;
       rnum++;
       SUBS=CREATE(NSUBS,int);
       for ( k=0;k<NSUBS;k++ ) {
         *(SUBS+k)=C->colind[C->rowptr[j]+k];
       }
       --SUBS;
       inxywb(ROWNUM,TYPE,NSUBS,SUBS, T);
       SUBS++;
       free((int *)SUBS);  
     }
   }     
   orcolb(T);
   xrnzv=CREATE(spbcon_.NCOLS+1,int);
   xnzsubv=CREATE(spbcon_.NCOLS,int);
   nzsubv=CREATE(spbcon_.NOFNZ,int);
   permv=CREATE(spbcon_.NCOLS,int);
   invpv=CREATE(spbcon_.NCOLS,int);
   --permv;
   --invpv;
   --xrnzv;
   --xnzsubv;
   --nzsubv;
   j = spbcon_.NCOLS+1;
   for (i = 1; i <= j; i++) {
     xrnzv[i] = T[xrnz+i];
     xnzsubv[i] = T[xnzsub+i];
   }
   
   for (i = 1; i <= spbcon_.NOFNZ; i++) {
     nzsubv[i] = T[nzsub+i];
   }
   j = spbcon_.NCOLS;
   for (i = 1; i <= j; i++) {
     permv[i] = T[perm+i];
     invpv[i] = T[invp+i];
   }
   
/* garbage collection */

   ++T;
   free((int *)(T));        

/* begin the solution */

#ifdef XERROR
   xnormv=v_get(xexact->dim);
   for (i=0;i<xnormv->dim;i++)
     xnormv->ve[i]=x->ve[i]-xexact->ve[i];
   err=norm2(xnormv);
#else
   r_norm=r_zero_norm;
   err=r_norm/r_zero_norm;   
#endif   

   converged=0;
   omega=1.0;
   counter=0;

   D=CREATE(ngroup,VECTOR *);
   --D;

   DS=CREATE(ngroup,VECTOR *);
   --DS;

printf("\n before loop..");       
#pragma parallel shared(converged,sub_ptr,DS,D,residual,C,k) 
{ 

#pragma one processor
{   
   time1=clock();
}  
   while (converged==0) {

#pragma pfor local(i,AA,sol,s) affinity (i) = thread (i)
{  
     for (i=1;i<=ngroup;i++) {
       AA=sub_ptr[i-1];
       sol=v_get(AA->m);
       s=v_get(AA->n);
      
       s=multAtransr3(AA,residual,s);
       s=permdv(AA->n,s,invpptr[i]);
       s=gsslv(AA->n,xrnzptr[i],rnzptr[i],xnzsubptr[i],nzsubptr[i],ddiagptr[i],s);
       
/* permute the solution to its original order of columns */

       s=permdv(AA->n,s,permptr[i]);      
       sol=m_fulvectormlt(AA,s,sol);
             
       DS[i]=sol;
       D[i]=s;    
     }
}     

#pragma one processor
{
                    
/* compute C^T*C and initialize lnz and diag */

     ddiag=CREATE(spbcon_.NCOLS,double);
     lnzv=CREATE(spbcon_.NOFNZ,double);
     --ddiag;
     --lnzv;
     
/* get the necessary storage for C^TC */

     i=spbcon_.NOFNZ+spbcon_.NCOLS;
     isub=CREATE(i,int);
     jsub=CREATE(i,int);
     values=CREATE(i,double);
     --isub;
     --jsub;
     --values;
     createC22(C,DS,ngroup);
     sym_mult3(C,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values);  
  
/* garbage collection */

     isub++;
     jsub++;
     values++;
     free((int *)(isub));
     free((int *)(jsub));
     free((double *)(values));
   
/* the next three vectors are initialized in gsfct */
   
     llink=CREATE(spbcon_.NCOLS,int);
     temp=CREATE(spbcon_.NCOLS,double);
     first=CREATE(spbcon_.NCOLS,int);
     --llink;
     --temp;
     --first;

     gsfct(spbcon_.NCOLS,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,llink,first,temp,flag__); 

/* garbage collection */

     ++llink;
     ++temp;
     ++first;
     free((int *)(llink));
     free((double *)(temp));
     free((int *)(first));
     
/* do the solution */

     s=v_get(C->n);
     s=multAtransr3(C,residual,s);
     s=permdv(C->n,s,invpv);
     s=gsslv(C->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,s); 
     s=permdv(C->n,s,permv);  

/* update x vector */
  
     k=0;
     for (i=1;i<=ngroup;i++) {
       for (j=1;j<=D[i]->dim;j++) {
         x->ve[k]+=D[i]->ve[j-1]*s->ve[i-1];
         k++;
       }
     }
      
/* update the residual vector */

       for (i=1;i<=ngroup;i++) {
       sind=C->Jind[i];
       eind=C->Jind[i+1]-1;
       for (j=sind;j<=eind;j++) {
         residual->ve[C->Iind[j]-1]+=DS[i]->ve[C->Iind[j]-1]*s->ve[i-1];
       }  
     } 
/* find the error */

     r_norm=norm2(residual);

/* printf("\n residual norm: %16.14e\n",r_norm);          */
       
     err=r_norm/r_zero_norm;
                 
     counter++;

     if ((err < epsilon)||(r_norm==0.0)||(counter>15000)) {
       converged=1;
     }
     else {
     
/* garbage collection */

       for (i=1;i<=ngroup;i++)
         free((double *)(DS[i]->ve));
       free((double *)(C->Valuerow));
       free((double *)(C->Valuecol));
       free((double *)(s->ve));
       free((VECTOR *)(s));
       ++ddiag;
       ++lnzv;
       free((double *)(ddiag));
       free((double *)(lnzv));    
       for (i=1;i<=ngroup;i++)
         free((double *)(D[i]->ve));         
     } 
}  /* end of one processor */
   }  /* end of while */
} /* end of parallel */

   time2=clock();
   dtime=(double) (time2)- (double) (time1);
   dtime=dtime/CLOCKS_PER_SEC;
   printf("\n Converged after %d loops.",counter);
#ifdef XERROR
     printf("\n The error is: %16.14e \n",err   );
#else   
     printf("\n Two norm of residual is: %16.14e \n",r_norm); 
#endif    
   printf("\n Elapsed time: %lf seconds \n",dtime);    
   exit(0); 
}    
       
