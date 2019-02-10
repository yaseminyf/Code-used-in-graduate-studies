/* last revised 24.11.1998                                             */
/* iter_seq_jacobi.c = routine for solving the LSQ problem iteratively */
/* and using jacobi updates                                            */
/* it is possible to choose either the QR or the Cholesky approach     */
/* by changing the 'define' switch                                     */
/* this version computes Jacobi with p=DS                              */
/* the equations are reformulated to save computation                  */
/***********************************************************************/
/* #include "mpi.h" */
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include "data_struct.h"

#define CHOL  
  
VECTOR **D,**Ad,**AD;
MATRIX *hatA,*tildeA,**tildeAptr;
VECTOR **Apptr,**tmpApptr;
VECTOR *pvec,**pvecptr,*pd;
MATRIX *AA;
double *deltav,*sndvector;
VECTOR *vvec,*zvec,**zvecptr,*littles;
int antall; 

#define PONTUS

#define XERROR   /* to check convergence with error on x */

#define FMSTART  /* starting values of p are F_M values */

#ifdef FMSTART
  VECTOR *evec;
#endif 

#ifdef XERROR
  VECTOR *xnormv;
#endif  
 
#define DIFFGROUPSIZE       
                                 
/**************************************************************************/
int main()

{  

/* local variables */
   
   static int i,j,k,l,t,tt;
   int rnum, count, counter;
   char choice1,c1;
   double verdi;
   double number,num,*ip;
   int nc;
   clock_t time1,time2;
   double dtime;
   int *nonzeroptr;
#ifdef DIFFGROUPSIZE   
   int choice2;
#endif

   printf("\n Enter the number of the matrix you want to use..");
   printf("\n              [1,5]\n");
   choice1=getchar(); 

#ifdef DIFFGROUPSIZE   
   printf("\n Enter the number of groups..");
   scanf("%d",&choice2); 
#endif   
   antall=1;
    
/*   choice1='1'; */
         
   switch(choice1) {
     case '1':
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
     case '2':
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
     case '3':
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
     case '4':
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
     case '5':
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
   }
if ( (choice1 >= '1') && (choice1<= '5') ){
#ifdef XERROR
       epsilon=1.0e-6;
#else
       epsilon=1.0e-5;      
#endif            
       readmatrix(fp,fp1); 
       xexact=v_get(A->n);  

#ifdef PONTUS 
     switch(choice1) {
       case '1':
         fp=fopen("ash219x","r");
         break;
       case '2':
         fp=fopen("ash958x","r");
         break;
       case '3':
         fp=fopen("ash331x","r");
         break;
       case '4':
         fp=fopen("ash608x","r");
         break;
       case '5':
         fp=fopen("abb313x","r");
         break;
     }
#else              
     switch(choice1) {
       case '1':
         fp=fopen("matrix1x","r");
         break;
       case '2':
         fp=fopen("matrix2x","r");
         break;
       case '3':
         fp=fopen("matrix3x","r");
         break;
       case '4':
         fp=fopen("matrix4x","r");
         break;
       case '5':
         fp=fopen("matrix5x","r");
         break;
     }
#endif       
     for (i=0;i<xexact->dim;i++) {
       fscanf(fp,"%le\n",&verdi);
       xexact->ve[i]=verdi;
     }
     fclose(fp); 

     fp1=fopen("data","w");
   
      rhs=v_get(A->m);     
      rhs=m_fulvectormlt(A,xexact,rhs); 

/* garbage collection */

#ifndef XERROR
      free((double *)(xexact->ve));
      free((VECTOR *)(xexact)); 
#endif            
   }
   else if ( (choice1 >= '6') && (choice1<= '9') ) {
      epsilon=1.0e-01;
      fp1=fopen("data","w");
      readnewmatrix(fp,fp1);
   }

   x=v_get(A->n);
   residual=v_get(rhs->dim);   
   for ( i=0;i<rhs->dim;i++ )
     residual->ve[i]=-1.0*rhs->ve[i];

/* garbage collection */

   free((double *)(rhs->ve));
   free((VECTOR *)(rhs));
     
/* record the initial data in a file */

   fprintf(fp1,"\n\nThe initial residual:\n\n");
   for ( i=0;i<residual->dim;i++ ) {
     fprintf(fp1,"%16.15e ",residual->ve[i]);
     if ( ((i+1)%5)==0 )
       fprintf(fp1,"\n");
   } 

/* form the non-overlapping column groups of A */

   if ( (choice1 >= '1') && (choice1<= '5') ) {
#ifdef DIFFGROUPSIZE
     ngroup=choice2;
#else   
     ngroup=floor((A->n)/10)+1;
#endif      
   }
   else if ( (choice1 >= '6') && (choice1<= '9') ) {
#ifdef DIFFGROUPSIZE
     ngroup=choice2;
#else    
     ngroup=floor((A->n)/25)+1;
#endif
   }      
   sub_ptr=CREATE(ngroup,MATRIX *);
   zeroptr=CREATE(ngroup,int *);
   if ( (choice1 >= '1') && (choice1<= '5') ) {

#ifdef DIFFGROUPSIZE
     sub_ptr=create_nover_diff(A,sub_ptr,ngroup,zeroptr); 
     ip=CREATE(1,double);
     number=(double)(A->n)/(double)(ngroup);
     num=modf(number,ip);
     if (num < 0.5)
       nc=floor(number);
     else if (num >= 0.5) 
       nc=ceil(number); 
     free((double *)(ip));  
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
/*     free((MATRIX *)(A));  */
     
     r_zero_norm=norm2(residual);
     printf("\n initial norm: %16.14e",r_zero_norm);
   }
   else if ( (choice1 >= '6') && (choice1<= '9') ){
     sub_ptr=create_nover25(A,sub_ptr,ngroup,zeroptr);
     nc=25;
     rvec2=v_get(A->n);
     rvec2=multAtransr2(A,residual,rvec2);
     r_zero_norm=norm2(rvec2);
   }

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

   spbusr_.MAXSB=6000;
   TOL=0.0000000001;
   TYPTOL=1;

   spksys_.MAXINT=SHRT_MAX;
    
   spbusr_.IERRB=0;
   spbusr_.MSGLVB=1;     

/* for each group do the factorization and store the necessary vectors */

   time1=clock();
   pvecptr=CREATE(ngroup,VECTOR *);
   for ( l=0; l<ngroup; l++ ) {
     AA=sub_ptr[l];
     pvec=v_get(sub_ptr[l]->n);        
#ifdef FMSTART
     evec=v_get(sub_ptr[l]->n);
     evec=vones(evec);
     sol=v_get(sub_ptr[l]->m);
     sol=m_fulvectormlt(sub_ptr[l],evec,sol);
     s=v_get(sub_ptr[l]->n);
     s=multAtransr2(sub_ptr[l],sol,s);
#endif   
     for (j=0;j<pvec->dim;j++) {         
#ifdef FMSTART
         pvec->ve[j]=1.0/s->ve[j];
#else
       pvec->ve[j]=1.0;
#endif       
     }         
#ifdef FMSTART     
     free((double *)(evec->ve));
     free((VECTOR *)(evec));
     free((double *)(sol->ve));
     free((VECTOR *)(sol));
     free((double *)(s->ve));
     free((VECTOR *)(s)); 
#endif 
     pvecptr[l]=pvec;
   }
   Apptr=CREATE(ngroup,VECTOR *);
   for ( l=0; l<ngroup; l++ ) {
     AA=sub_ptr[l];
     pvec=pvecptr[l];
     sol=v_get(AA->m);
     sol=m_fulvectormlt(AA,pvec,sol); 
     Apptr[l]=sol;
   }

/* form tildeA pointer */
   
   tildeAptr=CREATE(ngroup,MATRIX *);
   tildeAptr=createtildeAptr(tildeAptr,sub_ptr,Apptr,ngroup);
             
   for ( l=0; l<ngroup; l++ ) {
     AA=tildeAptr[l];
     
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
     spbusr_.MSEQNS=AA->m;
    
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
     
     xrnzptr[l+1]=xrnzv;
     nzsubptr[l+1]=nzsubv;
     xnzsubptr[l+1]=xnzsubv;
     permptr[l+1]=permv;
     invpptr[l+1]=invpv;
     nonzeroptr[l+1]=spbcon_.NOFNZ;     
   } 
   
/* end of for loop for all groups */

   time2=clock();
   dtime=(double) (time2) - (double) (time1);
   dtime=dtime/CLOCKS_PER_SEC;
   printf("\n Factorization time: %lf seconds \n",dtime);    

/* do the symmetric factorization of hatA */

   hatA=CREATE(1,MATRIX);
   hatA->m=residual->dim;     
   hatA->n=ngroup;
   createC1(sub_ptr,hatA,ngroup);

/* factorize hatA */

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
   spbusr_.MCOLS=hatA->n;
   spbusr_.MSEQNS=hatA->m;
    
   spbcon_.NCOLS=spbusr_.MCOLS;  
   spbusr_.MAXSB=6000;
   T=CREATE(spbusr_.MAXSB,int);
printf("\n created T...");   
   for (i=0;i<spbusr_.MAXSB;i++) {
     T[i]=-(i+1);
   }
   --T;
   
/* read in the matrix into SPARSPAK */

   TYPE=1;   
   rnum=1;
   for ( j=1;j<=hatA->m;j++ ) {
     if ( hatA->rowptr[j] != 0 ) { 
       count=1;
       while ( hatA->rowptr[j+count] == 0 )
         count++;
       NSUBS=hatA->rowptr[j+count]-hatA->rowptr[j];
       spbcon_.NOFNZ+=NSUBS;
       ROWNUM=rnum;
       rnum++;
       SUBS=CREATE(NSUBS,int);
       for ( k=0;k<NSUBS;k++ ) {
         *(SUBS+k)=hatA->colind[hatA->rowptr[j]+k];
       }
       --SUBS;
       inxywb(ROWNUM,TYPE,NSUBS,SUBS,T);
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
   counter=0;

   i=residual->dim+nc+ngroup;
   sndvector=CREATE(i,double);
     
   D=CREATE(ngroup,VECTOR *);
   --D;
   
   Ad=CREATE(ngroup,VECTOR *);
   --Ad;

   AD=CREATE(ngroup,VECTOR *);
   --AD;
   
   tmpApptr=CREATE(ngroup,VECTOR *);
   
printf("\n before loop...");   
   time1=clock();
   
   while (err >= epsilon) {
     if ( counter > 0 ) {
       for (t=1;t<=antall;t++) {
         for (i=1;i<=ngroup;i++) {
           tildeA=tildeAptr[i-1];
           A=sub_ptr[i-1];
           s=v_get(tildeA->n);
           s=multAtransr3(tildeA,vvec,s);
           s=permdv(tildeA->n,s,invpptr[i]);
           s=gsslv(tildeA->n,xrnzptr[i],rnzptr[i],xnzsubptr[i],nzsubptr[i],ddiagptr[i],s);
           s=permdv(tildeA->n,s,permptr[i]);   
           sndvector[0]=tildeA->n;
       
           for (j=0;j<s->dim;j++)
             sndvector[tildeA->m+1+j]=s->ve[j];
           littles=v_get(A->n);
           for (j=0;j<A->n;j++) {
             littles->ve[j]=s->ve[i-1+j];
           }            
           free((double *)(s->ve));
           free((VECTOR *)(s));
           sol=v_get(A->m);
           sol=m_fulvectormlt(A,littles,sol);    
           for (j=0;j<sol->dim;j++)
             sndvector[j+1]=sol->ve[j];
           free((double *)(sol->ve));
           free((double *)(littles->ve));
           free((VECTOR *)(sol));
           free((VECTOR *)(littles));
           sol=v_get(residual->dim);
           for (j=0;j<sol->dim;j++)
             sol->ve[j]=sndvector[j+1];
           s=v_get((int)(sndvector[0]));
           for (j=0;j<s->dim;j++)
             s->ve[j]=sndvector[sol->dim+1+j];
           D[i]=s;
           Ad[i]=sol;
         }

         ddiag=CREATE(spbcon_.NCOLS,double);
         lnzv=CREATE(spbcon_.NOFNZ,double);
         --ddiag;
         --lnzv;
         for (j=1;j<=spbcon_.NOFNZ;j++)
           lnzv[j]=0.0;
         for (j=1;j<=spbcon_.NCOLS;j++)
           ddiag[j]=0.0;    
         i=spbcon_.NOFNZ+spbcon_.NCOLS;
         isub=CREATE(i,int);
         jsub=CREATE(i,int);
         values=CREATE(i,double);
         --isub;
         --jsub;
         --values;
         pd=v_get(x->dim);
         for (i=1;i<=ngroup;i++) {
           count=0;
           k=0;
           for (j=1;j<=ngroup;j++) {
             if (i!=j) {
               for (l=1;l<=pvecptr[j-1]->dim;l++) {
                 pd->ve[k]+=pvecptr[j-1]->ve[l-1]*D[i]->ve[count];
                 k++;
               }
               count++;
             }
             else {
               tt=D[i]->dim-ngroup+1;
               for (l=1;l<=tt;l++) {
                 pd->ve[k]+=D[i]->ve[count];
                 count++;
                 k++;
               }
             }
           }
         }
         deltav=CREATE(ngroup,double);
         --deltav;
         for (i=1;i<=ngroup;i++) {
           deltav[i]=0.0;
           for (j=1;j<=ngroup;j++) {
             if (i!=j) {
               if (j>i) {
                 deltav[i]+=D[j]->ve[i-1];
               }
               if (j<i) {
                 deltav[i]+=D[j]->ve[j-1+nc+i-j-1];
               }
             }
           }
         }   
         for (i=1;i<=ngroup;i++) {
           for (j=0;j<Ad[i]->dim;j++) {
             Ad[i]->ve[j]+=Apptr[i-1]->ve[j]*deltav[i];
           }
           AD[i]=Ad[i];
         }
         createC22(hatA,AD);  
         sym_mult(hatA,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values);
         isub++;
         jsub++;
         values++;
         free((int *)(isub));
         free((int *)(jsub));
         free((double *)(values));
         llink=CREATE(spbcon_.NCOLS,int);
         temp=CREATE(spbcon_.NCOLS,double);
         first=CREATE(spbcon_.NCOLS,int);
         --llink;
         --temp;
         --first;
         gsfct(spbcon_.NCOLS,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,llink,first,temp,flag__); 
         ++llink;
         ++temp;
         ++first;
         free((int *)(llink));
         free((double *)(temp));
         free((int *)(first));
         s=v_get(hatA->n);
         s=multAtransr3(hatA,residual,s);
         s=permdv(hatA->n,s,invpv);
         s=gsslv(hatA->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,s);   
         s=permdv(hatA->n,s,permv);  
         k=0;
         for (i=1;i<=ngroup;i++) {
           tt=D[i]->dim-ngroup+1;
           for (j=1;j<=tt;j++) {
             verdi=pd->ve[k]*s->ve[i-1];
             zvecptr[i-1]->ve[j-1]+=verdi;
             k++;
           }
         }
         for (i=1;i<=ngroup;i++) {
           for ( j=0;j<residual->dim;j++ ) {
             verdi=s->ve[i-1]*AD[i]->ve[j];
             vvec->ve[j]+=verdi;
             tmpApptr[i-1]->ve[j]+=verdi; 
           }
         }
         free((double *)(s->ve));
         free((VECTOR *)(s));
         ++ddiag;
         ++lnzv;
         free((double *)(ddiag));
         free((double *)(lnzv));
         free((double *)(pd->ve));
         free((VECTOR *)(pd));
         ++deltav;
         free((double *)(deltav));
         for (i=1;i<=ngroup;i++) {
           free((double *)(AD[i]->ve));
           free((double *)(Ad[i]->ve));
           free((double *)(D[i]->ve));
           free((VECTOR *)(AD[i]));
           free((VECTOR *)(Ad[i]));
           free((VECTOR *)(D[i])); 
         }
         free((double *)(hatA->Valuerow));
         free((double *)(hatA->Valuecol));
       }
       for (i=0;i<ngroup;i++) {
         for (j=0;j<pvecptr[i]->dim;j++)
           pvecptr[i]->ve[j]=zvecptr[i]->ve[j]; 
         free((double *)(zvecptr[i]->ve));
         free((VECTOR *)(zvecptr[i]));
       } 
       for (i=0;i<ngroup;i++) {
         for (j=0;j<residual->dim;j++)
           Apptr[i]->ve[j]=tmpApptr[i]->ve[j];
         free((double *)(tmpApptr[i]->ve));
         free((VECTOR *)(tmpApptr[i]));
       }
       
       free((VECTOR **)(zvecptr));
       free((double *)(vvec->ve));
       free((VECTOR *)(vvec)); 
       for (i=1;i<=ngroup;i++) {
         tildeA=tildeAptr[i-1];
         A=sub_ptr[i-1];
         putvaluesintildeA3(tildeA,Apptr,A,ngroup,i); 
         tildeAptr[i-1]=tildeA;
       }        
 
     } 
     
     for (i=1;i<=ngroup;i++) {
       AA=tildeAptr[i-1];

/* factorize tildeA */

       ddiag=CREATE(AA->n,double);
       lnzv=CREATE(nonzeroptr[i],double);
       --ddiag;
       --lnzv;
       for (j=1;j<=nonzeroptr[i];j++)
         lnzv[j]=0.0;
       for (j=1;j<=AA->n;j++)
         ddiag[j]=0.0;  
       j=nonzeroptr[i]+AA->n;
       isub=CREATE(j,int);
       jsub=CREATE(j,int);
       values=CREATE(j,double);
       --isub;
       --jsub;
       --values;  
       sym_mult(AA,lnzv,ddiag,invpptr[i],xnzsubptr[i],nzsubptr[i],xrnzptr[i],isub,jsub,values);

/* garbage collection */

       isub++;
       jsub++;
       values++;
       free((int *)(isub));
       free((int *)(jsub));
       free((double *)(values));
   
/* the next three vectors are initialized in gsfct */
   
       llink=CREATE(AA->n,int);
       temp=CREATE(AA->n,double);
       first=CREATE(AA->n,int);
       
       --llink;
       --temp;
       --first;

       gsfct(AA->n,xrnzptr[i],lnzv,xnzsubptr[i],nzsubptr[i],ddiag,llink,first,temp,flag__); 

/* garbage collection */
      
       ++llink;
       ++temp;
       ++first;
       free((int *)(llink));
       free((double *)(temp));
       free((int *)(first));
                       
       s=v_get(AA->n);
       s=multAtransr3(AA,residual,s);
       s=permdv(AA->n,s,invpptr[i]);
       s=gsslv(AA->n,xrnzptr[i],lnzv,xnzsubptr[i],nzsubptr[i],ddiag,s);
         
/* permute the solution to its original order of colums */

       s=permdv(AA->n,s,permptr[i]);   

/* copy this s to the vector to be sent to the master */ 
     
       for (j=0;j<s->dim;j++)
         sndvector[AA->m+1+j]=s->ve[j];    

/* take the littles out of s */

       littles=v_get(sub_ptr[i-1]->n);
       for (j=0;j<sub_ptr[i-1]->n;j++) {
         littles->ve[j]=s->ve[i-1+j];
       }                

       free((double *)(s->ve));
       free((VECTOR *)(s));
       
/* compute A*d */

       sol=v_get(sub_ptr[i-1]->m);
       sol=m_fulvectormlt(sub_ptr[i-1],littles,sol);    
       
/* copy this sol to the vector to be sent to the master */
       
       sndvector[0]=tildeAptr[i-1]->n;
       for (j=0;j<sol->dim;j++)
         sndvector[j+1]=sol->ve[j];       

/* garbage collection */

       free((double *)(sol->ve));
       free((double *)(littles->ve));
       free((VECTOR *)(sol));
       free((VECTOR *)(littles));
       rnzptr[i]=lnzv;
       ddiagptr[i]=ddiag;
       sol=v_get(residual->dim);
       for (j=0;j<sol->dim;j++)
         sol->ve[j]=sndvector[j+1];
       s=v_get((int)(sndvector[0]));
       for (j=0;j<s->dim;j++)
         s->ve[j]=sndvector[sol->dim+1+j];
       D[i]=s;
       Ad[i]=sol;
     }       
       
          
/* compute hatA^T*hatA and initialize lnz and diag */

     ddiag=CREATE(spbcon_.NCOLS,double);
     lnzv=CREATE(spbcon_.NOFNZ,double);
     --ddiag;
     --lnzv;
     for (j=1;j<=spbcon_.NOFNZ;j++)
       lnzv[j]=0.0;
     for (j=1;j<=spbcon_.NCOLS;j++)
       ddiag[j]=0.0;    
     
/* get the necessary storage for hatA^ThatA */

     i=spbcon_.NOFNZ+spbcon_.NCOLS;
     isub=CREATE(i,int);
     jsub=CREATE(i,int);
     values=CREATE(i,double);
     --isub;
     --jsub;
     --values;

/* compute D matrix first */

     pd=v_get(x->dim);
     for (i=1;i<=ngroup;i++) {
       count=0;
       k=0;
       for (j=1;j<=ngroup;j++) {
         if (i!=j) {
           for (l=1;l<=pvecptr[j-1]->dim;l++) {
             pd->ve[k]+=pvecptr[j-1]->ve[l-1]*D[i]->ve[count];
             k++;
           }
           count++;
         }
         else {
           for (l=1;l<=sub_ptr[j-1]->n;l++) {
             pd->ve[k]+=D[i]->ve[count];
             count++;
             k++;
           }
         }
       }
     }

/* compute the delta values for each group */

     deltav=CREATE(ngroup,double);
     --deltav;
     for (i=1;i<=ngroup;i++) {
       deltav[i]=0.0;
       for (j=1;j<=ngroup;j++) {
         if (i!=j) {
           if (j>i) {
             deltav[i]+=D[j]->ve[i-1];
           }
           if (j<i) {
             deltav[i]+=D[j]->ve[j-1+nc+i-j-1];
           }
         }
       }
     }
     
/* compute the columns of A*D=hatA matrix */
   
     for (i=1;i<=ngroup;i++) {
       for (j=0;j<Ad[i]->dim;j++) {
         Ad[i]->ve[j]+=Apptr[i-1]->ve[j]*deltav[i];
       }
       AD[i]=Ad[i];
     }
          
/* put in the values of hatA */

     createC22(hatA,AD);  
     sym_mult(hatA,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values);

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

     s=v_get(hatA->n);
     s=multAtransr3(hatA,residual,s);
     s=permdv(hatA->n,s,invpv);
     s=gsslv(hatA->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,s);   
     s=permdv(hatA->n,s,permv);  

/* update x vector and zvecptr at the same time */
     
     zvecptr=CREATE(ngroup,VECTOR *);   
          
     k=0;
     for (i=1;i<=ngroup;i++) {         
       t=D[i]->dim-ngroup+1;
       zvec=v_get(t);              
       for (j=1;j<=t;j++) {
         verdi=pd->ve[k]*s->ve[i-1];
         x->ve[k]+=verdi;
         zvec->ve[j-1]=verdi;
         k++;       
       }
       zvecptr[i-1]=zvec;              
     }
     
/* update the residual vector and vvec at the same time */

printf("\n");    
for (i=0;i<s->dim;i++)
printf("\n %16.14e",s->ve[i]); 
     vvec=v_get(residual->dim);   
     for (i=1;i<=ngroup;i++) {
       sol=v_get(residual->dim);
       for ( j=0;j<residual->dim;j++ ) {
         verdi=s->ve[i-1]*AD[i]->ve[j];
         residual->ve[j]+=verdi;
         vvec->ve[j]=residual->ve[j]+verdi;
         sol->ve[j]=verdi;  
       }
       tmpApptr[i-1]=sol;
     }
     
/* find the error */

#ifdef XERROR
     for (i=0;i<xnormv->dim;i++)
       xnormv->ve[i]=x->ve[i]-xexact->ve[i];
     err=norm2(xnormv);
     r_norm=norm2(residual);

printf("\n error: %16.14e\n",err);  
printf("\n residual norm: %16.14e\n",r_norm);        

#else
     if ( choice1 <= '5' ) {
       r_norm=norm2(residual);

printf("\n residual norm: %16.14e\n",r_norm);       
       
     }
     else {
       rvec2=multAtransr2(A,sol,rvec2);

/* if choice is bigger than 5 form a sol vector while updating residual */         
       free((double *)(sol->ve));
       free((VECTOR *)(sol));
       r_norm=norm2(rvec2);
     }
     err=r_norm/r_zero_norm;
#endif     
     
/* garbage collection */

     free((double *)(s->ve));
     free((VECTOR *)(s));
     
     for (i=1;i<=ngroup;i++) {
       free((double *)(AD[i]->ve));
       free((double *)(Ad[i]->ve));
       free((VECTOR *)(AD[i]));
       free((VECTOR *)(Ad[i]));           
       free((double *)(D[i]->ve));
       free((VECTOR *)(D[i])); 
     }
     
     ++ddiag;
     ++lnzv;
     free((double *)(ddiag));
     free((double *)(lnzv));
     free((double *)(pd->ve));
     free((VECTOR *)(pd));
     
     free((double *)(hatA->Valuerow));
     free((double *)(hatA->Valuecol));   
            
     counter++;
   }

   time2=clock();
   dtime=(double) (time2) - (double) (time1);
   dtime=dtime/CLOCKS_PER_SEC;
   printf("\n Converged after %d loops.",counter);   
   printf("\n Two norm of residual is: %16.14e \n",r_norm);    
   printf("\n Elapsed time: %lf seconds \n",dtime);    

/* check the final results */

   fprintf(fp1,"\n The solution vector x found: \n");
   for ( i=0;i<x->dim;i++ ) {
     fprintf(fp1,"%16.14e ",x->ve[i]);
     if ( ((i+1)%5)==0 )
       fprintf(fp1,"\n");
    }
    fprintf(fp1,"\n The final residual: \n");
    for ( i=0;i<residual->dim;i++ ) {
      fprintf(fp1,"%16.14e ",residual->ve[i]);
      if ( ((i+1)%5)==0 )
        fprintf(fp1,"\n");
    } 
    fclose(fp1);
    


   exit(0); 
}    
       
