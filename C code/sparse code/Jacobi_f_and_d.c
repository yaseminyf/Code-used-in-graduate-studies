/* last revised 13.11.1998                                             */
/* iter_seq_jacobi.c = routine for solving the LSQ problem iteratively */
/* and using jacobi updates                                            */
/* it is possible to choose either the QR or the Cholesky approach     */
/* by changing the 'define' switch                                     */
/* this version computes Jacobi with p=1, p=DS and prediction          */
/* iterations                                                          */
/***********************************************************************/
/* #include "mpi.h" */
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include "data_struct.h"

/* #define DEBUG   */

#define PONTUS  

#define CHOL  

/* #define PERM    */
/* #define MMDP    */
/* #define COLP      */
/* #define NEWPERM    */
  
VECTOR **D;
MATRIX *hatA,*tildeA,**tildeAptr;
VECTOR **Apptr;
VECTOR *pvec,**pvecptr,*pd;
MATRIX *AA;

/* #define PISONES      */

 #define XERROR   /* to check convergence with error on x */

#define FMSTART   /* starting values of p are F_M values */

#ifdef FMSTART
  VECTOR *evec;
#endif
  
#define PREDICT        

#ifdef PREDICT
  VECTOR *vvec,*zvec,**zvecptr;
  int antall,t;
  double **lnzvptr,**ddiagptr;
#endif
  
#ifdef XERROR
  VECTOR *xnormv;
#endif
  
#ifdef PERM
  int *p,*ivp;
  double *tmp;
#endif
 
#define DIFFGROUPSIZE        
                                     
/**************************************************************************/
int main()

{  

/* local variables */
   
   static int i,j,k,l;
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

#ifdef PREDICT
   printf("\n Enter the number of prediction iterations..");
   scanf("%d",&antall);
#endif
      
         
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
      case '6':
        if ((fp=fopen("matrix6","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case '7':
        if ((fp=fopen("matrix7","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case '8':
        if ((fp=fopen("matrix8","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case '9':
        if ((fp=fopen("matrix9","r"))==NULL) {
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

#ifdef DEBUG    
      fprintf(fp1,"\n got the xexact vector..");
#endif

#ifdef PONTUS
/*      xexact=vones(xexact); */
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
        fscanf(fp,"%le\n",&verdi,&c1);
        xexact->ve[i]=verdi;
      }
      fclose(fp);

#ifdef DEBUG    
      fprintf(fp1,"\n initialized the xexact vector..");
#endif   

      fp1=fopen("data","w");

#ifdef PERM
/* permute the matrix A */

      switch(choice1) {
        case '1':
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
        case '2':
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
        case '3':
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
        case '4':
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
        case '5':
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
        fscanf(fp,"%d\n",&j,&c1);
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
      
#ifdef DEBUG    
      fprintf(fp1,"\n got the rhs vector...");
#endif 

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

/* garbage collection */

#ifndef XERROR
      free((double *)(xexact->ve));
      free((VECTOR *)(xexact));      
#endif
      
#ifdef DEBUG    
      fprintf(fp1,"\n initialized the rhs vector..");
#endif   
   }
   else if ( (choice1 >= '6') && (choice1<= '9') ) {
      epsilon=1.0e-01;
      fp1=fopen("data","w");
      readnewmatrix(fp,fp1);
   }

#ifdef DEBUG
   for (i=1;i<=20;i++) 
     fprintf(fp1,"%16.15e \n",A->Valuecol[i]);    
#endif

   x=v_get(A->n);

#ifdef DEBUG    
   fprintf(fp1,"\n got the initial solution x..");
#endif   

   residual=v_get(rhs->dim);

#ifdef DEBUG    
   fprintf(fp1,"\n got the residual vector..");
#endif   

   for ( i=0;i<rhs->dim;i++ )
     residual->ve[i]=-1.0*rhs->ve[i];

/* garbage collection */

   free((double *)(rhs->ve));
   free((VECTOR *)(rhs));
   
#ifdef DEBUG    
   fprintf(fp1,"\n computed the residual..");
#endif   

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
     free((MATRIX *)(A)); 
     
     r_zero_norm=norm2(residual);
     printf("\n initial norm: %16.14e",r_zero_norm);
   }
   else if ( (choice1 >= '6') && (choice1<= '9') ){

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
     sub_ptr=create_nover25(A,sub_ptr,ngroup,zeroptr);
     nc=25;
#endif

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
   
#ifdef DEBUG    
   fprintf(fp1,"\n sparspak initialization ends..");
#endif    

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
printf("\n %lf",pvec->ve[j]);
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

   for ( i=0; i<ngroup; i++ ) {
     free((double *)(Apptr[i]->ve));
     free((VECTOR *)(Apptr[i]));
   }
   free((VECTOR **)(Apptr));  
              
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

#ifdef DEBUG    
   fprintf(fp1,"\n T vector is allocated..");
#endif   

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

#ifdef DEBUG       
     fprintf(fp1,"\n after inxywb...");   
     fprintf(fp1,"\n before orcolb...");
#endif

     orcolb(T);

#ifdef DEBUG       
     fprintf(fp1,"\n after orcolb...");

#endif   

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

#ifdef DEBUG    
   fprintf(fp1,"\n before solution...");
#endif

#ifdef XERROR
   xnormv=v_get(xexact->dim);
   for (i=0;i<xnormv->dim;i++)
     xnormv->ve[i]=x->ve[i]-xexact->ve[i];
   err=norm2(xnormv);
#else
   r_norm=r_zero_norm;
   err=r_norm/r_zero_norm;   
#endif      

   omega=1.0;
   counter=0;

#ifdef PREDICT
   lnzvptr=CREATE(ngroup,double *);
   ddiagptr=CREATE(ngroup,double *);
   --lnzvptr;
   --ddiagptr;
#endif
   
   D=CREATE(ngroup,VECTOR *);
   --D;
printf("\n before loop...");   
   time1=clock();
   
   while (err >= epsilon) {
#ifdef PREDICT
     if ( counter > 0 ) {
 /*      vvec=v_get(residual->dim);
       for (i=0;i<ngroup;i++)
         vvec=m_fulvectormlt(sub_ptr[i],zvecptr[i],vvec);
       for (j=0;j<vvec->dim;j++)
         vvec->ve[j]+=residual->ve[j]; */
       for (t=0;t<antall;t++) {
         for (i=1;i<=ngroup;i++) {
           AA=tildeAptr[i-1];
           s=v_get(AA->n);
           sol=v_get(AA->m);
           s=multAtransr3(AA,vvec,s);

           s=permdv(AA->n,s,invpptr[i]);
           s=gsslv(AA->n,xrnzptr[i],lnzvptr[i],xnzsubptr[i],nzsubptr[i],ddiagptr[i],s);     
           s=permdv(AA->n,s,permptr[i]);   

           D[i]=s;
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
               for (l=1;l<=sub_ptr[j-1]->n;l++) {
                 pd->ve[k]+=D[i]->ve[count];
                 count++;
                 k++;
               }
             }
           }
         }
         createC2(hatA,pd,sub_ptr,ngroup);  
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
         s=multAtransr3(hatA,vvec,s);
         s=permdv(hatA->n,s,invpv);
         s=gsslv(hatA->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,s);   
         s=permdv(hatA->n,s,permv);  
         k=0;
         for (i=1;i<=ngroup;i++) {
           for (j=1;j<=sub_ptr[i-1]->n;j++) {
             zvecptr[i-1]->ve[j-1]+=pd->ve[k]*s->ve[i-1];
             k++;
           }
         }
         sol=v_get(vvec->dim);
         sol=m_fulvectormlt(hatA,s,sol);
         for ( j=0;j<vvec->dim;j++ ) {
           vvec->ve[j]+=sol->ve[j];
         }
         free((double *)(s->ve));
         free((VECTOR *)(s));
         free((double *)(sol->ve));
         free((VECTOR *)(sol));
         ++ddiag;
         ++lnzv;
         free((double *)(ddiag));
         free((double *)(lnzv));
         free((double *)(pd->ve));
         free((VECTOR *)(pd));
       }
       free((double *)(vvec->ve));
       free((VECTOR *)(vvec));
       for (i=0;i<ngroup;i++) {
         free((double *)(pvecptr[i]->ve));
         free((VECTOR *)(pvecptr[i]));
         pvecptr[i]=zvecptr[i];
       }
     }
#endif   

     if ( counter > 0 ) {
#ifndef PISONES     
       putvaluesintildeA(tildeAptr,sub_ptr,pvecptr,ngroup);  
#endif 
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

#ifdef PREDICT
       lnzvptr[i]=lnzv;
       ddiagptr[i]=ddiag;
#endif
                       
       s=v_get(AA->n);
       sol=v_get(AA->m);
       s=multAtransr3(AA,residual,s);

       s=permdv(AA->n,s,invpptr[i]);
       s=gsslv(AA->n,xrnzptr[i],lnzv,xnzsubptr[i],nzsubptr[i],ddiag,s);
         
/* permute the solution to its original order of colums */

       s=permdv(AA->n,s,permptr[i]);   
    
       D[i]=s;
   
#ifdef DEBUG
       fprintf(fp1,"\n computed sol..");
#endif  

#ifndef PREDICT
       ++lnzv;
       ++ddiag;
       free((double *)(lnzv));
       free((double *)(ddiag));
#endif       
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

/* put in the values of hatA */

     createC2(hatA,pd,sub_ptr,ngroup);  
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

/* update x vector and the pvecptr at the same time */
     
#ifndef PISONES
#ifndef PREDICT
     for (i=0;i<ngroup;i++) {
       free((double *)(pvecptr[i]->ve));
       free((VECTOR *)(pvecptr[i]));
     }
     free((VECTOR **)(pvecptr));
     
     pvecptr=CREATE(ngroup,VECTOR *);
#endif     
#endif

#ifdef PREDICT
     zvecptr=CREATE(ngroup,VECTOR *);   
#endif
       
     k=0;
     for (i=1;i<=ngroup;i++) {
#ifndef PISONES 
#ifdef PREDICT
       zvec=v_get(sub_ptr[i-1]->n);
#else
       pvec=v_get(sub_ptr[i-1]->n);       
#endif       
#endif       
       for (j=1;j<=sub_ptr[i-1]->n;j++) {
         x->ve[k]+=pd->ve[k]*s->ve[i-1];
#ifndef PISONES 
#ifdef PREDICT
         zvec->ve[j-1]=pd->ve[k]*s->ve[i-1];
#else
         pvec->ve[j-1]=pd->ve[k]*s->ve[i-1];
#endif

#endif  
         k++;       
       }
#ifndef PISONES
#ifdef PREDICT
       zvecptr[i-1]=zvec;
#else
       pvecptr[i-1]=pvec;       
#endif       
#endif       
     }
     
/* update the residual vector */

     sol=v_get(residual->dim);
printf("\n");    
for (i=0;i<s->dim;i++)
printf("\n %16.14e",s->ve[i]);     
     sol=m_fulvectormlt(hatA,s,sol);
#ifdef PREDICT
     vvec=v_get(residual->dim);
#endif          
     for ( j=0;j<residual->dim;j++ ) {
       residual->ve[j]+=sol->ve[j];
#ifdef PREDICT
       vvec->ve[j]=residual->ve[j]+sol->ve[j];
#endif              
     }
     
/* garbage collection */

     free((double *)(s->ve));
     free((VECTOR *)(s));
     if ( choice1 > '5' ) {     
       free((double *)(sol->ve));
       free((VECTOR *)(sol));
     }
     ++ddiag;
     ++lnzv;
     free((double *)(ddiag));
     free((double *)(lnzv));
     free((double *)(pd->ve));
     free((VECTOR *)(pd));
     
     free((double *)(hatA->Valuerow));
     free((double *)(hatA->Valuecol));   
       
/* find the error */

#ifdef XERROR
     for (i=0;i<xnormv->dim;i++)
       xnormv->ve[i]=x->ve[i]-xexact->ve[i];
     err=norm2(xnormv);
     r_norm=norm2(residual);
printf("\n error: %16.14e",err);  
printf("\n residual norm: %16.14e",r_norm);        
#else
     if ( choice1 <= '5' ) {
       r_norm=norm2(residual);
printf("\n residual norm: %16.14e",r_norm);       
       
     }
     else {
       rvec2=multAtransr2(A,sol,rvec2);
       free((double *)(sol->ve));
       free((VECTOR *)(sol));
       r_norm=norm2(rvec2);
     }
     err=r_norm/r_zero_norm;
#endif     
           
     counter++;
   }

   time2=clock();
   dtime=(double) (time2) - (double) (time1);
   dtime=dtime/CLOCKS_PER_SEC;
   printf("\n Converged after %d loops.",counter);

#ifdef XERROR
   printf("\n The error is: %16.14e \n",err   );
#else   
   printf("\n Two norm of residual is: %16.14e \n",r_norm); 
#endif
   
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
    
/* garbage collection */

   ++xrnzv;
   ++nzsubv;
   ++xnzsubv;
   ++permv;
   ++invpv;
   ++ddiagptr;
   ++xrnzptr;
   ++nzsubptr;
   ++xnzsubptr;
   ++permptr;
   ++invpptr;
   ++rnzptr;

   free((int *)(hatA->colind));
   free((int *)(hatA->rowptr));
   free((int *)(hatA->Iind));
   free((int *)(hatA->Jind));
   free((MATRIX *)(hatA));   
   free((int *)(xrnzv));
   free((int *)(nzsubv));
   free((int *)(xnzsubv));
   free((int *)(permv));
   free((int *)(invpv));
   free((double **)(ddiagptr));
   free((int **)(xrnzptr));
   free((int **)(nzsubptr));
   free((int **)(xnzsubptr));
   free((int **)(permptr));
   free((int **)(invpptr));
   free((double **)(rnzptr));

   exit(0); 
}    
       
