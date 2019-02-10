/* last revised 20.09.1998                                         */
/* iter_seq_GS.c = routine for solving the LSQ problem iteratively */
/* and using GS updates                                            */
/* it is possible to choose either the QR or the Cholesky approach */
/* by changing the 'define' switch in 'data_struct.h'              */
/*******************************************************************/
/* #include "mpi.h" */
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include "data_struct.h"

/* #define DEBUG   */
#define CHOL  
/* #define QR   */

#define PERM  
/* #define MMDP  */ 
/* #define COLP     */
#define NEWPERM  

#ifdef PERM
  int *p,*ivp;
  double *tmp;
#endif
 
 #define DIFFGROUPSIZE  

/* int myid, numprocs;    */                     /* to identify the processors */
                                     
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
   double time1,time2;
   double dtime;
   
#ifdef DIFFGROUPSIZE   
   int choice2;
#endif

   printf("\n Enter the number of the matrix you want to use..");
   printf("\n              [1,6]\n");
   choice1=getchar(); 

#ifdef DIFFGROUPSIZE   
   printf("\n Enter the number of groups..");
   scanf("%d",&choice2); 
#endif
   
/*   choice1='1'; */
         
   switch(choice1) {
      case '1':
        if ((fp=fopen("../matrix1","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        if ((fp1=fopen("matrix1dat","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case '2':
        if ((fp=fopen("../matrix2","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        if ((fp1=fopen("matrix2dat","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case '3':
        if ((fp=fopen("../matrix3","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        if ((fp1=fopen("matrix3dat","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case '4':
        if ((fp=fopen("../matrix4","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        if ((fp1=fopen("matrix4dat","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case '5':
        if ((fp=fopen("../matrix5","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        if ((fp1=fopen("matrix5dat","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case '6':
        if ((fp=fopen("../matrix6","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case '7':
        if ((fp=fopen("../matrix7","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case '8':
        if ((fp=fopen("../matrix8","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
      case '9':
        if ((fp=fopen("../matrix9","r"))==NULL) {
          printf(" Error (main)        :  open \n");
          return(0);
        }
        break;
   }
   if ( (choice1 >= '1') && (choice1<= '5') ){
      epsilon=1.0e-5;
      readmatrix(fp,fp1); 
      xexact=v_get(A->n);

#ifdef DEBUG    
      fprintf(fp1,"\n got the xexact vector..");
#endif   
      
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
      free(xexact->ve);
      xexact->ve=tmp;
#endif
      
      rhs=m_fulvectormlt(A,xexact,rhs); 
      
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
   sub_ptr=CREATE(ngroup,int);
   zeroptr=CREATE(ngroup,int);
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
#else      
     sub_ptr=create_nover(A,sub_ptr,ngroup,zeroptr);
     nc=10;
#endif
     r_zero_norm=norm2(residual);
     printf("\n initial norm: %lf",r_zero_norm);
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
#else  
     sub_ptr=create_nover25(A,sub_ptr,ngroup,zeroptr);
     nc=25;
#endif
     rvec2=v_get(A->n);
     rvec2=multAtransr2(A,residual,rvec2);
     r_zero_norm=norm2(rvec2);
   }

/* allocate the necessary storage to keep the vectors for each group */

   rnzptr=CREATE(ngroup,double);
   ddiagptr=CREATE(ngroup,double);
   xrnzptr=CREATE(ngroup,int);
   nzsubptr=CREATE(ngroup,int);
   xnzsubptr=CREATE(ngroup,int);
   permptr=CREATE(ngroup,int);
   invpptr=CREATE(ngroup,int);
   --rnzptr;
   --ddiagptr;
   --xrnzptr;
   --nzsubptr;
   --xnzsubptr;
   --permptr;
   --invpptr;
   
/* initialize for SPARSPAK */

#ifdef DEBUG    
   fprintf(fp1,"\n sparspak initialization begins..");
#endif   

   spbusr_.MAXSB=4000;
   TOL=0.0000000001;
   TYPTOL=1;

   spksys_.MAXINT=SHRT_MAX;
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
    
   spbusr_.IERRB=0;
   spbusr_.MSGLVB=1;
   
#ifdef DEBUG    
   fprintf(fp1,"\n sparspak initialization ends..");
#endif    

/* for each group do the factorization and store the necessary vectors */
    
   for ( l=0; l<ngroup; l++ ) {
     A=sub_ptr[l];
     spbusr_.MCOLS=A->n;
     spbusr_.MSEQNS=A->m-zeroptr[l][0];
    
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
    
     for ( j=1;j<=A->m;j++ ) {
       if ( A->rowptr[j] != 0 ) { 
         count=1;
         while ( A->rowptr[j+count] == 0 )
           count++;
         NSUBS=A->rowptr[j+count]-A->rowptr[j];
         spbcon_.NOFNZ+=NSUBS;
         ROWNUM=rnum;
         rnum++;
         SUBS=CREATE(NSUBS,int);
         for ( k=0;k<NSUBS;k++ ) {
           *(SUBS+k)=A->colind[A->rowptr[j]+k];
         }
         --SUBS;
         inxywb(ROWNUM,TYPE,NSUBS,SUBS, T);
         SUBS++;
         free(SUBS); 
       }
     }

#ifdef DEBUG       
     fprintf(fp1,"\n after inxywb...");   
     fprintf(fp1,"\n before orcolb...");
#endif
  
     orcolb(T);
   
#ifdef DEBUG       
     fprintf(fp1,"\n after orcolb...");

#ifdef QR   
     fprintf(fp1,"\n before lsqslvsn3...");
#endif   

#endif   

#ifdef QR
     xrnzv=CREATE(spbcon_.NCOLS+1,int);
     xnzsubv=CREATE(spbcon_.NCOLS,int);
     permv=CREATE(spbcon_.NCOLS,int);
     invpv=CREATE(spbcon_.NCOLS,int);
     nzsubv=CREATE(spbcon_.NOFNZ,int);
     --xrnzv;
     --xnzsubv;
     --permv;
     --invpv;
     --nzsubv;

     ddiag=CREATE(spbcon_.NCOLS,double);
     rnzv=CREATE(spbcon_.NOFNZ,double);
     --ddiag;
     --rnzv;
     
     zerorv(spbcon_.NCOLS, ddiag);
     zerols(spbcon_.NCOLS, T, xrnz, rnzv);
     
     lsqslvsn3(TOL,TYPTOL,T,xnzsub,nzsub,xrnz,invp,perm);

/* garbage collection */

     T++;
     SUBS++;
     VALUES++;
     free(T);
     free(SUBS);
     free(VALUES);
     
#ifdef DEBUG 
     fprintf(fp1,"\n after lsqslvsn3...");
#endif     

#endif

#ifdef CHOL
 
/* compute A^T*A and initialize lnz and diag */

     ddiag=CREATE(spbcon_.NCOLS,double);
     lnzv=CREATE(spbcon_.NOFNZ,double);
     --ddiag;
     --lnzv;
     for (j=1;j<=spbcon_.NOFNZ;j++)
       lnzv[j]=0.0;
/*    zerols(spbcon_.NOFNZ, T, xrnz, lnzv); */
     for (j=1;j<=spbcon_.NCOLS;j++)
       ddiag[j]=0.0;
/*       zerorv(spbcon_.NCOLS, ddiag); */
     
#ifdef DEBUG
     fprintf(fp1,"\n before sym_mult......");
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

/* get the necessary storage for A^TA */

     i=spbcon_.NOFNZ+spbcon_.NCOLS;
     isub=CREATE(i,int);
     jsub=CREATE(i,int);
     values=CREATE(i,double);
     --isub;
     --jsub;
     --values;
     sym_mult(A,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values);

#ifdef DEBUG
     fprintf(fp1,"\n after sym_mult......");
#endif      

/* the next three vectors are initialized in gsfct */
   
     llink=CREATE(spbcon_.NCOLS,int);
     temp=CREATE(spbcon_.NCOLS,double);
     first=CREATE(spbcon_.NCOLS,int);
     --llink;
     --temp;
     --first;

/* garbage collection */

     T++;
     isub++;
     jsub++;
     values++;
     free(T);
     free(isub);
     free(jsub);
     free(values);
   
/* begin Cholesky factorization */ 

#ifdef DEBUG
     fprintf(fp1,"\n before gsfct..");
#endif

     gsfct(spbcon_.NCOLS,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,llink,first,temp,flag__); 

#ifdef DEBUG
     fprintf(fp1,"\n after gsfct..");
#endif

#endif
   
     ddiagptr[l+1]=ddiag;
     xrnzptr[l+1]=xrnzv;
     nzsubptr[l+1]=nzsubv;
     xnzsubptr[l+1]=xnzsubv;
     permptr[l+1]=permv;
     invpptr[l+1]=invpv;

#ifdef QR
     rnzptr[l+1]=rnzv;
#endif

#ifdef CHOL
     rnzptr[l+1]=lnzv;
#endif
     
   } /* end of for loop for all groups */

/* begin the solution */

#ifdef DEBUG    
   fprintf(fp1,"\n before solution...");
#endif

   i=1;
   A=sub_ptr[i-1];
   s=v_get(A->n);
   sol=v_get(A->m);
   r_norm=r_zero_norm;
   err=r_norm/r_zero_norm;
   omega=1.0;
   counter=0;
/*   time1=MPI_Wtime(); */
   while (err >= epsilon) {
     s=multAtransr3(A,residual,s);

#ifdef QR   
     s=solution(A->n,s,rnzptr[i],ddiagptr[i],xnzsubptr[i],nzsubptr[i],xrnzptr[i],permptr[i],invpptr[i]);
#endif

#ifdef CHOL

     s=permdv(A->n,s,invpptr[i]);
     s=gsslv(A->n,xrnzptr[i],rnzptr[i],xnzsubptr[i],nzsubptr[i],ddiagptr[i],s);
       
/* permute the solution to its original order of colums */

     s=permdv(A->n,s,permptr[i]);   
    
#endif   

     sol=m_fulvectormlt(A,s,sol);
   
#ifdef DEBUG
     fprintf(fp1,"\n computed sol..");
#endif  

     for ( j=0;j<s->dim;j++ ) {       
       x->ve[(i-1)*nc+j]=x->ve[(i-1)*nc+j]+omega*s->ve[j];
     } 
     
     for ( j=0;j<residual->dim;j++ ) {
       if (sol->ve[j]!=0.0) {
         residual->ve[j]=residual->ve[j]+omega*sol->ve[j];
       }
     }   

/* find the norm of the residual */

     if ( choice1 <= '5' ) {
       r_norm=norm2(residual);
     }
     else {
       rvec2=multAtransr2(A,sol,rvec2);
       r_norm=norm2(rvec2);
     }
     err=r_norm/r_zero_norm;
     i++;
     if (i==ngroup+1) i=1;
     A=sub_ptr[i-1];
     s=v_get(A->n);
     sol=v_get(A->m);
     counter++;
   }
/*   time2=MPI_Wtime();
   dtime=time2-time1; */
   printf("\n Converged after %d loops.",counter);
   printf("\n Two norm of residual is: %16.14e \n",r_norm); 
/*   printf("\n Elapsed time: %16.14e seconds\n",dtime);    */
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

   ++ddiagptr;
   ++xrnzptr;
   ++nzsubptr;
   ++xnzsubptr;
   ++permptr;
   ++invpptr;
   ++rnzptr;

   free(ddiagptr);
   free(xrnzptr);
   free(nzsubptr);
   free(xnzsubptr);
   free(permptr);
   free(invpptr);
   free(rnzptr);

   exit(0); 
}    
       
