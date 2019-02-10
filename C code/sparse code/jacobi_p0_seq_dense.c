/* jacobi_p0_seq_dense.c                                          */
/* written 25.08.1999                                             */
/* last revised 25.08.1999                                        */
/* this program is for checking the dense programming model       */
/* routines from LAPACK package are used for computing solutions  */

#include "mpi.h"
#include <unistd.h>
#include <limits.h>
#include <float.h> 
#include <time.h>
#include "datastruct_dense.h"

/* variables needed for LAPACK routines */

int lda,ldc,ldb;                         /* leading dimension of the matrix */
char uplo,trans;
double alpha,beta;
int info,incx,incy,nrhs;
   
/* variables needed for MPI */

int myid, numprocs;                         /* to identify the processors */

/**************************************************************************/
int main(argc,argv)
int argc;
char *argv[];
{  

/* local variables */
   
   static int i,j,k;
   int counter;
   int choice1,choice2;
   double verdi;
   double number,num,*ip;
   int nc,converged;
   double starttime, endtime, difftime;

   
#ifdef LINESEARCH
   double check,dotnorm;
#endif
   
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   
   sscanf(argv[1],"%d",&choice1); 
   sscanf(argv[2],"%d",&choice2); /* get the number of groups */
            
   switch(choice1) {
     case 1:
       if ((fp=fopen("dense1","r"))==NULL) {
         printf(" Error (main)        :  open \n");
         return(0);
       }
       break;
     case 2:
       if ((fp=fopen("dense2","r"))==NULL) {
         printf(" Error (main)        :  open \n");
         return(0);
       }
       break;
     case 3:
       if ((fp=fopen("dense3","r"))==NULL) {
         printf(" Error (main)        :  open \n");
         return(0);
       }
       break;
     case 4:
       if ((fp=fopen("dense4","r"))==NULL) {
         printf(" Error (main)        :  open \n");
         return(0);
       }
       break;
     case 5:
       if ((fp=fopen("dense5","r"))==NULL) {
         printf(" Error (main)        :  open \n");
         return(0);
       }
       break;
     case 6:
       if ((fp=fopen("dense6","r"))==NULL) {
         printf(" Error (main)        :  open \n");
         return(0);
       }
       break;     
   }
   epsilon=1.0e-5; 
   readdense(fp);
printf("read the matrix\n");   
printf("n %d, m %d, A[1]: %le \n",A->n,A->m,A->Value[0]);
   xexact=v_get(A->n);   
   switch(choice1) {
     case 1:
       fp=fopen("../ash219x","r");
       break;
     case 2:
       fp=fopen("../ash958x","r");
       break;
     case 3:
       fp=fopen("../ash331x","r");
       break;
     case 4:
       fp=fopen("../ash608x","r");
       break;
     case 5:
       fp=fopen("../abb313x","r");
       break; 
     case 6:
       fp=fopen("dense6x","r");
       break;     
   }
   for (i=0;i<xexact->dim;i++) {
     fscanf(fp,"%le\n",&verdi);
     xexact->ve[i]=verdi;
   }
   fclose(fp); 
printf("read xexact\n");   
   residual=v_get(A->m);
   trans='N';
   alpha=-1.0;
   beta=0.0;
   lda=A->m;
   incx=1;
   incy=1;
printf("before computing residual..\n");   
   dgemv_(&trans,&A->m,&A->n,&alpha,A->Value,&lda,xexact->ve,&incx,&beta,residual->ve,&incy);   
   x=v_get(A->n);                   /* initial x vector, all zeros */
printf("initialized residual..r[1]: %le, r[m]: %le\n",residual->ve[0],residual->ve[residual->dim-1]);
   
/* form the non-overlapping column groups of A */

   ngroup=choice2;        
   sub_ptr=CREATE(ngroup,FMATRIX *);
   ip=CREATE(1,double);
   number=(double)(A->n)/(double)(ngroup);
   num=modf(number,ip);
   if (num < 0.5)
     nc=floor(number);
   else if (num >= 0.5) 
     nc=ceil(number); 
   free((double *)(ip));
printf("before creating blocks\n");   
   sub_ptr=create_dense_blocks(A,sub_ptr,nc,ngroup);
printf("created blocks\n");   
   r_zero_norm=norm2(residual);
   printf("\n initial norm: %16.14e\n",r_zero_norm);

/* do the factorization of each block */
/* first compute A_i^T A_i */
/* then do the cholesky factorization for each of them */
/* and save the resulting R matrices */   

   Rptr=CREATE(ngroup,FMATRIX *);
printf("created Rptr..\n");   
   for (i=0;i<ngroup;i++) {
   
/* call symmetric multiplication BLAS routine */
   
     alpha=1.0;
     beta=0.0;
     uplo='L';
     trans='T';
     j=sub_ptr[i]->n;
printf("j %d\n",j);     
     B=CREATE(1,FMATRIX);
printf("created B\n");     
     B->m=j;
     B->n=j;
     numnz=B->m*B->n;
     B->Value=CREATE(numnz,double); 
printf("allocated space for B\n");         
     lda=sub_ptr[i]->m;
     k=sub_ptr[i]->m;
     ldc=sub_ptr[i]->n;
     dsyrk_(&uplo,&trans,&j,&k,&alpha,sub_ptr[i]->Value,&lda,&beta,B->Value,&ldc);
printf("computed AtransA\n");
/* do the Cholesky factorization */

     uplo='L';
     lda=B->m;
     k=B->n;
     spotf2_(&uplo,&k,B->Value,&lda,&info);
printf("did cholesky factorization..info %d\n",info);
     
     Rptr[i]=B; 
printf("R[1] %le R[last] %le\n",Rptr[i]->Value[0],Rptr[i]->Value[numnz-1]);     
   }

   free((double *)(A->Value));
   free((FMATRIX *)(A));
   
/* create hatA matrix */
   
   C=CREATE(1,FMATRIX);
   C->m=residual->dim;     
   C->n=ngroup; 
   C->Value=CREATE(C->m*C->n,double);
   CC=CREATE(1,FMATRIX); 
   CC->m=C->n;
   CC->n=C->n;
   CC->Value=CREATE(CC->m*CC->n,double);
     
   r_norm=r_zero_norm;
   err=r_norm/r_zero_norm;   
   
   converged=0;
   counter=0;

   D=CREATE(ngroup,VECTOR *);
   Ad=CREATE(ngroup,VECTOR *);

printf("\n before loop..\n");       

   starttime=MPI_Wtime();  
   while (converged==0) {
     for (i=0;i<ngroup;i++) {
       AA=sub_ptr[i];
       s=v_get(AA->n);
       alpha=1.0;
       beta=0.0;
       trans='T';
       lda=AA->m;
       incx=1;
       incy=1;
       dgemv_(&trans,&AA->m,&AA->n,&alpha,AA->Value,&lda,residual->ve,&incx,&beta,s->ve,&incy);
       B=Rptr[i];
       uplo='L';
       nrhs=1;
       lda=B->m;
       ldb=B->n;
       spotrs_(&uplo,&B->n,&nrhs,B->Value,&lda,s->ve,&ldb,&info);
       sol=v_get(AA->m);
       trans='N';
       alpha=1.0;
       beta=0.0;
       lda=AA->m;
       incx=1;
       incy=1;
       dgemv_(&trans,&AA->m,&AA->n,&alpha,AA->Value,&lda,s->ve,&incx,&beta,sol->ve,&incy);
       Ad[i]=sol;
       D[i]=s;
printf("%le %le\n",s->ve[0],s->ve[s->dim-1]);       
     }

/* put the nonzero values of hatA matrix */
     
     k=0;
     for (i=0;i<ngroup;i++) {
       for (j=0;j<C->m;j++) {
         C->Value[k]=Ad[i]->ve[j];
         k++;
       }
     }
     alpha=1.0;
     beta=0.0;
     uplo='L';
     trans='T';
     j=CC->n;     
     lda=C->m;
     k=C->m;
     ldc=j;
     dsyrk_(&uplo,&trans,&j,&k,&alpha,&C->Value,&lda,&beta,CC->Value,&ldc);
     lda=CC->n;
     j=CC->n;
     spotf2_(&uplo,&j,CC->Value,&lda,&info);
     s=v_get(C->n);
     alpha=1.0;
     beta=0.0;
     trans='T';
     lda=C->m;
     incx=1;
     incy=1;
     dgemv_(&trans,&C->m,&C->n,&alpha,C->Value,&lda,residual->ve,&incx,&beta,s->ve,&incy);
     nrhs=1;
     lda=CC->m;
     ldb=CC->n;
     spotrs_(&uplo,&CC->n,&nrhs,CC->Value,&lda,s->ve,&ldb,&info);

/* update x vector */

     k=0;
     for (i=0;i<ngroup;i++) {
       for (j=0;j<D[i]->dim;j++) {
         x->ve[k]+=D[i]->ve[j]*s->ve[i];
         k++;
       }
     }
     
/* update the residual vector */       

     for (i=0;i<ngroup;i++) {
       for (j=0;j<residual->dim;j++)
         residual->ve[j]+=Ad[i]->ve[j]*s->ve[i];
     }

/* find the error */

     r_norm=norm2(residual);
 printf("\n residual norm: %16.14e\n",r_norm);         
     
     err=r_norm/r_zero_norm;     
     
     counter++;
     if ((err < epsilon)||(r_norm==0.0)||(counter>15000)) {
       converged=1;
     }
     else {

/* garbage collection */

       free((double *)(s->ve));
       free((VECTOR *)(s));
       for (i=0;i<ngroup;i++) {
         free((double *)(D[i]->ve));
         free((double *)(Ad[i]->ve));
       }
     }     
   }
   endtime=MPI_Wtime();     
   difftime=endtime-starttime;  
   printf("\n Converged after %d loops.",counter);
   printf("\n Two norm of residual is: %16.14e \n",r_norm);
   printf("\n Elapsed time: %16.14e seconds \n",difftime); 
   MPI_Finalize();
}    
