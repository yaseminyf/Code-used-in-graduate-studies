/* jacobi_p0_seq_dense.c                                          */
/* written 25.08.1999                                             */
/* last revised 14.09.1999                                        */
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
   double starttime, endtime, totaltime;

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
   residual=v_get(A->m);
   trans='N';
   alpha=-1.0;
   beta=0.0;
   lda=A->m;
   incx=1;
   incy=1;
   dgemv_(&trans,&A->m,&A->n,&alpha,A->Value,&lda,xexact->ve,&incx,&beta,residual->ve,&incy);   
   x=v_get(A->n);                   /* initial x vector, all zeros */
   
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
   sub_ptr=create_dense_blocks(A,sub_ptr,nc,ngroup);
   r_zero_norm=norm2(residual);
   printf("\n initial norm: %16.14e\n",r_zero_norm);

   totaltime=0.0;
   for (i=0;i<ngroup;i++) {
       AA=sub_ptr[i];
       s=v_get(AA->n);
       alpha=-1.0;
       beta=0.0;
       trans='T';
       lda=AA->m;
       incx=1;
       incy=1;
       starttime=MPI_Wtime();
       dgemv_(&trans,&AA->m,&AA->n,&alpha,AA->Value,&lda,residual->ve,&incx,&beta,s->ve,&incy);
       endtime=MPI_Wtime();
       totaltime+=endtime-starttime;
     }

   printf("\n Elapsed time: %16.14e seconds \n",totaltime); 
   MPI_Finalize();
}    
