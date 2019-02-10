/* jacobi_predict_seq_dense_with_comm.c                           */
/* written 15.09.1999                                             */
/* last revised 15.09.1999                                        */
/* this program is for checking the dense programming model       */
/* routines from LAPACK package are used for computing solutions  */
/* this sequential version includes also communication            */
/* in this one predictions                                        */

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
int dest, tag;                   /* node to send/receive msg, and msg tag */
int buf_count;                       /* the number of elements to be sent */
double *sndvector;
double **sndptr;
int position;                                      
char *buffer; 

/**************************************************************************/
int main(argc,argv)
int argc;
char *argv[];
{  

/* local variables */
   
   static int i,j,k,t,l;
   int counter;
   int choice1,choice2,choice3,count;
   double verdi,*deltav;
   double number,num,*ip;
   int nc,converged,sind,eind;
   double starttime, endtime, difftime;
   int buf_countmain;
   int sfdouble;
   
   MPI_Request *req,*req2;
   MPI_Status status;
   
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   
   sscanf(argv[1],"%d",&choice1); 
   sscanf(argv[2],"%d",&choice2); /* get the number of groups */
   sscanf(argv[3],"%d",&choice3); /* get the number of predictions */

   sfdouble=sizeof(double);
        
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
     case 7:
       if ((fp=fopen("dense7","r"))==NULL) {
         printf(" Error (main)        :  open \n");
         return(0);
       }
       break; 
     case 8:
       if ((fp=fopen("dense8","r"))==NULL) {
         printf(" Error (main)        :  open \n");
         return(0);
       }
       break; 
     case 9:
       if ((fp=fopen("dense9","r"))==NULL) {
         printf(" Error (main)        :  open \n");
         return(0);
       }
       break; 
     case 10:
       if ((fp=fopen("dense10","r"))==NULL) {
         printf(" Error (main)        :  open \n");
         return(0);
       }
       break; 
     case 11:
       if ((fp=fopen("dense11","r"))==NULL) {
         printf(" Error (main)        :  open \n");
         return(0);
       }
       break; 
        case 12:
         if ((fp=fopen("dense12","r"))==NULL) {
           printf(" Error (main)        :  open \n");
           return(0);
         }
         break; 
       case 13:
         if ((fp=fopen("dense13","r"))==NULL) {
           printf(" Error (main)        :  open \n");
           return(0);
         }
         break;
       case 14:
         if ((fp=fopen("dense14","r"))==NULL) {
           printf(" Error (main)        :  open \n");
           return(0);
         }
         break;
       case 15:
         if ((fp=fopen("dense15","r"))==NULL) {
           printf(" Error (main)        :  open \n");
           return(0);
         }
         break;
       case 16:
         if ((fp=fopen("dense16","r"))==NULL) {
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
     case 7:
       fp=fopen("dense7x","r");
       break; 
     case 8:
       fp=fopen("dense8x","r");
       break; 
     case 9:
       fp=fopen("dense9x","r");
       break; 
     case 10:
       fp=fopen("dense10x","r");
       break; 
     case 11:
       fp=fopen("dense11x","r");
       break;  
      case 12:
         fp=fopen("dense12x","r");
         break;
       case 13:
         fp=fopen("dense13x","r");
         break;
       case 14:
         fp=fopen("dense14x","r");
         break;
       case 15:
         fp=fopen("dense15x","r");
         break;
       case 16:
         fp=fopen("dense16x","r");
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
   fp=fopen("res_predict_dense_seq","a");
   fprintf(fp,"\n file: %d , m: %d , n: %d ",choice1,A->m,A->n);
   fprintf(fp,"\n initial norm: %16.14e",r_zero_norm);

/* form A-i*p-i  */

   pvecptr=CREATE(ngroup,VECTOR *);
   for ( l=0; l<ngroup; l++ ) {
     AA=sub_ptr[l];
     pvec=v_get(AA->n);
     pvec=vones(pvec);
     
#ifdef FMSTART
     evec=v_get(AA->n);
     evec=vones(evec);
     sol=v_get(AA->m);
     trans='N';
     alpha=1.0;
     beta=0.0;
     lda=AA->m;
     incx=1;
     incy=1;
     dgemv_(&trans,&AA->m,&AA->n,&alpha,AA->Value,&lda,evec->ve,&incx,&beta,sol->ve,&incy);   
     s=v_get(AA->n);
     trans='T';
     dgemv_(&trans,&AA->m,&AA->n,&alpha,AA->Value,&lda,sol->ve,&incx,&beta,s->ve,&incy);   
     for (j=0;j<pvec->dim;j++)
       pvec->ve[j]=pvec->ve[j]/s->ve[j];    
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
   for (i=0;i<ngroup;i++) {
     AA=sub_ptr[i];
     sol=v_get(AA->m);
     trans='N';
     alpha=1.0;
     beta=0.0;
     lda=AA->m;
     incx=1;
     incy=1;
     dgemv_(&trans,&AA->m,&AA->n,&alpha,AA->Value,&lda,pvecptr[i]->ve,&incx,&beta,sol->ve,&incy);   
     Apptr[i]=sol;
   }

/* form the tildeA matrices */ 

   tildeAptr=CREATE(ngroup,FMATRIX *);
   for (i=0;i<ngroup;i++) {
     AA=sub_ptr[i];
     tildeA=CREATE(1,FMATRIX);
     tildeA=createtildeA(tildeA,AA,Apptr,ngroup,i); 
     tildeAptr[i]=tildeA;   
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

   Rptr=CREATE(ngroup,FMATRIX *);  
   tmpApptr=CREATE(ngroup,VECTOR *);
   
   r_norm=r_zero_norm;
   err=r_norm/r_zero_norm;   
   
   converged=0;
   counter=0;

   D=CREATE(ngroup,VECTOR *);
   Ad=CREATE(ngroup,VECTOR *);
   --D;
   --Ad;
   for (i=1;i<=ngroup;i++) {
     D[i]=v_get(sub_ptr[i-1]->n);
     Ad[i]=v_get(residual->dim);
   }

   buf_count=nc+ngroup+residual->dim;
   sndptr=CREATE(ngroup,double *);
   --sndptr;
   for (i=1;i<=ngroup;i++) {
     sndvector=CREATE(buf_count,double);
     sndptr[i]=sndvector;
   }
   buf_countmain=nc+ngroup+residual->dim;
   sndvector=CREATE(buf_countmain,double); 
   req=CREATE(ngroup,MPI_Request);
   req2=CREATE(ngroup,MPI_Request);       
   starttime=MPI_Wtime();  
   while (converged==0) {
     for (i=1;i<=ngroup;i++ ) {
       tag=i*100;
       MPI_Irecv(sndptr[i],buf_countmain,MPI_DOUBLE,MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,&req[i-1]);           
     }         
     deltav=CREATE(ngroup,double);
     --deltav;
     for (i=1;i<=ngroup;i++) {
       tildeA=tildeAptr[i-1];
       AA=sub_ptr[i-1];

/* call symmetric multiplication BLAS routine */
   
       alpha=1.0;
       beta=0.0;
       uplo='L';
       trans='T';
       j=tildeA->n;
       B=CREATE(1,FMATRIX);
       B->m=j;
       B->n=j;
       numnz=B->m*B->n;
       B->Value=CREATE(numnz,double); 
       lda=tildeA->m;
       k=tildeA->m;
       ldc=tildeA->n;
       dsyrk_(&uplo,&trans,&j,&k,&alpha,tildeA->Value,&lda,&beta,B->Value,&ldc);

/* do the Cholesky factorization */

       uplo='L';
       lda=B->m;
       k=B->n;
       dpotf2_(&uplo,&k,B->Value,&lda,&info);
       Rptr[i-1]=B;
     }
     for (i=1;i<=ngroup;i++) {
       tildeA=tildeAptr[i-1];
       AA=sub_ptr[i-1]; 
       B=Rptr[i-1];    
       s=v_get(tildeA->n);
       alpha=-1.0;
       beta=0.0;
       trans='T';
       lda=tildeA->m;
       incx=1;
       incy=1;
       dgemv_(&trans,&tildeA->m,&tildeA->n,&alpha,tildeA->Value,&lda,residual->ve,&incx,&beta,s->ve,&incy);
       uplo='L';
       nrhs=1;
       lda=B->m;
       ldb=B->n;
       dpotrs_(&uplo,&B->n,&nrhs,B->Value,&lda,s->ve,&ldb,&info);

/* take the littles out of s */

       littles=v_get(AA->n);
       k=AA->n;
       for (j=0;j<k;j++) {
         littles->ve[j]=s->ve[i+j-1];
       }                      
       sol=v_get(AA->m);
       trans='N';
       alpha=1.0;
       beta=0.0;
       lda=AA->m;
       incx=1;
       incy=1;
       dgemv_(&trans,&AA->m,&AA->n,&alpha,AA->Value,&lda,littles->ve,&incx,&beta,sol->ve,&incy);

/* send both s and sol to master */

       k=s->dim;
       for (j=0;j<k;j++)
         sndvector[j]=s->ve[j];
       count=sol->dim;
       for (j=0;j<count;j++) 
         sndvector[j+k]=sol->ve[j];      
       buf_count=k+count;
       tag=i*100;
       MPI_Isend(&sndvector[0],buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&req2[i-1]); 
       free((double *)(s->ve));
       free((double *)(sol->ve));
       free((VECTOR *)(s));
       free((VECTOR *)(sol)); 
       free((double *)(littles->ve));
       free((VECTOR *)(littles));
       MPI_Request_free(&req2[i-1]);       
     }
     for (i=1;i<=ngroup;i++ ) {
       MPI_Waitany(ngroup,req,&j,&status);
       dest=j+1;
       for (j=1;j<dest;j++)
         deltav[j]+=sndptr[dest][j-1];
       count=D[dest]->dim;
       for (j=0;j<count;j++) 
         D[dest]->ve[j]=sndptr[dest][j+dest-1];
       count+=dest;
       l=D[dest]->dim+ngroup-1;
       k=1;
       for (j=count;j<=l;j++) {
         deltav[k+dest]+=sndptr[dest][j-1];
         k++;
       }
       sind=D[dest]->dim+ngroup-1;
       eind=Ad[dest]->dim;
       for (j=0;j<eind;j++)
         Ad[dest]->ve[j]=sndptr[dest][j+sind];
     }
     
/* compute the columns of A*D=hatA matrix */
   
     for (i=1;i<=ngroup;i++) {
       eind=Ad[i]->dim;
       for (j=0;j<eind;j++) {
         Ad[i]->ve[j]+=Apptr[i-1]->ve[j]*deltav[i];
       }
     } 
     
/* put the nonzero values of hatA matrix */
     
     k=0;
     for (i=1;i<=ngroup;i++) {
       count=C->m;
       for (j=0;j<count;j++) {
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
     dsyrk_(&uplo,&trans,&j,&k,&alpha,C->Value,&lda,&beta,CC->Value,&ldc);
     lda=CC->n;
     j=CC->n;
     dpotf2_(&uplo,&j,CC->Value,&lda,&info);
     s=v_get(C->n);
     alpha=-1.0;
     beta=0.0;
     trans='T';
     lda=C->m;
     incx=1;
     incy=1;
     dgemv_(&trans,&C->m,&C->n,&alpha,C->Value,&lda,residual->ve,&incx,&beta,s->ve,&incy);
     nrhs=1;
     lda=CC->m;
     ldb=CC->n;
     dpotrs_(&uplo,&CC->n,&nrhs,CC->Value,&lda,s->ve,&ldb,&info);

/* update the residual vector, vvec and tmpApptr at the same time */     
 
     for (i=0;i<ngroup;i++) 
       tmpApptr[i]=v_get(residual->dim);
     vvec=v_get(residual->dim);
     for (i=0;i<ngroup;i++) {
       count=residual->dim;
       for ( j=0;j<count;j++ ) {
         verdi=s->ve[i]*Ad[i+1]->ve[j];
         residual->ve[j]+=verdi;
         vvec->ve[j]+=verdi;
         tmpApptr[i]->ve[j]+=verdi;
       }
     }
     k=residual->dim;
     for (i=0;i<k;i++)
       vvec->ve[i]+=residual->ve[i];
       
/* update x vector and zvec */

     zvec=v_get(x->dim);
     k=0;
     for (i=0;i<ngroup;i++) {
       count=D[i+1]->dim;
       for (j=0;j<count;j++) {
         verdi=s->ve[i]*(D[i+1]->ve[j]+pvecptr[i]->ve[j]*deltav[i+1]);
         x->ve[k]+=verdi;
         zvec->ve[k]+=verdi;
         k++;
       }
     }  
     
/* garbage collection */

     free((double *)(s->ve));
     free((VECTOR *)(s)); 
     deltav++;
     free((double *)(deltav));   
     
/* find the error */

     r_norm=norm2(residual);
/* printf("\n residual norm: %16.14e\n",r_norm);           */
     
     err=r_norm/r_zero_norm;     
     
     counter++;
     if ((err < epsilon)||(r_norm==0.0)||(counter>15000)) {
       converged=1;
       vvec->ve[0]=10000.0;       
       buf_count=1;
       MPI_Bcast(&vvec->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD);
     }
     else {  /* do predictions */
       buf_count=residual->dim;                
       MPI_Bcast(&vvec->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD);
       
       for (t=1;t<=choice3;t++) {
         for (i=1;i<=ngroup;i++ ) {
           tag=i*100;
           MPI_Irecv(sndptr[i],buf_countmain,MPI_DOUBLE,MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,&req[i-1]);           
         }         
         deltav=CREATE(ngroup,double);
         --deltav;
         for (i=1;i<=ngroup;i++) {
           tildeA=tildeAptr[i-1];
           AA=sub_ptr[i-1]; 
           B=Rptr[i-1];    
           s=v_get(tildeA->n);
           alpha=-1.0;
           beta=0.0;
           trans='T';
           lda=tildeA->m;
           incx=1;
           incy=1;
           dgemv_(&trans,&tildeA->m,&tildeA->n,&alpha,tildeA->Value,&lda,vvec->ve,&incx,&beta,s->ve,&incy);
           uplo='L';
           nrhs=1;
           lda=B->m;
           ldb=B->n;
           dpotrs_(&uplo,&B->n,&nrhs,B->Value,&lda,s->ve,&ldb,&info);

/* take the littles out of s */

           littles=v_get(AA->n);
           k=AA->n;
           for (j=0;j<k;j++) {
             littles->ve[j]=s->ve[i+j-1];
           }                      
           sol=v_get(AA->m);
           trans='N';
           alpha=1.0;
           beta=0.0;
           lda=AA->m;
           incx=1;
           incy=1;
           dgemv_(&trans,&AA->m,&AA->n,&alpha,AA->Value,&lda,littles->ve,&incx,&beta,sol->ve,&incy);

/* send both s and sol to master */

           k=s->dim;
           for (j=0;j<k;j++)
             sndvector[j]=s->ve[j];
           count=sol->dim;
           for (j=0;j<count;j++) 
             sndvector[j+k]=sol->ve[j];      
           buf_count=k+count;
           tag=i*100;
           MPI_Isend(&sndvector[0],buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&req2[i-1]); 
           free((double *)(s->ve));
           free((double *)(sol->ve));
           free((VECTOR *)(s));
           free((VECTOR *)(sol)); 
           free((double *)(littles->ve));
           free((VECTOR *)(littles));
           MPI_Request_free(&req2[i-1]);       
         }
         for (i=1;i<=ngroup;i++ ) {
           MPI_Waitany(ngroup,req,&j,&status);
           dest=j+1;
           for (j=1;j<dest;j++)
             deltav[j]+=sndptr[dest][j-1];
           count=D[dest]->dim;
           for (j=0;j<count;j++) 
             D[dest]->ve[j]=sndptr[dest][j+dest-1];
           count+=dest;
           l=D[dest]->dim+ngroup-1;
           k=1;
           for (j=count;j<=l;j++) {
             deltav[k+dest]+=sndptr[dest][j-1];
             k++;
           }
           sind=D[dest]->dim+ngroup-1;
           eind=Ad[dest]->dim;
           for (j=0;j<eind;j++)
             Ad[dest]->ve[j]=sndptr[dest][j+sind];
         }
     
/* compute the columns of A*D=hatA matrix */
   
         for (i=1;i<=ngroup;i++) {
           eind=Ad[i]->dim;
           for (j=0;j<eind;j++) {
             Ad[i]->ve[j]+=Apptr[i-1]->ve[j]*deltav[i];
           }
         } 
     
/* put the nonzero values of hatA matrix */
     
         k=0;
         for (i=1;i<=ngroup;i++) {
           count=C->m;
           for (j=0;j<count;j++) {
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
         dsyrk_(&uplo,&trans,&j,&k,&alpha,C->Value,&lda,&beta,CC->Value,&ldc);
         lda=CC->n;
         j=CC->n;
         dpotf2_(&uplo,&j,CC->Value,&lda,&info);
         s=v_get(C->n);
         alpha=-1.0;
         beta=0.0;
         trans='T';
         lda=C->m;
         incx=1;
         incy=1;
         dgemv_(&trans,&C->m,&C->n,&alpha,C->Value,&lda,vvec->ve,&incx,&beta,s->ve,&incy);
         nrhs=1;
         lda=CC->m;
         ldb=CC->n;
         dpotrs_(&uplo,&CC->n,&nrhs,CC->Value,&lda,s->ve,&ldb,&info);

/* update the v vector and tmpApptr at the same time*/       

         for (i=0;i<ngroup;i++) {
           count=residual->dim;
           for (j=0;j<count;j++) {
             verdi=Ad[i+1]->ve[j]*s->ve[i];
             vvec->ve[j]+=verdi;
             tmpApptr[i]->ve[j]+=verdi;
           }
         }

/* update z vector */

         k=0;
         for (i=0;i<ngroup;i++) {
           count=D[i+1]->dim;
           for (j=0;j<count;j++) {
             verdi=s->ve[i]*(D[i+1]->ve[j]+pvecptr[i]->ve[j]*deltav[i+1]);
             zvec->ve[k]+=verdi;
             k++;
           }
         }                
         if (t!=choice3) {
           buf_count=residual->dim;
           MPI_Bcast(&vvec->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD);
         }
         else {
           k=0;
           for (i=0;i<ngroup;i++) {
             count=pvecptr[i]->dim;
             for (j=0;j<count;j++) {
               pvecptr[i]->ve[j]=zvec->ve[k]; 
               k++;
             }
           }
           for (i=0;i<ngroup;i++) { 
             count=residual->dim;
             for (j=0;j<count;j++)  
               Apptr[i]->ve[j]=tmpApptr[i]->ve[j]; 
             free((double *)(tmpApptr[i]->ve));
             free((VECTOR *)(tmpApptr[i]));
             free((double *)(Rptr[i]->Value));
             free((FMATRIX *)(Rptr[i]));
           } 
           buf_count=sfdouble*residual->dim*(ngroup+1);
           buffer=CREATE(buf_count,char);
           position=0;
           j=residual->dim; 
           MPI_Pack(&residual->ve[0],j,MPI_DOUBLE,buffer,buf_count,&position,MPI_COMM_WORLD);
           for (i=0;i<ngroup;i++) 
             MPI_Pack(&Apptr[i]->ve[0],j,MPI_DOUBLE,buffer,buf_count,&position,MPI_COMM_WORLD);
           MPI_Bcast(buffer,position,MPI_PACKED,0,MPI_COMM_WORLD);
           j=residual->dim;
           position=0;
           MPI_Unpack(buffer,buf_count,&position,&residual->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
           for (i=0;i<ngroup;i++) 
             MPI_Unpack(buffer,buf_count,&position,&Apptr[i]->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
           for (i=0;i<ngroup;i++) {
             tildeA=tildeAptr[i];
             AA=sub_ptr[i];
             putvaluesintildeA(tildeA,Apptr,AA,ngroup,i);
           } 
           free((double *)(zvec->ve));
           free((VECTOR *)(zvec));
           free((double *)(vvec->ve));
           free((VECTOR *)(vvec));
         } 
         free((double *)(s->ve));
         free((VECTOR *)(s));
         ++deltav;
         free((double *)(deltav));
       }       
     }     
   }
   endtime=MPI_Wtime();     
   difftime=endtime-starttime;  
   fprintf(fp,"\n Converged after %d loops.",counter);
   fprintf(fp,"\n Two norm of residual is: %16.14e ",r_norm);
   fprintf(fp,"\n Number of predictions: %d",choice3);
   fprintf(fp,"\n Elapsed time: %16.14e seconds \n",difftime); 
   fclose(fp);
   MPI_Finalize();
}    

