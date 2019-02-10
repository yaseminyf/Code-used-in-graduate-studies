/* jacobi_p0_dense_diffslave.c                                          */
/* written 14.09.1999                                             */
/* last revised 19.09.1999                                        */
/* this program is for checking the dense programming model       */
/* routines from LAPACK package are used for computing solutions  */
/* this version is with p=0 and using slave number 1 to ngroup */

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
double *sndvector,**sndptr;
int position;                                      
char *buffer; 
int kcount;

/**************************************************************************/
int main(argc,argv)
int argc;
char *argv[];
{  

/* local variables */
   
   static int i,j,k;
   int counter;
   int choice1,choice2,count;
   double verdi;
   double number,num,*ip;
   int nc,converged;
   double starttime, endtime, difftime;
   int numslave;
   int sfdouble,buf_countmain;
   
   MPI_Request *req;
   MPI_Status status;
   
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   sfdouble=sizeof(double);
   
   if ( myid == 0 ) {    /* master */
     sscanf(argv[1],"%d",&choice1); 
     sscanf(argv[2],"%d",&choice2); /* get the number of groups */
     sscanf(argv[3],"%d",&numslave); /* get the number of slaves */   
         
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
       case 17:
         if ((fp=fopen("bigdat","r"))==NULL) {
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
       case 17:
         fp=fopen("bigdatx","r");
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
     fp=fopen("res_p0_dense_diffslave","a");
     fprintf(fp,"\n file: %d , m: %d , n: %d, numslave: %d ",choice1,A->m,A->n,numslave);
     fprintf(fp,"\n initial norm: %16.14e ",r_zero_norm);

     kcount=ngroup/numslave;

/* broadcast how many groups each slave will have */

     MPI_Bcast(&kcount,1,MPI_INT,0,MPI_COMM_WORLD);
     
/* send each column group to a different processor */
     
     i=0;
     for (dest=1;dest<=numslave;dest++) {
       tag=dest;
       for ( k=1;k<=kcount;k++ ) {
         position=0;
         buf_count=2+sub_ptr[i]->n*sub_ptr[i]->m;
         buf_count=buf_count*sfdouble;
         buffer=CREATE(buf_count,char);
         MPI_Send(&buf_count,1,MPI_INT,dest,50,MPI_COMM_WORLD); 
         position=0;
         j=1;
         MPI_Pack(&sub_ptr[i]->m,j,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
         MPI_Pack(&sub_ptr[i]->n,j,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
         j=sub_ptr[i]->m*sub_ptr[i]->n;
         MPI_Pack(&sub_ptr[i]->Value[0],j,MPI_DOUBLE,buffer,buf_count,&position,MPI_COMM_WORLD);
         MPI_Send(buffer,position,MPI_PACKED,dest,tag,MPI_COMM_WORLD); 
         free((char *)(buffer));
         i++;
       }
     }       
          
/* broadcast the residual and the number of columns in each block */     

     buf_count=residual->dim+1;
     buf_count=buf_count*sfdouble; 
     buffer=CREATE(buf_count,char);
     position=0;
     MPI_Pack(&nc,1,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
     j=residual->dim;
     MPI_Pack(&residual->ve[0],j,MPI_DOUBLE,buffer,buf_count,&position,MPI_COMM_WORLD);
     MPI_Bcast(buffer,position,MPI_PACKED,0,MPI_COMM_WORLD);
     free((char *)(buffer));     

     D=CREATE(ngroup,VECTOR *);
     Ad=CREATE(ngroup,VECTOR *);
     --D;
     --Ad;
     for (i=0;i<ngroup;i++ ) {
       D[i+1]=v_get(sub_ptr[i]->n);
       Ad[i+1]=v_get(residual->dim);
     }
       
/* garbage collection */
     
     for (i=0;i<ngroup;i++) {
       free((double *)(sub_ptr[i]->Value));
       free((FMATRIX *)(sub_ptr[i])); 
     }
     free((FMATRIX **)(sub_ptr));

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
     buf_countmain=nc+residual->dim;
     sndptr=CREATE(ngroup,double *);
     --sndptr;
     for (i=1;i<=ngroup;i++) {
       sndvector=CREATE(buf_countmain,double);
       sndptr[i]=sndvector;
     }
     req=CREATE(ngroup,MPI_Request);
     MPI_Barrier(MPI_COMM_WORLD);   
     starttime=MPI_Wtime();  
     while (converged==0) {
       i=1;
       for (dest=1;dest<=numslave;dest++) {
         for (k=1;k<=kcount;k++ ) {
           MPI_Irecv(sndptr[i],buf_countmain,MPI_DOUBLE,dest,500,MPI_COMM_WORLD,&req[i-1]);           
           i++;
         }
       }        
       for (i=1;i<=ngroup;i++ ) {
         MPI_Waitany(ngroup,req,&j,&status);
         dest=j+1;
         count=D[dest]->dim;
         for (j=0;j<count;j++) 
           D[dest]->ve[j]=sndptr[dest][j];
         for (j=0;j<Ad[dest]->dim;j++)
           Ad[dest]->ve[j]=sndptr[dest][j+count];
       }       

/* put the nonzero values of hatA matrix */
     
       k=0;
       for (i=1;i<=ngroup;i++) {
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
     
/* update the residual vector */       

       for (i=1;i<=ngroup;i++) {
         count=residual->dim;
         for (j=0;j<count;j++)
           residual->ve[j]+=Ad[i]->ve[j]*s->ve[i-1];
       }

/* find the error */

       r_norm=norm2(residual);
/* printf("\n residual norm: %16.14e\n",r_norm);           */
     
       err=r_norm/r_zero_norm;     
     
       counter++;
       if ((err < epsilon)||(r_norm==0.0)||(counter>15000)) {
         converged=1;
         residual->ve[0]=10000.0;
         buf_count=1; 
         MPI_Bcast(&residual->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD); 
       }
       else {
         buf_count=residual->dim;                
         MPI_Bcast(&residual->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD);
         if (counter==1) endtime=MPI_Wtime();     

       } 

/* update x vector */

       k=0;
       for (i=1;i<=ngroup;i++) {
         count=D[i]->dim;
         for (j=0;j<count;j++) {
           x->ve[k]+=D[i]->ve[j]*s->ve[i-1];
           k++;
         }
       }

/* garbage collection */

       free((double *)(s->ve));
       free((VECTOR *)(s));                  
     }
/*     endtime=MPI_Wtime();      */
     difftime=endtime-starttime;  
     fprintf(fp,"\n Converged after %d loops.",counter);
     fprintf(fp,"\n Two norm of residual is: %16.14e",r_norm);
     fprintf(fp,"\n Elapsed time: %16.14e seconds \n",difftime);
     fclose(fp); 
   }
   else { /* slaves */

/* get the number of groups to handle */
   
     MPI_Bcast(&kcount,1,MPI_INT,0,MPI_COMM_WORLD);
     
     sub_ptr=CREATE(kcount,FMATRIX *);
     for (i=0;i<kcount;i++) 
       sub_ptr[i]=CREATE(1,FMATRIX);
     tag=myid;

     for (i=0;i<kcount;i++) {
       MPI_Recv(&buf_count,1,MPI_INT,0,50,MPI_COMM_WORLD,&status);
       buffer=CREATE(buf_count,char);
       position=0;
       MPI_Recv(buffer,buf_count,MPI_PACKED,0,tag,MPI_COMM_WORLD,&status);
       MPI_Unpack(buffer,buf_count,&position,&sub_ptr[i]->m,1,MPI_INT,MPI_COMM_WORLD);
       MPI_Unpack(buffer,buf_count,&position,&sub_ptr[i]->n,1,MPI_INT,MPI_COMM_WORLD);
       j=sub_ptr[i]->m*sub_ptr[i]->n;
       sub_ptr[i]->Value=CREATE(j,double);
       MPI_Unpack(buffer,buf_count,&position,&sub_ptr[i]->Value[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
       free((char *)(buffer));
     }
     
     buf_count=1+sub_ptr[0]->m;
     buf_count=buf_count*sfdouble; 
     buffer=CREATE(buf_count,char);
     position=0;
     MPI_Bcast(buffer,buf_count,MPI_PACKED,0,MPI_COMM_WORLD);
     MPI_Unpack(buffer,buf_count,&position,&nc,1,MPI_INT,MPI_COMM_WORLD);
     residual=v_get(sub_ptr[0]->m);
     j=residual->dim;
     MPI_Unpack(buffer,buf_count,&position,&residual->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
     free((char *)(buffer));

/* do the factorization of each block */
/* first compute A_i^T A_i */
/* then do the cholesky factorization for each of them */
/* and save the resulting R matrices */   

     Rptr=CREATE(kcount,FMATRIX *);
     for (i=0;i<kcount;i++) {
   
/* call symmetric multiplication BLAS routine */
   
       alpha=1.0;
       beta=0.0;
       uplo='L';
       trans='T';
       j=sub_ptr[i]->n;
       B=CREATE(1,FMATRIX);
       B->m=j;
       B->n=j;
       numnz=B->m*B->n;
       B->Value=CREATE(numnz,double); 
       lda=sub_ptr[i]->m;
       k=sub_ptr[i]->m;
       ldc=sub_ptr[i]->n;
       dsyrk_(&uplo,&trans,&j,&k,&alpha,sub_ptr[i]->Value,&lda,&beta,B->Value,&ldc);

/* do the Cholesky factorization */

       uplo='L';
       lda=B->m;
       k=B->n;
       dpotf2_(&uplo,&k,B->Value,&lda,&info);
       Rptr[i]=B; 
     }

/* begin the solution */

     tag=500; 
     
     converged=0;
     buf_count=nc+residual->dim;
     sndvector=CREATE(buf_count,double);
     MPI_Barrier(MPI_COMM_WORLD); 

     while ( converged != 1 ) {
       for (i=0;i<kcount;i++) {   
         AA=sub_ptr[i];
         s=v_get(AA->n);
         alpha=-1.0;
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
         dpotrs_(&uplo,&B->n,&nrhs,B->Value,&lda,s->ve,&ldb,&info);
         sol=v_get(AA->m);
         trans='N';
         alpha=1.0;
         beta=0.0;
         lda=AA->m;
         incx=1;
         incy=1;
         dgemv_(&trans,&AA->m,&AA->n,&alpha,AA->Value,&lda,s->ve,&incx,&beta,sol->ve,&incy);
         count=sub_ptr[i]->n;
         for (j=0;j<count;j++)
           sndvector[j]=s->ve[j];
         count=sol->dim;
         k=s->dim;
         for (j=0;j<count;j++) 
           sndvector[j+k]=sol->ve[j]; 
         buf_count=k+count;
         MPI_Send(&sndvector[0],buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD); 
         free((double *)(s->ve));
         free((double *)(sol->ve));
         free((VECTOR *)(s));
         free((VECTOR *)(sol));  
       }

/* probe for a final signal or new residual from master */

       buf_count=residual->dim;
       MPI_Bcast(&residual->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD);  
    
       if (residual->ve[0]==10000.0) {
         converged=1;
       }
     }
   }   /* end of slave */       
         
   MPI_Finalize();
}

