/* last revised 24.09.1998                                         */
/* asynch_GS.c = routine for solving the LSQ problem asynchronously*/
/* and using GS updates                                            */
/* it is possible to choose either the QR or the Cholesky approach */
/* by changing the 'define' switch in 'data_struct.h'              */
/*******************************************************************/
#include "mpi.h"

#include <unistd.h>
#include <limits.h>
#include "data_struct.h"

/* #define DEBUG   */

#define CHOL 

/* #define QR  */

/* #define PERM  */
/* #define MMDP  */ 
/* #define COLP     */
/* #define NEWPERM  */

#ifdef PERM
  int *p,*ivp;
  double *tmp;
#endif
 
#define DIFFGROUPSIZE 

/* variables used for parallel operations */

int myid, numprocs;                         /* to identify the processors */
int dest, tag;                   /* node to send/receive msg, and msg tag */
int buf_count;                       /* the number of elements to be sent */
                                      
/**************************************************************************/
int main(argc,argv)
int argc;
char *argv[];
{  

/* local variables */
   
   static int i,j,k,l;
   int rnum, count, counter;
   char choice1,c1;
   double *largex,verdi;
   MATRIX *AA;
   int converged;
   int sizex;
   double number,num,*ip;
   int nc;
   double starttime, endtime, difftime;

#ifdef DIFFGROUPSIZE   
   int choice2;
#endif   

   MPI_Request req;
   MPI_Status status;
   
/* initialize MPI and identify yourself */

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
     
   if ( myid == 0 ) {    /* master */
     sscanf(argv[1],"%c",&choice1); 
/*     choice1='2';  */

#ifdef DIFFGROUPSIZE   
   choice2=6;
#endif                

     switch(choice1) {
        case '1':
          if ((fp=fopen("matrix1","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          if ((fp1=fopen("matrix1dat","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          break;
        case '2':
          if ((fp=fopen("matrix2","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          if ((fp1=fopen("matrix2dat","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          break;
        case '3':
          if ((fp=fopen("matrix3","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          if ((fp1=fopen("matrix3dat","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          break;
        case '4':
          if ((fp=fopen("matrix4","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          if ((fp1=fopen("matrix4dat","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          break;
        case '5':
          if ((fp=fopen("matrix5","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          if ((fp1=fopen("matrix5dat","r"))==NULL) {
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
        epsilon=1.0e-5;
        readmatrix(fp,fp1);
        xexact=v_get(A->n);
        
#ifdef DEBUG    
        fprintf(fp1,"\n got the xexact vector..");
#endif   
printf("\n got the xexact vector..");
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
          fscanf(fp,"%le\n",&verdi);
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
        epsilon=1.0e-03;
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

     for ( i=0;i<rhs->dim;i++ ) {
       residual->ve[i]=-1.0*rhs->ve[i];
     }
     
#ifdef DEBUG    
     fprintf(fp1,"\n computed the residual..");
/* record the initial data in a file */

     fprintf(fp1,"\n\nThe initial residual:\n\n");
     for ( i=0;i<residual->dim;i++ ) {
       fprintf(fp1,"%16.15e ",residual->ve[i]);
       if ( ((i+1)%5)==0 )
         fprintf(fp1,"\n");
     } 
#endif   

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

/* send each column group to a different processor */
printf("\n before sending...");
     for ( i=0;i<ngroup;i++ ) {
       AA=sub_ptr[i];
       dest=i+1;
       
       buf_count=1;
       tag=1;
       MPI_Send(&AA->m,buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);
      
       tag++;
       MPI_Send(&AA->n,buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);
      
       tag++;
       MPI_Send(&AA->nonzeros,buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);
 
       tag++;
       buf_count=AA->nonzeros+1;
       MPI_Send(&AA->Iind[0],buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);

       tag++;
       buf_count=AA->n+2;
       MPI_Send(&AA->Jind[0],buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);

       tag++;
       buf_count=AA->m+2;
       MPI_Send(&AA->rowptr[0],buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);
      
       tag++;
       buf_count=AA->nonzeros+1;
       MPI_Send(&AA->colind[0],buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);

       tag++;
       MPI_Send(&AA->Valuecol[0],buf_count,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);

       tag++;
       MPI_Send(&AA->Valuerow[0],buf_count,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
      
       tag++;
       buf_count=AA->m+1;
       zeroptr_par=CREATE(AA->m+1,int);
       for ( j=0;j<=AA->m;j++ ) 
         zeroptr_par[j]=zeroptr[i][j];
       MPI_Send(&zeroptr_par[0],buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);
       free(zeroptr_par);

       tag++;
       buf_count=1;
       sizex=A->n;
       MPI_Send(&sizex,buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);
     }
printf("\n finished sending..begin broadcasts..");
/* broadcast the residual */     

     tag++;
     buf_count=residual->dim;
     MPI_Bcast(&residual->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD);

/* broadcast the matrix number */

     tag++;
     buf_count=1;
     MPI_Bcast(&choice1,buf_count,MPI_CHAR,0,MPI_COMM_WORLD);

/* broadcast the number of columns in each block */
     
     tag++;
     buf_count=1;
     MPI_Bcast(&nc,buf_count,MPI_CHAR,0,MPI_COMM_WORLD);
     
printf("\n finished broadcasts..");     
/* begin the solution */
     
     r_norm=r_zero_norm;
     err=r_norm/r_zero_norm;
     omega=1.0;
     counter=0;
     sol=v_get(residual->dim);
     buf_count=sol->dim;
     converged=0;
printf("\n just before while loop.."); 
     MPI_Barrier(MPI_COMM_WORLD);  
     starttime=MPI_Wtime();  
     while ((converged==0)&&(counter<1000)) {
       MPI_Recv(&sol->ve[0],buf_count,MPI_DOUBLE,MPI_ANY_SOURCE,50,MPI_COMM_WORLD,&status);
printf("\n received from %d, r_norm %16.14e",status.MPI_SOURCE,r_norm);
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
       counter++;
       if ((err < epsilon)||(r_norm==0.0)||(counter>=1000)) { /* send final message */
         converged=1;     
         for ( dest=1;dest<=ngroup;dest++)
           MPI_Send(&tag,1,MPI_INT,dest,1000,MPI_COMM_WORLD);
         largex=CREATE(sizex,double);
         buf_count=sizex;
         MPI_Reduce(&x->ve[0],&largex[0],buf_count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
         endtime=MPI_Wtime();
       }
       else {
         dest=status.MPI_SOURCE;
         MPI_Isend(&residual->ve[0],buf_count,MPI_DOUBLE,dest,500,MPI_COMM_WORLD,&req);
       }
     }
     difftime=endtime-starttime;
     printf("\n Converged after %d updates.",counter);
     printf("\n Two norm of residual is: %16.14e \n",r_norm); 
     printf("\n Elapsed time: %f \n", difftime);       
   }    
   else  {  /* slave */
     A=CREATE(1,MATRIX);
     
     A->rowdefined=TRUE;
     zeroptr_par=CREATE(A->m+1,int);
     
     tag=1;
     buf_count=1;
     MPI_Recv(&A->m,buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag++;
     MPI_Recv(&A->n,buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag++;
     MPI_Recv(&A->nonzeros,buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag++;
     buf_count=A->nonzeros+1;
     A->Iind=CREATE(buf_count,int);
     MPI_Recv(&A->Iind[0],buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag++;
     buf_count=A->n+2;
     A->Jind=CREATE(buf_count,int);
     MPI_Recv(&A->Jind[0],buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag++;
     buf_count=A->m+2;
     A->rowptr=CREATE(buf_count,int);
     MPI_Recv(&A->rowptr[0],buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag++;
     buf_count=A->nonzeros+1;
     A->colind=CREATE(buf_count,int);
     MPI_Recv(&A->colind[0],buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag++;
     A->Valuecol=CREATE(buf_count,double);
     MPI_Recv(&A->Valuecol[0],buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);

     tag++;
     A->Valuerow=CREATE(buf_count,double);
     MPI_Recv(&A->Valuerow[0],buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);

     tag++;
     buf_count=A->m+1;
     zeroptr_par=CREATE(buf_count,int);
     MPI_Recv(&zeroptr_par[0],buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag++;
     buf_count=1;
     MPI_Recv(&sizex,buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag++;
     buf_count=A->m;
     residual=v_get(buf_count);
     MPI_Bcast(&residual->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD);

     tag++;
     buf_count=1;
     MPI_Bcast(&choice1,buf_count,MPI_CHAR,0,MPI_COMM_WORLD);

     tag++;
     buf_count=1;
     MPI_Bcast(&nc,buf_count,MPI_INT,0,MPI_COMM_WORLD);
     
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
printf("\n sparspak initialization ends..");   
     spbusr_.MCOLS=A->n;
     spbusr_.MSEQNS=A->m-zeroptr_par[0];
    
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
    
     orcolb(T);
   
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
#endif

#ifdef CHOL
 
/* compute A^T*A and initialize lnz and diag */

     ddiag=CREATE(spbcon_.NCOLS,double);
     lnzv=CREATE(spbcon_.NOFNZ,double);
     --ddiag;
     --lnzv;
     for (j=1;j<=spbcon_.NOFNZ;j++)
       lnzv[j]=0.0;
     for (j=1;j<=spbcon_.NCOLS;j++)
       ddiag[j]=0.0;    
/*     zerols(spbcon_.NOFNZ, T, xrnz, lnzv); */
/*     zerorv(spbcon_.NCOLS, ddiag); */

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
   
/* do Cholesky factorization */ 

     gsfct(spbcon_.NCOLS,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,llink,first,temp,flag__); 

#endif
   
/* begin the solution */

     s=v_get(A->n);
     sol=v_get(A->m);
     x=v_get(sizex);            /* initial x vector to be updated */
     omega=1.0;
     buf_count=sol->dim;
     tag=50;
     srand(32755);
     converged=0;
     MPI_Barrier(MPI_COMM_WORLD);  
     while ( converged != 1 ) {
       s=multAtransr3(A,residual,s);
#ifdef QR   
       s=solution(A->n,s,rnzv,ddiag,xnzsubv,nzsubv,xrnzv,permv,invpv);
#endif

#ifdef CHOL

       s=permdv(A->n,s,invpv);
       s=gsslv(A->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,s);

/* permute the solution to its original order of columns */

       s=permdv(A->n,s,permv);   

#endif   

       sol=m_fulvectormlt(A,s,sol);

/* send solution to master */

       MPI_Isend(&sol->ve[0],buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&req);
      
/*       for ( j=0;j<s->dim;j++ ) 
         x->ve[(myid-1)*nc+j]=x->ve[(myid-1)*nc+j]+omega*s->ve[j]; */
   
       j=rand();
       j=j%20;
       sleep(j); 
/*       MPI_Wait(&req,&status);  */
       
/* probe for a final signal or new residual from master */

       MPI_Probe(0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
       if (status.MPI_TAG == 500 ) {      /* receive residual */
         MPI_Recv(&residual->ve[0],buf_count,MPI_DOUBLE,0,500,MPI_COMM_WORLD,&status);     
         s=v_zero(s);
         sol=v_zero(sol);
       }
       else {    /* final message; gather the x vector */
          MPI_Recv(&tag,1,MPI_INT,0,1000,MPI_COMM_WORLD,&status);
          MPI_Reduce(&x->ve[0],&largex[0],buf_count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); 
          
/* garbage collection */

          ++xrnzv;
          ++ddiag; 
          ++invpv;
          ++permv;
          ++nzsubv;
          ++xnzsubv;
          free(invpv);
          free(permv);
          free(xrnzv);
          free(ddiag);
          free(nzsubv);
          free(xnzsubv);

#ifdef CHOL
          ++lnzv;
          free(lnzv);
#endif

#ifdef QR
          ++rnzv;
          free(rnzv);   
#endif    
          converged=1;
     
       }
     }
   }   /* end of slave */  
   MPI_Finalize(); 
}    
       
