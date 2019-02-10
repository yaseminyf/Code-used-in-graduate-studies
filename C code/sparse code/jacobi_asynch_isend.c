/* last revised 27.01.1999                                                  */
/* jacobi_asynch_isend = routine for solving the LSQ problem asynchronously */
/* has options linesearch, pure Jacobi or planar search after j=g,2g,3g     */
/* updates                                                                  */
/****************************************************************************/

#include "mpi.h"

#include <unistd.h>
#include <limits.h>
#include "data_struct.h"


#define PONTUS 

/* #define LINESEARCH    */

#ifdef LINESEARCH
  double check, dotnorm;
#endif

#define SUBSPACECOR      

  char *sndvector,**sndptr;
  int position;
  VECTOR **D,**DS;
  MATRIX *C;
  
/* #define PERM   */
/* #define MMDP  */ 
/* #define COLP     */
/* #define NEWPERM   */

#ifdef PERM
  int *p,*ivp;
  double *tmp;
#endif

/* #define XERROR */  /* to check convergence with error on x */ 

#define DIFFGROUPSIZE 

#ifdef XERROR
  VECTOR *xnormv;
#endif
  
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
   double *ip,number,num;
   int nc;
   double starttime, endtime, difftime;
   
#ifdef DIFFGROUPSIZE   
   int choice2;
#endif   

   MPI_Request request,*req;
   MPI_Status status;
   
/* initialize MPI and identify yourself */

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

/* initialize for SPARSPAK */ 

   spbusr_.MAXSB=6000;
   TOL=0.0000000001;
   TYPTOL=1;

   spksys_.MAXINT=SHRT_MAX;
    
   spbusr_.IERRB=0;
   spbusr_.MSGLVB=1;     

   if ( myid == 0 ) {    /* master */
     sscanf(argv[1],"%c",&choice1); 

#ifdef DIFFGROUPSIZE   
     sscanf(argv[2],"%d",&choice2); /* get the number of groups */
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
     /*  free((double *)(xexact->ve));
       free((VECTOR *)(xexact));    */  
#endif   

     }
     else if ( (choice1 >= '6') && (choice1<= '9') ) {
        epsilon=1.0e-03;
        fp1=fopen("data","w");
        readnewmatrix(fp,fp1);
     }
     x=v_get(A->n);

     residual=v_get(rhs->dim);  
     residual=m_fulvectormlt(A,x,residual);
     for ( i=0;i<rhs->dim;i++ ) {
       residual->ve[i]=residual->ve[i]-rhs->ve[i];
     }

/* garbage collection */

 /*    free((double *)(rhs->ve));
     free((VECTOR *)(rhs)); */
     
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

  /*     free((int *)(A->Jind));
       free((int *)(A->Iind));
       free((int *)(A->colind));
       free((int *)(A->rowptr));
       free((double *)(A->Valuerow));
       free((double *)(A->Valuecol));
       free((MATRIX *)(A));        */
       
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
        printf("\n initial norm: %16.14e",r_zero_norm);
     }

/* send each column group to a different processor */

     for ( i=0;i<ngroup;i++ ) {
       AA=sub_ptr[i];
       dest=i+1;
       
       buf_count=1;
       tag=1;
       MPI_Send(&AA->m,buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);
      
       tag=2;
       MPI_Send(&AA->n,buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);
      
       tag=3;
       MPI_Send(&AA->nonzeros,buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);
 
       tag=4;
       buf_count=AA->nonzeros+1;
       MPI_Send(&AA->Iind[0],buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);

       tag=5;
       buf_count=AA->n+2;
       MPI_Send(&AA->Jind[0],buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);

       tag=6;
       buf_count=AA->m+2;
       MPI_Send(&AA->rowptr[0],buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);
      
       tag=7;
       buf_count=AA->nonzeros+1;
       MPI_Send(&AA->colind[0],buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);

       tag=8;
       MPI_Send(&AA->Valuecol[0],buf_count,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);

       tag=9;
       MPI_Send(&AA->Valuerow[0],buf_count,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
       
       tag=11;
       buf_count=AA->m+1;
       zeroptr_par=CREATE(AA->m+1,int);
       for ( j=0;j<=AA->m;j++ ) 
         zeroptr_par[j]=zeroptr[i][j];
       MPI_Send(&zeroptr_par[0],buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);
       free((int *)(zeroptr_par));
     }

/* broadcast the residual */     

     buf_count=residual->dim;
     MPI_Bcast(&residual->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD);

#ifdef SUBSPACECOR

/* do the symmetric factorization of the C matrix */

     C=CREATE(1,MATRIX);
     C->m=residual->dim;     
     C->n=ngroup;
     createC1(sub_ptr,C,ngroup);

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
     spbcon_.NZCONS=0;spbusr_.MCOLS=C->n;
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
#endif

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

     omega=1.0;
     counter=0;
     converged=0;

     D=CREATE(ngroup,VECTOR *);
     --D;
     DS=CREATE(ngroup,VECTOR *);
     --DS;
     sndptr=CREATE(ngroup,char *);
     for (i=0;i<ngroup;i++) { 
       sndvector=CREATE(100000,char);
       sndptr[i]=sndvector;
     }  

     req=CREATE(ngroup+1,MPI_Request);
printf("\n before loop...\n");
     fp1=fopen("res","w");
     fp=fopen("check","w");
     MPI_Barrier(MPI_COMM_WORLD);       
     starttime=MPI_Wtime();  

     while (converged==0) { 
       buf_count=100000;
       for (i=1;i<=ngroup;i++) {
         position=0;
         MPI_Irecv(&sndptr[i-1][0],buf_count,MPI_PACKED,i,50,MPI_COMM_WORLD,&req[i]);
       }     
       for (i=1;i<=ngroup;i++) {
         MPI_Waitany(ngroup+1,req,&j,&status);
         dest=status.MPI_SOURCE;
         buf_count=100000;
         position=0;
         MPI_Unpack(&sndptr[dest-1][0],buf_count,&position,&verdi,1,MPI_DOUBLE,MPI_COMM_WORLD);
         sol=v_get(residual->dim);
         MPI_Unpack(&sndptr[dest-1][0],buf_count,&position,&sol->ve[0],residual->dim,MPI_DOUBLE,MPI_COMM_WORLD);
         s=v_get((int)(verdi));
         MPI_Unpack(&sndptr[dest-1][0],buf_count,&position,&s->ve[0],s->dim,MPI_DOUBLE,MPI_COMM_WORLD);
         D[dest]=s;        
         DS[dest]=sol;   
#ifdef LINESEARCH
         check=0.0;
         for ( j=0;j<residual->dim;j++ )
           check=check+sol->ve[j]*residual->ve[j];
         dotnorm=0.0;
         for ( j=0;j<residual->dim;j++)
           dotnorm+=sol->ve[j]*sol->ve[j];
         omega=-check/dotnorm;
#endif 
         k=(dest-1)*nc;
         for (j=0;j<D[dest]->dim;j++) {
           x->ve[k+j]=x->ve[k+j]+D[dest]->ve[j]*omega;
         }
      
         for ( j=0;j<residual->dim;j++ ) {
           if (DS[dest]->ve[j]!=0.0) { 
             residual->ve[j]+=DS[dest]->ve[j]*omega;
           } 
         }
         
#ifdef SUBSPACECOR
         j=2*ngroup;
         if ( ((counter+1) % j) == 0) { /* time to do planar search */

/* compute C^T*C and initialize lnz and diag */

           ddiag=CREATE(spbcon_.NCOLS,double);
           lnzv=CREATE(spbcon_.NOFNZ,double);
           --ddiag;
           --lnzv; 

/* get the necessary storage for C^TC */

           j=spbcon_.NOFNZ+spbcon_.NCOLS;
           isub=CREATE(j,int);
           jsub=CREATE(j,int);
           values=CREATE(j,double);
           --isub;
           --jsub;
           --values;
           createC22(C,DS,ngroup);
           sym_mult(C,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values);

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
          for (l=1;l<=ngroup;l++) {
            for (j=1;j<=D[l]->dim;j++) {
              x->ve[k]+=D[l]->ve[j-1]*s->ve[l-1];
              k++;
            }
          }
     
/* update the residual vector */

          if ( choice1 >= '6' ) {
            sol=v_get(residual->dim);
          }
          for (l=1;l<=ngroup;l++) {
            for ( j=0;j<residual->dim;j++ ) {
              if (DS[l]->ve[j]!=0.0) { 
                residual->ve[j]+=DS[l]->ve[j]*s->ve[l-1];
                if ( choice1 >= '6' ) { 
                  sol->ve[j]+=DS[l]->ve[j]*s->ve[l-1];
                }
              }  
            }  
          }
            
/* garbage collection */

          free((double *)(s->ve));
          free((VECTOR *)(s));         
          ++ddiag;
          ++lnzv;
          free((double *)(C->Valuerow));
          free((double *)(C->Valuecol));
          free((double *)(ddiag));
          free((double *)(lnzv)); 
          for (l=1;l<=ngroup;l++) {
            free((double *)(D[l]->ve)); 
            free((double *)(DS[l]->ve));
          }   
        }                             
#endif
           
/* find the error */

        if ( choice1 <= '5' ) {
          r_norm=norm2(residual);

printf("\n residual norm: %16.14e\n",r_norm); 
 fprintf(fp1,"\n  %16.14e",r_norm);     
       
        }
        else {
          rvec2=multAtransr2(A,sol,rvec2);
          r_norm=norm2(rvec2);
printf("\n residual norm: %16.14e\n",r_norm); 
          free((double *)(sol->ve));
          free((VECTOR *)(sol));
        }
        err=r_norm/r_zero_norm;            

/* check ..... */
        sol=v_get(residual->dim);
        sol=m_fulvectormlt(A,x,sol);
        for (j=0;j<sol->dim;j++) {
          sol->ve[j]=sol->ve[j]-rhs->ve[j];
        }
        verdi=norm2(sol);
printf("\n checked residual norm %16.14e\n",verdi);        
fprintf(fp,"\n %16.14e",verdi); 
free((double *)(sol->ve));
free((VECTOR *)(sol));
        counter++;
        if ((err < epsilon)||(r_norm==0.0)) { /* send final message */
          converged=1;     
          for (j=1;j<=ngroup;j++) {
            if (req[j]!=MPI_REQUEST_NULL) {
              MPI_Cancel(&req[j]);
              MPI_Request_free(&req[j]);
            }
          }
          for (j=1;j<=ngroup;j++) 
            MPI_Isend(&tag,1,MPI_INT,j,1000,MPI_COMM_WORLD,&req[j]);          
          break;
        }
        else {
          buf_count=residual->dim;
          MPI_Isend(&residual->ve[0],buf_count,MPI_DOUBLE,dest,100,MPI_COMM_WORLD,&request);
        } 
      }                
    }
    endtime=MPI_Wtime();
    difftime=endtime-starttime;     
    printf("\n Converged after %d loops.",counter);
    printf("\n Two norm of residual is: %16.14e \n",r_norm);  
    printf("\n Elapsed time: %lf seconds \n",difftime);   
    sol=v_get(x->dim);
    for (j=0;j<sol->dim;j++) 
          sol->ve[j]=x->ve[j]-xexact->ve[j];
verdi=norm2(sol);
printf("\n norm in the difference of xes %16.14e\n",verdi);      
          
  }    
  else  {  /* slave */
     A=CREATE(1,MATRIX);
     
     A->rowdefined=TRUE;
     zeroptr_par=CREATE(A->m+1,int);
     
     tag=1;
     buf_count=1;
     MPI_Recv(&A->m,buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag=2;
     MPI_Recv(&A->n,buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag=3;
     MPI_Recv(&A->nonzeros,buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag=4;
     buf_count=A->nonzeros+1;
     A->Iind=CREATE(buf_count,int);
     MPI_Recv(&A->Iind[0],buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag=5;
     buf_count=A->n+2;
     A->Jind=CREATE(buf_count,int);
     MPI_Recv(&A->Jind[0],buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag=6;
     buf_count=A->m+2;
     A->rowptr=CREATE(buf_count,int);
     MPI_Recv(&A->rowptr[0],buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag=7;
     buf_count=A->nonzeros+1;
     A->colind=CREATE(buf_count,int);
     MPI_Recv(&A->colind[0],buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     tag=8;
     A->Valuecol=CREATE(buf_count,double);
     MPI_Recv(&A->Valuecol[0],buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);

     tag=9;
     A->Valuerow=CREATE(buf_count,double);
     MPI_Recv(&A->Valuerow[0],buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);

     tag=11;
     buf_count=A->m+1;
     zeroptr_par=CREATE(buf_count,int);
     MPI_Recv(&zeroptr_par[0],buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

     buf_count=A->m;
     residual=v_get(buf_count);
     MPI_Bcast(&residual->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD);
     
/* initialize for SPARSPAK */
   
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
    
     spbusr_.MCOLS=A->n;
     spbusr_.MSEQNS=A->m-zeroptr_par[0];
    
     spbcon_.NCOLS=spbusr_.MCOLS;  
        
     T=CREATE(spbusr_.MAXSB,int);
     for (i=0;i<spbusr_.MAXSB;i++) {
       T[i]=-(i+1);
     }
     --T;

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
         free((int *)(SUBS));
       }
     }
    
     orcolb(T);
 
/* compute A^T*A and initialize lnz and diag */

     ddiag=CREATE(spbcon_.NCOLS,double);
     lnzv=CREATE(spbcon_.NOFNZ,double);
     --ddiag;
     --lnzv;

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
     
/* get the necessary storage for A^TA */

     i=spbcon_.NOFNZ+spbcon_.NCOLS;
     isub=CREATE(i,int);
     jsub=CREATE(i,int);
     values=CREATE(i,double);
     --isub;
     --jsub;
     --values;
     
     sym_mult(A,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values);   

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

/* do Cholesky factorization */ 

     gsfct(spbcon_.NCOLS,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,llink,first,temp,flag__); 

/* garbage collection */

     ++llink;
     ++temp;
     ++first;
     free((int *)(llink));
     free((double *)(temp));
     free((int *)(first));
   
/* begin the solution */

     s=v_get(A->n);
     sol=v_get(A->m);

     verdi=(double)(A->n);  
     tag=50;
     sndvector=CREATE(100000,char);

     converged=0;

     MPI_Barrier(MPI_COMM_WORLD);  
     while ( converged != 1 ) {
       s=multAtransr3(A,residual,s);
       s=permdv(A->n,s,invpv);
       s=gsslv(A->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,s);

/* permute the solution to its original order of columns */

       s=permdv(A->n,s,permv);    
       sol=m_fulvectormlt(A,s,sol);
  
/* send both s and sol to master */

       buf_count=100000;
       position=0; 
       MPI_Pack(&verdi,1,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
       j=sol->dim;
       MPI_Pack(&sol->ve[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
       j=s->dim;
       MPI_Pack(&s->ve[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
       MPI_Isend(sndvector,position,MPI_PACKED,0,tag,MPI_COMM_WORLD,&request); 
       free((double *)(s->ve));
       free((double *)(sol->ve));
       free((VECTOR *)(s));
       free((VECTOR *)(sol));         
       MPI_Wait(&request,&status);
          
/* probe for a final signal or new residual from master */

       MPI_Probe(0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
               
       if (status.MPI_TAG==100) {
         MPI_Recv(&residual->ve[0],buf_count,MPI_DOUBLE,0,100,MPI_COMM_WORLD,&status);     
         s=v_get(A->n);
         sol=v_get(A->m);
       }
       else {    /* final message; gather the x vector */
         MPI_Recv(&tag,1,MPI_INT,0,1000,MPI_COMM_WORLD,&status);           
         converged=1;
       }
     }
   }   /* end of slave */  
   MPI_Finalize(); 
}    
       
