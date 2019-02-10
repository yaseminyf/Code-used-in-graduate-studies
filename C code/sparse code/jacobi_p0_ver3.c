/* last revised 16.04.1999                                         */
/* jacobi.c = routine for solving the LSQ problem synchronously    */
/* and using Jacobi updates and planar search                      */
/* uses All-reduce command on g+1 processors                       */
/*******************************************************************/
#include "mpi.h"

#include <unistd.h>
#include <limits.h>
#include "data_struct.h"

#define CHOL  

#define PONTUS 

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

double *sndvector;
int position;
MATRIX *C;
VECTOR **D;
double *largex;
  
/* variables used for parallel operations */

int myid, numprocs,sid;                     /* to identify the processors */
int dest, tag;                   /* node to send/receive msg, and msg tag */
int buf_count;                       /* the number of elements to be sent */
VECTOR **DS;
double ss;
int flag;
                                      
/**************************************************************************/
int main(argc,argv)
int argc;
char *argv[];
{  

/* local variables */
   
   static int i,j,k,l;
   int rnum, count, counter;
   char c1;
   double *largex,verdi;
   MATRIX *AA;
   int converged;
   int sizex;
   double *ip,number,num;
   int nc,sind,eind,choice1;
   double starttime, endtime, difftime;
   int *indices,*ranks;
   
#ifdef DIFFGROUPSIZE   
   int choice2;
#endif   

   MPI_Request req;
   MPI_Status status;
   MPI_Group group1,group2;
   MPI_Comm slaves;
   
/* initialize MPI and identify yourself */

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
/*   MPI_Comm_group(MPI_COMM_WORLD,&group1);
   ranks=CREATE(1,int);
   ranks[0]=0;
   MPI_Group_excl(group1,1,&ranks[0],&group2);
   MPI_Comm_create(MPI_COMM_WORLD,group2,&slaves); */
   
/* initialize for SPARSPAK */ 

   spbusr_.MAXSB=12000;
   TOL=0.0000000001;
   TYPTOL=1;

   spksys_.MAXINT=SHRT_MAX;
    
   spbusr_.IERRB=0;
   spbusr_.MSGLVB=1;     
   if ( myid == 0 ) {    /* master */
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
        case 13:
          if ((fp=fopen("wl5016dat","r"))==NULL) {
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
     }
     else if ( (choice1 >= 6) && (choice1<= 13) ) {
        epsilon=1.0e-03; 
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
           case 13:
             fp=fopen("wl5016x","r");
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

     for ( i=0;i<rhs->dim;i++ ) {
       residual->ve[i]=-1.0*rhs->ve[i];
     }

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

     x=v_get(A->n);
        
/* garbage collection */

     free((int *)(A->Jind));
     free((int *)(A->Iind));
     free((int *)(A->colind));
     free((int *)(A->rowptr));
     free((double *)(A->Valuerow));
     free((double *)(A->Valuecol));
       
     r_zero_norm=norm2(residual);
     printf("\n initial norm: %16.14e",r_zero_norm);

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

/* broadcast the number of columns in each block */
     
     buf_count=1;
     MPI_Bcast(&nc,buf_count,MPI_INT,0,MPI_COMM_WORLD);
     
     buf_count=1;
     MPI_Bcast(&ngroup,buf_count,MPI_INT,0,MPI_COMM_WORLD);
     
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
     spbcon_.NZCONS=0;spbusr_.MCOLS=C->n;
     spbusr_.MSEQNS=C->m;
    
     spbcon_.NCOLS=spbusr_.MCOLS;  
     if ( choice1<13 )
       spbusr_.MAXSB=6000;
     else
       spbusr_.MAXSB=10000;  
     T=CREATE(spbusr_.MAXSB,int);
     for (i=0;i<spbusr_.MAXSB;i++) {
       T[i]=-(i+1);
     }
     --T;

/* read in the matrix into SPARSPAK */

     TYPE=1;   
     rnum=1;
     for ( j=1;j<=C->m;j++ ) {
       NSUBS=C->rowptr[j+1]-C->rowptr[j];
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

/* send the nonzero indices of A_id_i's to the slaves */
     
     for (i=1;i<=ngroup;i++ ) {
       tag=12;
       sind=C->Jind[i];
       eind=C->Jind[i+1]-1;
       buf_count=eind-sind+1;
       indices=CREATE(buf_count,int);
       count=0;
       for (j=sind;j<=eind;j++) {
         indices[count]=C->Iind[j];
         count++;
       }
       MPI_Send(&indices[0],buf_count,MPI_INT,i,tag,MPI_COMM_WORLD);
     }
     
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

     D=CREATE(ngroup,VECTOR *);
     --D; 

     DS=CREATE(ngroup,VECTOR *);
     --DS;
     
printf("\n before loop...\n");
     omega=1.0;
     counter=0;
     converged=0;
     position=residual->dim;
     sndvector=CREATE(position,double); 
     
     MPI_Barrier(MPI_COMM_WORLD);   
     starttime=MPI_Wtime();  
     while (converged==0) {
       for (i=1;i<=ngroup;i++ ) {
         MPI_Recv(&sndvector[0],position,MPI_DOUBLE,MPI_ANY_SOURCE,500,MPI_COMM_WORLD,&status); 
         MPI_Get_count(&status,MPI_DOUBLE,&count);
         dest=status.MPI_SOURCE;
         sol=v_get(count);
         for (j=0;j<count;j++) 
           sol->ve[j]=sndvector[j];
         DS[dest]=sol;
       }
       
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
       
       createC222(C,DS,ngroup);
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

/* send the s_i values to each processor */

       for (i=1;i<=ngroup;i++) {
         MPI_Isend(&s->ve[i-1],1,MPI_DOUBLE,i,300,MPI_COMM_WORLD,&req); 
         MPI_Request_free(&req); 
       }
       
/* do a global reduce for the residual */

       buf_count=residual->dim;
       sol=v_get(residual->dim);
       for (i=0;i<residual->dim;i++)
         sol->ve[i]=residual->ve[i];
       MPI_Allreduce(&sol->ve[0],&residual->ve[0],buf_count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);                  

/* find the error */

       r_norm=norm2(residual);

/* printf("\n residual norm: %16.14e\n",r_norm);            */
       
       err=r_norm/r_zero_norm;

       counter++;
       if ((err < epsilon)||(r_norm==0.0)||(counter>15000)) { /* send final message */
         buf_count=1;
         for (i=1;i<=ngroup;i++) {
           MPI_Send(&tag,buf_count,MPI_INT,i,900,MPI_COMM_WORLD);
     /*      MPI_Request_free(&req); */
         }            
         endtime=MPI_Wtime();
         converged=1;
       }
       else {
         
/* garbage collection */

         free((double *)(s->ve));
         free((VECTOR *)(s));         
         ++ddiag;
         ++lnzv;

         free((double *)(ddiag));
         free((double *)(lnzv)); 
/*         free((double *)(C->Valuerow)); */
         free((double *)(C->Valuecol));
         for (i=1;i<=ngroup;i++) {
           free((double *)(DS[i]->ve)); 
         }      
       }                 
     }
     difftime=endtime-starttime;     
     printf("\n Converged after %d loops.",counter);
#ifdef XERROR
     printf("\n The error is: %16.14e \n",err   );
#else   
     printf("\n Two norm of residual is: %16.14e \n",r_norm); 
#endif 
     printf("\n Elapsed time: %lf seconds \n",difftime);    
   }    
   else  {  /* slave */
 /*    MPI_Comm_rank(slaves,&sid); */
   
     A=CREATE(1,MATRIX);
     
     A->rowdefined=TRUE;
     
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

     buf_count=1;
     MPI_Bcast(&nc,buf_count,MPI_INT,0,MPI_COMM_WORLD);
     
     buf_count=1;
     MPI_Bcast(&ngroup,buf_count,MPI_INT,0,MPI_COMM_WORLD);
     
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
     
     sym_mult3(A,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values);   

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
   
/* get the nonzero indices of A_id-i's */

     buf_count=A->m-zeroptr_par[0];
     tag=12;
     indices=CREATE(buf_count,int);
     MPI_Recv(&indices[0],buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);
       
/* begin the solution */

     s=v_get(A->n);
     sol=v_get(A->m);

     tag=500; 
     
     converged=0;
     buf_count=A->m-zeroptr_par[0];
     sndvector=CREATE(buf_count,double);
     x=v_get(nc);
  /*   if ( sid==0) 
       largex=CREATE(nc*ngroup,double); */
       
     MPI_Barrier(MPI_COMM_WORLD); 

     while ( converged != 1 ) {
       s=multAtransr3(A,residual,s);
       s=permdv(A->n,s,invpv);
       s=gsslv(A->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,s);

/* permute the solution to its original order of columns */

       s=permdv(A->n,s,permv);   
       sol=m_fulvectormlt(A,s,sol);
   
/* send both s and sol to master */
       
       eind=A->m-zeroptr_par[0];
       for (j=0;j<eind;j++) 
         sndvector[j]=sol->ve[indices[j]-1];      
       buf_count=eind;
       MPI_Iprobe(0,900,MPI_COMM_WORLD,&flag,&status);
       if ( flag==FALSE ) {
         MPI_Send(&sndvector[0],buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD); 
         
/* get the s_i value */

         MPI_Recv(&ss,1,MPI_DOUBLE,0,300,MPI_COMM_WORLD,&status);

/* compute the update on the residual */

         for (j=0;j<eind;j++)
           sol->ve[indices[j]-1]=sol->ve[indices[j]-1]*ss;
         
/* do a global update on the residual */

         buf_count=residual->dim;
         MPI_Allreduce(&sol->ve[0],&residual->ve[0],buf_count,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);                  
       
/* update the x_i */

         for (j=0;j<x->dim;j++)
           x->ve[j]+=ss*s->ve[j];
         free((double *)(s->ve));
         free((double *)(sol->ve));
         free((VECTOR *)(s));
         free((VECTOR *)(sol));         
         s=v_get(A->n);
         sol=v_get(A->m);
       }
       else {  /* final message; gather the x_vector */
         buf_count=nc*ngroup;
/*         MPI_Gather(&x->ve[0],nc,MPI_DOUBLE,&largex[0],buf_count,MPI_DOUBLE,0,slaves); */
printf("finallll..\n");
         converged=1;
       }
     }
   }   /* end of slave */
   MPI_Finalize();
}    
       
