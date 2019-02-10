/* last revised 20.07.1999                                            */
/* jacobi_pDs.c = routine for solving the LSQ problem iteratively      */
/* and using jacobi updates                                            */
/* this version computes Jacobi with p=DS                              */
/* iterations on a parallel computer synchronously                     */
/***********************************************************************/

/* tag = 100 : used for sending/receiving computed data from slaves    */

#include "mpi.h" 
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include "data_struct.h"
  
VECTOR **D,**Ad,**AD;

/* Ad is used for multiplying A-i with d-i */
/* AD is used for computing the nonzeros of hatA */

MATRIX *hatA,*tildeA;
VECTOR **Apptr;
VECTOR *pvec,**pvecptr,*pd;
MATRIX *AA;
double *deltav,*sndvector;
VECTOR *littles;

/* #define XERROR */  /* to check convergence with error on x */

/* #define FMSTART */ /* starting values of p are F_M values */

#ifdef FMSTART
  VECTOR *evec;
#endif
  
#ifdef XERROR
  VECTOR *xnormv;
#endif  

int totalapsize;
int **indicesptr;

/* variables used for parallel operations */

int myid, numprocs;                         /* to identify the processors */
int dest, tag;                   /* node to send/receive msg, and msg tag */
int buf_count;                       /* the number of elements to be sent */
double *sndvector;
double **sndptr;       
char *buffer;
int position;                                                                    
/**************************************************************************/
int main(argc,argv)
int argc;
char *argv[];
{  

/* local variables */
   
   static int i,j,k,l;
   int rnum, count, counter;
   int choice1,choice2;
   double verdi;
   double number,num,*ip;
   int nc,t,converged,sind,eind;
   double starttime, endtime, difftime;

   MPI_Request *req;
   MPI_Status status,*stats;

/* initialize MPI and identify yourself */

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

/* initialize for SPARSPAK */ 

   spbusr_.MAXSB=12000;
   TOL=0.0000000001;
   TYPTOL=1;
   spksys_.MAXINT=SHRT_MAX;    
   spbusr_.IERRB=0;
   spbusr_.MSGLVB=1;
      
   if ( myid == 0 ) {    /* master */
     sscanf(argv[1],"%d",&choice1);  /* get the test matrix number */
     sscanf(argv[2],"%d",&choice2);  /* get the number of groups */

     switch(choice1) {
        case 1:
          if ((fp=fopen("matrix1","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          if ((fp1=fopen("ash219.d","r"))==NULL) {                 
            printf(" Error (main)        :  open \n");
            return(0);
          }
          break;
        case 2:
          if ((fp=fopen("matrix2","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          if ((fp1=fopen("ash958.d","r"))==NULL) {         
            printf(" Error (main)        :  open \n");
            return(0);
          }
          break;
        case 3:
          if ((fp=fopen("matrix3","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          if ((fp1=fopen("ash331.d","r"))==NULL) {         
            printf(" Error (main)        :  open \n");
            return(0);
          }
          break;
        case 4:
          if ((fp=fopen("matrix4","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          if ((fp1=fopen("ash608.d","r"))==NULL) {          
            printf(" Error (main)        :  open \n");
            return(0);
          }
          break;
        case 5:
          if ((fp=fopen("matrix5","r"))==NULL) {
            printf(" Error (main)        :  open \n");
            return(0);
          }
          if ((fp1=fopen("abb313.d","r"))==NULL) {       
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
     }
     else if ( (choice1 >= 6) && (choice1<= 13) ) {
       if (( choice1 == 6) || (choice1 == 10))
         epsilon=1.0e-03;
       else
         epsilon=1.0e-05;
       readbigmatrix(fp);
     }  
         
     xexact=v_get(A->n);  
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
#ifndef XERROR
     free((double *)(xexact->ve));
     free((VECTOR *)(xexact));      
#endif 
     x=v_get(A->n);  /* starting vector x, all zeros */
     residual=v_get(rhs->dim);  
     for ( i=0;i<rhs->dim;i++ )
       residual->ve[i]=-1.0*rhs->ve[i];

/* garbage collection */

     free((double *)(rhs->ve));
     free((VECTOR *)(rhs)); 
     
/* form the non-overlapping column groups of A */

     ngroup=choice2;
     sub_ptr=CREATE(ngroup,MATRIX *);
     zeroptr=CREATE(ngroup,int *);
     ip=CREATE(1,double);
     number=(double)(A->n)/(double)(ngroup);
     num=modf(number,ip);
     if (num < 0.5)
       nc=floor(number);
     else if (num >= 0.5) 
       nc=ceil(number); 
     free((double *)(ip));
     sub_ptr=create_nover_diff(A,sub_ptr,ngroup,zeroptr,nc);           

/* garbage collection */

     free((int *)(A->Jind));
     free((int *)(A->Iind));
     free((int *)(A->colind));
     free((int *)(A->rowptr));
     free((double *)(A->Valuerow));
     free((double *)(A->Valuecol));
     r_zero_norm=norm2(residual);
     printf("\n initial norm: %16.14e \n",r_zero_norm);
       
/* send each column group to a different processor */

     for ( i=0;i<ngroup;i++ ) {
       AA=sub_ptr[i];
       dest=i+1;
       position=0;
       buf_count=12+4*AA->nonzeros+AA->n+AA->m;
       buf_count=buf_count*sizeof(double);
       buffer=CREATE(buf_count,char);
       j=1;
       MPI_Pack(&AA->m,j,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
       MPI_Pack(&AA->n,j,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
       MPI_Pack(&AA->nonzeros,j,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
       j=AA->nonzeros+1;
       MPI_Pack(&AA->Iind[0],j,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
       j=AA->n+2;
       MPI_Pack(&AA->Jind[0],j,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
       j=AA->m+2;
       MPI_Pack(&AA->rowptr[0],j,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
       j=AA->nonzeros+1;
       MPI_Pack(&AA->colind[0],j,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
       MPI_Pack(&AA->Valuecol[0],j,MPI_DOUBLE,buffer,buf_count,&position,MPI_COMM_WORLD);
       MPI_Pack(&AA->Valuerow[0],j,MPI_DOUBLE,buffer,buf_count,&position,MPI_COMM_WORLD);
       tag=1;
       MPI_Send(buffer,position,MPI_PACKED,dest,tag,MPI_COMM_WORLD);        
       free((char *)(buffer));       
     }

/* broadcast the residual and ngroup */     

     buf_count=1+residual->dim;
     buf_count=buf_count*sizeof(double); 
     buffer=CREATE(buf_count,char);
     position=0;
     MPI_Pack(&ngroup,1,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
     j=residual->dim;
     MPI_Pack(&residual->ve[0],j,MPI_DOUBLE,buffer,buf_count,&position,MPI_COMM_WORLD);
     MPI_Bcast(buffer,position,MPI_PACKED,0,MPI_COMM_WORLD);
     free((char *)(buffer));     
     
/* form A-i*p-i  */

     pvecptr=CREATE(ngroup,VECTOR *);
     for ( l=0; l<ngroup; l++ ) {
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

/* create hatA */

     hatA=CREATE(1,MATRIX);
     hatA->m=residual->dim;     
     hatA->n=ngroup;
     indicesptr=CREATE(ngroup,int *);
     createC1(hatA,ngroup,zeroptr,indicesptr);
     
     Ad=CREATE(ngroup,VECTOR *);
     --Ad;
     Apptr=CREATE(ngroup,VECTOR *);
     totalapsize=0;
     for ( l=0; l<ngroup; l++ ) {
       pvec=pvecptr[l];
       sol=v_get(sub_ptr[l]->m);
       sol=m_fulvectormlt(sub_ptr[l],pvec,sol); 
       count=residual->dim-zeroptr[l][0];
       Apptr[l]=v_get(count);
       totalapsize+=count;
       Ad[l+1]=v_get(count);
       for (j=0;j<count;j++) 
         Apptr[l]->ve[j]=sol->ve[indicesptr[l][j]-1];
       free((double *)(sol->ve));
       free((VECTOR *)(sol));
     }    
     
/* broadcast the size of the Apptr, Apptr and indices pointer */      

     buf_count=ngroup*sizeof(int);
     buffer=CREATE(buf_count,char);
     position=0;
     j=1;
     for (i=0;i<ngroup;i++) {
       MPI_Pack(&Apptr[i]->dim,j,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
     }
     MPI_Bcast(buffer,position,MPI_PACKED,0,MPI_COMM_WORLD);
     free((char *)(buffer));
     
     buf_count=totalapsize*(sizeof(double)+sizeof(int));
     buffer=CREATE(buf_count,char);
     position=0;
     for (i=0;i<ngroup;i++) {
       j=Apptr[i]->dim;
       MPI_Pack(&Apptr[i]->ve[0],j,MPI_DOUBLE,buffer,buf_count,&position,MPI_COMM_WORLD);
       MPI_Pack(&indicesptr[i][0],j,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
     }
     MPI_Bcast(buffer,position,MPI_PACKED,0,MPI_COMM_WORLD);
     free((char *)(buffer));
     
     D=CREATE(ngroup,VECTOR *);
     --D;          
     for (i=1;i<=ngroup;i++ ) {
       j=sub_ptr[i-1]->n;
       D[i]=v_get(j);
     }
     
/* garbage collection */

     for (i=0;i<ngroup;i++) {
       free((int *)(sub_ptr[i]->Jind));
       free((int *)(sub_ptr[i]->Iind));
       free((int *)(sub_ptr[i]->colind));
       free((int *)(sub_ptr[i]->rowptr));
       free((double *)(sub_ptr[i]->Valuerow));
       free((double *)(sub_ptr[i]->Valuecol));
       free((MATRIX *)(sub_ptr[i])); 
     }
     free((MATRIX **)(sub_ptr));     
     
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
     if (choice1<13) 
       spbusr_.MAXSB=6000;
     else
       spbusr_.MAXSB=10000;
     spbmap_.LLFREE=spbusr_.MAXSB;
     spbcon_.NEDGES=0;
     spbcon_.NZMAX=0;
     spbcon_.NZEQNS=0;
     spbcon_.NZCONS=0;
     spbusr_.MCOLS=hatA->n;
     spbusr_.MSEQNS=hatA->m;
     spbcon_.NCOLS=spbusr_.MCOLS;  
     T=CREATE(spbusr_.MAXSB,int);
     for (i=0;i<spbusr_.MAXSB;i++) {
       T[i]=-(i+1);
     }
     --T;
   
/* read in the matrix into SPARSPAK */

     TYPE=1;   
     rnum=1;
     for ( j=1;j<=hatA->m;j++ ) {
       NSUBS=hatA->rowptr[j+1]-hatA->rowptr[j];
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

     counter=0;

#ifdef XERROR
     xnormv=v_get(xexact->dim);
     for (i=0;i<xnormv->dim;i++)
       xnormv->ve[i]=x->ve[i]-xexact->ve[i];
     err=norm2(xnormv);
#else
     r_norm=r_zero_norm;
     err=r_norm/r_zero_norm;   
#endif   
printf("\n before loop...master");   

     sndptr=CREATE(ngroup,double *);
     --sndptr;
     for (i=1;i<=ngroup;i++) {
       buf_count=residual->dim+nc+ngroup;
       sndvector=CREATE(buf_count,double);
       sndptr[i]=sndvector;
     }
     converged=0;   
     stats=CREATE(ngroup,MPI_Status);
     req=CREATE(ngroup,MPI_Request);
     MPI_Barrier(MPI_COMM_WORLD);  
     starttime=MPI_Wtime();  

     while (converged==0) {
       for (i=1;i<=ngroup;i++) {
         MPI_Irecv(sndptr[i],buf_count,MPI_DOUBLE,i,500,MPI_COMM_WORLD,&req[i-1]);           
       }
       deltav=CREATE(ngroup,double);
       --deltav;  
       for (i=1;i<=ngroup;i++) {
         MPI_Waitany(ngroup+1,req,&j,&status);
         dest=status.MPI_SOURCE;
         for (j=1;j<dest;j++)
           deltav[j]+=sndptr[dest][j-1];
         k=0;
         count=D[dest]->dim+dest;
         for (j=dest;j<count;j++) {
           D[dest]->ve[k]=sndptr[dest][j-1];
           k++;
         }
         l=D[dest]->dim+ngroup-1;
         k=1;
         for (j=count;j<=l;j++) {
           deltav[k+dest]+=sndptr[dest][j-1];
           k++;
         }
         t=l+Ad[dest]->dim;
         k=0;
         for (j=l;j<t;j++) {
           Ad[dest]->ve[k]=sndptr[dest][j];
           k++;
         }
       }
          
/* compute hatA^T*hatA and initialize lnz and diag */

       ddiag=CREATE(hatA->n,double);
       lnzv=CREATE(spbcon_.NOFNZ,double);
       --ddiag;
       --lnzv;  
     
/* get the necessary storage for hatA^ThatA */

       i=spbcon_.NOFNZ+hatA->n;
       isub=CREATE(i,int);
       jsub=CREATE(i,int);
       values=CREATE(i,double);
       --isub;
       --jsub;
       --values;

/* compute the columns of A*D=hatA matrix */
   
       for (i=1;i<=ngroup;i++) {
         for (j=0;j<Ad[i]->dim;j++) {
           Ad[i]->ve[j]+=Apptr[i-1]->ve[j]*deltav[i];
         }
       }
       
/* put in the values of hatA */

       createC222(hatA,Ad,ngroup);  
       sym_mult3(hatA,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values);

/* garbage collection */

       isub++;
       jsub++;
       values++;
       free((int *)(isub));
       free((int *)(jsub));
       free((double *)(values));

/* the next three vectors are initialized in gsfct */

       llink=CREATE(hatA->n,int);
       temp=CREATE(hatA->n,double);
       first=CREATE(hatA->n,int);
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
       
/* update the residual and Apptr at the same time */

       for (i=1;i<=ngroup;i++) {
         sind=hatA->Jind[i];
         eind=hatA->Jind[i+1]-1;
         count=0;
         for (j=sind;j<=eind;j++) {
           verdi=Ad[i]->ve[count]*s->ve[i-1];         
           residual->ve[hatA->Iind[j]-1]+=verdi;
           Apptr[i-1]->ve[count]=verdi;
           count++;
         }
       }
       
/* find the error */

#ifdef XERROR
       for (j=0;j<xnormv->dim;j++)
         xnormv->ve[j]=x->ve[j]-xexact->ve[j];
       err=norm2(xnormv);
       r_norm=norm2(residual);

/* printf("\n error: %16.14e\n",err);  
printf("\n residual norm: %16.14e\n",r_norm);         */

#else
       r_norm=norm2(residual);

/* printf("\n residual norm: %16.14e\n",r_norm);        */
       
       err=r_norm/r_zero_norm;
#endif     

/* if not converged then broadcast residual and Apptr */
       
       counter++;
       if ((err < epsilon)||(r_norm==0.0)||(counter>15000)) { /* final */
         converged=1;     
         residual->ve[0]=10000.0;
         position=0;
         buf_count=residual->dim; 
         buf_count=buf_count*sizeof(double); 
         buffer=CREATE(buf_count,char);
         j=1; 
         MPI_Pack(&residual->ve[0],j,MPI_DOUBLE,buffer,buf_count,&position,MPI_COMM_WORLD);
         MPI_Bcast(buffer,position,MPI_PACKED,0,MPI_COMM_WORLD); 
         free((char *)(buffer));       
       }
       else {
         buf_count=residual->dim+totalapsize;
         position=0;
         j=residual->dim;
         buf_count=buf_count*sizeof(double); 
         buffer=CREATE(buf_count,char);
         MPI_Pack(&residual->ve[0],j,MPI_DOUBLE,buffer,buf_count,&position,MPI_COMM_WORLD);
         for (i=0;i<ngroup;i++) {
           j=Apptr[i]->dim;
           MPI_Pack(&Apptr[i]->ve[0],j,MPI_DOUBLE,buffer,buf_count,&position,MPI_COMM_WORLD);
         }
         MPI_Bcast(buffer,position,MPI_PACKED,0,MPI_COMM_WORLD);  

/* garbage collection */

         free((char *)(buffer));
         ++ddiag;
         ++lnzv;
         free((double *)(ddiag));
         free((double *)(lnzv));
       }

/* update x vector and the pvecptr at the same time */
     
       k=0;
       for (i=1;i<=ngroup;i++) {
         t=D[i]->dim;
         for (j=1;j<=t;j++) {
           verdi=s->ve[i-1]*(D[i]->ve[j-1]+deltav[i]*pvecptr[i-1]->ve[j-1]);
           x->ve[k]+=verdi;           
           pvecptr[i-1]->ve[j-1]=verdi;
           k++;       
         }
       }
       free((double *)(s->ve));
       free((VECTOR *)(s)); 
       deltav++;
       free((double *)(deltav));          
     }
     endtime=MPI_Wtime();
     difftime=endtime-starttime;
     fp=fopen("results_par_pDs","a");
     fprintf(fp,"\n Converged after %d loops.",counter);

#ifdef XERROR
     fprintf(fp,"\n The error is: %16.14e \n",err   );
#else   
     fprintf(fp,"\n Two norm of residual is: %16.14e \n",r_norm); 
#endif
   
     fprintf(fp,"\n Elapsed time: %lf seconds \n",difftime); 
     fclose(fp);   
   } /* end of master */
   else {  /* slaves */
     A=CREATE(1,MATRIX);
     A->rowdefined=TRUE;
     tag=1;
     buf_count=100000;
     buf_count=buf_count*sizeof(double); 
     buffer=CREATE(buf_count,char);
     MPI_Recv(buffer,buf_count,MPI_PACKED,0,tag,MPI_COMM_WORLD,&status);
     MPI_Unpack(buffer,buf_count,&position,&A->m,1,MPI_INT,MPI_COMM_WORLD);
     MPI_Unpack(buffer,buf_count,&position,&A->n,1,MPI_INT,MPI_COMM_WORLD);
     MPI_Unpack(buffer,buf_count,&position,&A->nonzeros,1,MPI_INT,MPI_COMM_WORLD);
     j=A->nonzeros+1;
     A->Iind=CREATE(j,int);
     MPI_Unpack(buffer,buf_count,&position,&A->Iind[0],j,MPI_INT,MPI_COMM_WORLD);
     j=A->n+2;
     A->Jind=CREATE(j,int);
     MPI_Unpack(buffer,buf_count,&position,&A->Jind[0],j,MPI_INT,MPI_COMM_WORLD);
     j=A->m+2;
     A->rowptr=CREATE(j,int);
     MPI_Unpack(buffer,buf_count,&position,&A->rowptr[0],j,MPI_INT,MPI_COMM_WORLD);
     j=A->nonzeros+1;
     A->colind=CREATE(j,int);
     MPI_Unpack(buffer,buf_count,&position,&A->colind[0],j,MPI_INT,MPI_COMM_WORLD);
     j=A->nonzeros+1;
     A->Valuecol=CREATE(j,double);
     MPI_Unpack(buffer,buf_count,&position,&A->Valuecol[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
     A->Valuerow=CREATE(j,double);
     MPI_Unpack(buffer,buf_count,&position,&A->Valuerow[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
     free((char *)(buffer));     

     buf_count=1+A->m;
     buf_count=buf_count*sizeof(double); 
     buffer=CREATE(buf_count,char);
     position=0;
     MPI_Bcast(buffer,buf_count,MPI_PACKED,0,MPI_COMM_WORLD);
     MPI_Unpack(buffer,buf_count,&position,&ngroup,1,MPI_INT,MPI_COMM_WORLD);
     residual=v_get(A->m);
     j=residual->dim;
     MPI_Unpack(buffer,buf_count,&position,&residual->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
     free((char *)(buffer));

/* broadcast the size of the Apptr, Apptr and indices pointer */      

     buf_count=ngroup*sizeof(int);
     buffer=CREATE(buf_count,char);
     MPI_Bcast(buffer,buf_count,MPI_PACKED,0,MPI_COMM_WORLD);
     Apptr=CREATE(ngroup,VECTOR *);
     indicesptr=CREATE(ngroup,int *);
     position=0;
     totalapsize=0;
     for (i=0;i<ngroup;i++) {
       MPI_Unpack(buffer,buf_count,&position,&count,1,MPI_INT,MPI_COMM_WORLD);
       Apptr[i]=v_get(count);
       indicesptr[i]=CREATE(count,int);
       totalapsize+=count;
     }
     free((char *)(buffer));
     
     buf_count=totalapsize*(sizeof(double)+sizeof(int));
     buffer=CREATE(buf_count,char);
     MPI_Bcast(buffer,buf_count,MPI_PACKED,0,MPI_COMM_WORLD);
     position=0;
     for (i=0;i<ngroup;i++) {
       j=Apptr[i]->dim;
       MPI_Unpack(buffer,buf_count,&position,&Apptr[i]->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
       MPI_Unpack(buffer,buf_count,&position,&indicesptr[i][0],j,MPI_INT,MPI_COMM_WORLD);
     }
     free((char *)(buffer));
               
/* form tildeA */

     tildeA=CREATE(1,MATRIX);
     tildeA=createtildeA(tildeA,A,Apptr,ngroup,myid,indicesptr);

/* do symbolic factorization on tildeA */
     
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
     spbusr_.MCOLS=tildeA->n;
     spbusr_.MSEQNS=tildeA->m;
     spbcon_.NCOLS=spbusr_.MCOLS;  
     T=CREATE(spbusr_.MAXSB,int);
     for (i=0;i<spbusr_.MAXSB;i++) {
       T[i]=-(i+1);
     }
     
     --T;

/* read in the matrix into SPARSPAK */

     TYPE=1;   
     rnum=1;
    
     for ( j=1;j<=tildeA->m;j++ ) {
       NSUBS=tildeA->rowptr[j+1]-tildeA->rowptr[j];
       spbcon_.NOFNZ+=NSUBS;
       ROWNUM=rnum;
       rnum++;
       SUBS=CREATE(NSUBS,int);
       for ( k=0;k<NSUBS;k++ ) {
         *(SUBS+k)=tildeA->colind[tildeA->rowptr[j]+k];
       }
       --SUBS;
       inxywb(ROWNUM,TYPE,NSUBS,SUBS,T);
       SUBS++;
       free((int *)(SUBS)); 
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
     j = tildeA->n+1;
     for (i = 1; i <= j; i++) {
       xrnzv[i] = T[xrnz+i];
       xnzsubv[i] = T[xnzsub+i];
     }
     for (i = 1; i <= spbcon_.NOFNZ; i++) {
	nzsubv[i] = T[nzsub+i];
     }
     j = tildeA->n;
     for (i = 1; i <= j; i++) {
       permv[i] = T[perm+i];
       invpv[i] = T[invp+i];
     }

/* garbage collection */

     ++T;
     free((int *)(T));  
     
/* begin iteration */     

     counter=0;
     count=Apptr[myid-1]->dim+tildeA->n;
     sndvector=CREATE(count,double);
     littles=v_get(A->n);
     buf_count=residual->dim+totalapsize;
     buf_count=buf_count*sizeof(double); 
     buffer=CREATE(buf_count,char);
          
     MPI_Barrier(MPI_COMM_WORLD);  

     while (counter==0) {

/* factorize tildeA */

       ddiag=CREATE(tildeA->n,double);
       lnzv=CREATE(spbcon_.NOFNZ,double);
       --ddiag;
       --lnzv;
       j=spbcon_.NOFNZ+tildeA->n;
       isub=CREATE(j,int);
       jsub=CREATE(j,int);
       values=CREATE(j,double);
       --isub;
       --jsub;
       --values;  
       sym_mult3(tildeA,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values);

/* garbage collection */

       isub++;
       jsub++;
       values++;
       free((int *)(isub));
       free((int *)(jsub));
       free((double *)(values));
   
/* the next three vectors are initialized in gsfct */
   
       llink=CREATE(tildeA->n,int);
       temp=CREATE(tildeA->n,double);
       first=CREATE(tildeA->n,int);
       --llink;
       --temp;
       --first;
       gsfct(tildeA->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,llink,first,temp,flag__); 

/* garbage collection */
      
       ++llink;
       ++temp;
       ++first;
       free((int *)(llink));
       free((double *)(temp));
       free((int *)(first));
                       
       s=v_get(tildeA->n);
       s=multAtransr3(tildeA,residual,s);
       s=permdv(tildeA->n,s,invpv);
       s=gsslv(tildeA->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,s);

/* permute the solution to its original order of colums */

       s=permdv(tildeA->n,s,permv);          

/* copy this s to the vector to be sent to the master */ 
     
       for (j=0;j<s->dim;j++)
         sndvector[j]=s->ve[j];

/* take the littles out of s */

       for (j=0;j<A->n;j++) {
         littles->ve[j]=s->ve[myid-1+j];
       }                
       
/* compute A*d */

       sol=v_get(A->m);
       sol=m_fulvectormlt(A,littles,sol);    
       
/* send s and the nonzero entries of sol to master */
/* get the nonzeros of sol */
         
       eind=Apptr[myid-1]->dim;
       for (j=0;j<eind;j++) 
         sndvector[j+s->dim]=sol->ve[indicesptr[myid-1][j]-1];
       tag=500;
       buf_count=eind+s->dim;
       MPI_Send(sndvector,buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);   
       
/* garbage collection */

       ++lnzv;
       ++ddiag;
       free((double *)(lnzv));
       free((double *)(ddiag));
       free((double *)(sol->ve));
       free((VECTOR *)(sol));
       free((double *)(s->ve));
       free((VECTOR *)(s));
       buf_count=residual->dim+totalapsize;
       buf_count=buf_count*sizeof(double); 
       MPI_Bcast(buffer,buf_count,MPI_PACKED,0,MPI_COMM_WORLD);
       j=residual->dim;
       position=0;
       MPI_Unpack(buffer,buf_count,&position,&residual->ve[0],1,MPI_DOUBLE,MPI_COMM_WORLD);      
       if (residual->ve[0]!=10000.0) {
         j=residual->dim-1;
         MPI_Unpack(buffer,buf_count,&position,&residual->ve[1],j,MPI_DOUBLE,MPI_COMM_WORLD);
         for (i=0;i<ngroup;i++) {
           j=Apptr[i]->dim;
           MPI_Unpack(buffer,buf_count,&position,&Apptr[i]->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
         }
         i=myid;
         putvaluesintildeA3(tildeA,Apptr,A,ngroup,i);  
       }
       else { /* final */ 
         counter=1;                  
       }
     }
   }
   MPI_Finalize(); 
}    
       
