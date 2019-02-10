/* last revised 12.07.1999                                              */
/* jacobi_predict_new.c = routine for solving the LSQ problem iteratively  */
/* and using jacobi updates                                            */
/* this version is to test the predictions with asynchronous communication */
/* and synchronization after predictions. in between two updates of the */
/* residual and x vector a temporary copy of the residual is used in the */
/* computations.                                                        */
/***********************************************************************/

/* tag = 100 : used for sending/receiving computed data from slaves    */

#include "mpi.h" 
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include "data_struct.h"

VECTOR **D,**Ad;

/* Ad is used for multiplying A-i with d-i */

MATRIX *hatA,*tildeA;
VECTOR **Apptr;
VECTOR *pvec,**pvecptr;
MATRIX *AA;
VECTOR *littles;
int *indices;

/* #define XERROR */   /* to check convergence with error on x */

/* #define FMSTART */ /* starting values of p are F_M values */

#ifdef FMSTART
  VECTOR *evec;
#endif

VECTOR *vvec;
int antall;
int totalapsize;
  
#ifdef XERROR
  VECTOR *xnormv;
#endif    

/* variables used for parallel operations */

int myid, numprocs;                         /* to identify the processors */
int dest, tag;                   /* node to send/receive msg, and msg tag */
int buf_count;                       /* the number of elements to be sent */
double *sndvector,**sndptr;
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
   int nc,t,converged;
   double starttime, endtime, difftime;
   int sind,eind;
   
   MPI_Request request,*req;
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
     sscanf(argv[3],"%d",&antall);  /* get the number of prediction steps */
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
       buf_count=13+4*AA->nonzeros+AA->n+2*AA->m;
       buf_count=buf_count*sizeof(double); 
       buffer=CREATE(buf_count,char);
       j=1;
       MPI_Pack(&AA->m,j,MPI_INT,sndvector,buf_count,&position,MPI_COMM_WORLD);
       MPI_Pack(&AA->n,j,MPI_INT,sndvector,buf_count,&position,MPI_COMM_WORLD);
       MPI_Pack(&AA->nonzeros,j,MPI_INT,sndvector,buf_count,&position,MPI_COMM_WORLD);
       j=AA->nonzeros+1;
       MPI_Pack(&AA->Iind[0],j,MPI_INT,sndvector,buf_count,&position,MPI_COMM_WORLD);
       j=AA->n+2;
       MPI_Pack(&AA->Jind[0],j,MPI_INT,sndvector,buf_count,&position,MPI_COMM_WORLD);
       j=AA->m+2;
       MPI_Pack(&AA->rowptr[0],j,MPI_INT,sndvector,buf_count,&position,MPI_COMM_WORLD);
       j=AA->nonzeros+1;
       MPI_Pack(&AA->colind[0],j,MPI_INT,sndvector,buf_count,&position,MPI_COMM_WORLD);
       MPI_Pack(&AA->Valuecol[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
       MPI_Pack(&AA->Valuerow[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
       zeroptr_par=CREATE(AA->m+1,int);
       for ( j=0;j<=AA->m;j++ ) 
         zeroptr_par[j]=zeroptr[i][j];
       j=AA->m+1;
       MPI_Pack(&zeroptr_par[0],j,MPI_INT,sndvector,buf_count,&position,MPI_COMM_WORLD);
       tag=1;
       MPI_Isend(buffer,position,MPI_PACKED,dest,tag,MPI_COMM_WORLD,&request);        
       free((int *)(zeroptr_par));
       free((double *)(sndvector));
     }

/* broadcast the number of columns in each block, number of groups */
/* number of predictions and the residual */
     
     buf_count=3+residual->dim;
     buf_count=buf_count*sizeof(double); 
     buffer=CREATE(buf_count,char);
     position=0;
     MPI_Pack(&nc,1,MPI_INT,sndvector,buf_count,&position,MPI_COMM_WORLD);
     MPI_Pack(&ngroup,1,MPI_INT,sndvector,buf_count,&position,MPI_COMM_WORLD);
     MPI_Pack(&antall,1,MPI_INT,sndvector,buf_count,&position,MPI_COMM_WORLD);
     j=residual->dim;
     MPI_Pack(&residual->ve[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
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
     Apptr=CREATE(ngroup,VECTOR *);
     for ( l=0; l<ngroup; l++ ) {
       pvec=pvecptr[l];
       sol=v_get(sub_ptr[l]->m);
       sol=m_fulvectormlt(sub_ptr[l],pvec,sol); 
       Apptr[l]=sol;
     }
     
/* broadcast the Apptr */

     buf_count=residual->dim*ngroup;
     buf_count=buf_count*sizeof(double); 
     buffer=CREATE(buf_count,char);
     position=0;
     j=residual->dim;
     for (i=0;i<ngroup;i++) {
       MPI_Pack(&Apptr[i]->ve[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
     } 
     MPI_Bcast(buffer,position,MPI_PACKED,0,MPI_COMM_WORLD);
     free((double *)(buffer));
     
/* do the symmetric factorization of hatA */

     hatA=CREATE(1,MATRIX);
     hatA->m=residual->dim;     
     hatA->n=ngroup;
     createC1(sub_ptr,hatA,ngroup,nc,zeroptr);
     D=CREATE(ngroup,VECTOR *);
     --D;
     for (i=1;i<=ngroup;i++ ) 
       D[i]=v_get(sub_ptr[i-1]->n);
       
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
     Ad=CREATE(ngroup,VECTOR *);
     --Ad;
     
/* send the nonzero indices of A_id_i's to the slaves */
     
     req=CREATE(ngroup+1,MPI_Request);
     for (i=1;i<=ngroup;i++ ) {
       tag=12;
       sind=hatA->Jind[i];
       eind=hatA->Jind[i+1]-1;
       buf_count=eind-sind+1;
       Ad[i]=v_get(buf_count);
       indices=CREATE(buf_count,int);
       count=0;
       for (j=sind;j<=eind;j++) {
         indices[count]=hatA->Iind[j];
         count++;
       }
       MPI_Isend(&indices[0],buf_count,MPI_INT,i,tag,MPI_COMM_WORLD,&req[i]);
       free((int *)(indices));
     }
     for (i=1;i<=ngroup;i++)
       MPI_Request_free(&req[i]);
       
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
     counter=0;
     sndptr=CREATE(ngroup,double *);
     --sndptr;
     for (i=1;i<=ngroup;i++) 
        sndptr[i]=CREATE(100000,double);
printf("before loop...\n");   
     converged=0;   

/* initialize the temporary residual */

     vvec=v_get(residual->dim);
     for (i=0;i<vvec->dim;i++)
       vvec->ve[i]=residual->ve[i];
     MPI_Barrier(MPI_COMM_WORLD);   
     starttime=MPI_Wtime();  
     while (converged==0) {

       for (t=1;t<=antall;t++) {
         buf_count=100000;
         for (i=1;i<=ngroup;i++) {
           MPI_Irecv(&sndptr[i][0],buf_count,MPI_DOUBLE,i,500,MPI_COMM_WORLD,&req[i]);
         }
         for (i=1;i<=ngroup;i++) {
           MPI_Waitany(ngroup+1,req,&j,&status);
           dest=status.MPI_SOURCE;         
printf(" %d\n",dest);           
           count=D[dest]->dim;
           for (j=0;j<count;j++)
             D[dest]->ve[j]+=sndptr[dest][j];
           sind=hatA->Jind[dest];
           eind=Ad[dest]->dim;
           for (j=0;j<eind;j++) {
             verdi=sndptr[dest][j+count];
             Ad[dest]->ve[j]+=verdi;
             vvec->ve[hatA->Iind[sind+j]-1]+=verdi;
           }
           if (t<antall) {
             buf_count=residual->dim;
             MPI_Isend(&vvec->ve[0],buf_count,MPI_DOUBLE,dest,600,MPI_COMM_WORLD,&request);   
             MPI_Request_free(&request);        
           }
           else {
             if (i<ngroup) {
               buf_count=residual->dim;
               MPI_Isend(&vvec->ve[0],buf_count,MPI_DOUBLE,dest,600,MPI_COMM_WORLD,&request); 
               MPI_Request_free(&request);  
             }      
           }
         }        
       }

/* do the planar search */
         
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
       gsfct(hatA->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,llink,first,temp,flag__); 

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
           Apptr[i-1]->ve[hatA->Iind[j]-1]=verdi;
           count++;
         }
       }
           
/* find the error */

#ifdef XERROR
       for (i=0;i<xnormv->dim;i++)
         xnormv->ve[i]=x->ve[i]-xexact->ve[i];
       err=norm2(xnormv);
       r_norm=norm2(residual);
printf("\n error: %16.14e\n",err);  
printf("\n residual norm: %16.14e\n",r_norm);        
#else
       r_norm=norm2(residual);
printf("\n residual norm: %16.14e\n",r_norm);       
       err=r_norm/r_zero_norm;
#endif     

/* if not converged then broadcast residual and A_ip_i*/
       
       counter++;      
       if ((err < epsilon)||(r_norm==0.0)||(counter>15000)) { /* final */
         converged=1;     
         residual->ve[0]=10000.0;
         position=0;
         buf_count=residual->dim;
         buf_count=buf_count*sizeof(double); 
         sndvector=CREATE(100000,double);
         j=residual->dim;
         MPI_Pack(&residual->ve[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
         MPI_Bcast(sndvector,position,MPI_PACKED,0,MPI_COMM_WORLD); 
         free((double *)(sndvector));
       }
       else {
         buf_count=100000;
         position=0;
         j=residual->dim;
         buf_count=buf_count*sizeof(double); 
         sndvector=CREATE(100000,double);
         MPI_Pack(&residual->ve[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
         for (i=0;i<ngroup;i++) {
           MPI_Pack(&Apptr[i]->ve[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
         }
         MPI_Bcast(sndvector,position,MPI_PACKED,0,MPI_COMM_WORLD);

/* garbage collection */

         free((double *)(sndvector));
         ++ddiag;
         ++lnzv;
         free((double *)(ddiag));
         free((double *)(lnzv));
         for (i=1;i<=ngroup;i++) {
           for (j=0;j<Ad[i]->dim;j++)
             Ad[i]->ve[j]=0.0;
         }

/* update the vvec */

         for (i=0;i<residual->dim;i++)
           vvec->ve[i]=residual->ve[i];            
       }

/* update x vector */

       k=0;
       for (i=1;i<=ngroup;i++) {
         for (j=0;j<D[i]->dim;j++) {
           x->ve[k]+=D[i]->ve[j]*s->ve[i-1];
           k++;
         }
       }

       free((double *)(s->ve));
       free((VECTOR *)(s)); 
       for (i=1;i<=ngroup;i++) {
         for (j=0;j<D[i]->dim;j++)
           D[i]->ve[j]=0.0;
       }         
     }
     endtime=MPI_Wtime();
     difftime=endtime-starttime;
     printf("\n Converged after %d loops.",counter);
#ifdef XERROR
     printf("\n The error is: %16.14e \n",err   );
#else   
     printf("\n Two norm of residual is: %16.14e \n",r_norm); 
#endif
     printf("\n Elapsed time: %lf seconds \n",difftime);    
   } /* end of master */
   else {  /* slaves */
     A=CREATE(1,MATRIX);
     A->rowdefined=TRUE;
     tag=1;
     buf_count=100000;
     sndvector=CREATE(buf_count,double);
     buf_count=buf_count*sizeof(double); 
     MPI_Recv(&sndvector[0],buf_count,MPI_PACKED,0,tag,MPI_COMM_WORLD,&status);
     MPI_Unpack(sndvector,buf_count,&position,&A->m,1,MPI_INT,MPI_COMM_WORLD);
     MPI_Unpack(sndvector,buf_count,&position,&A->n,1,MPI_INT,MPI_COMM_WORLD);
     MPI_Unpack(sndvector,buf_count,&position,&A->nonzeros,1,MPI_INT,MPI_COMM_WORLD);
     j=A->nonzeros+1;
     A->Iind=CREATE(j,int);
     MPI_Unpack(sndvector,buf_count,&position,&A->Iind[0],j,MPI_INT,MPI_COMM_WORLD);
     j=A->n+2;
     A->Jind=CREATE(j,int);
     MPI_Unpack(sndvector,buf_count,&position,&A->Jind[0],j,MPI_INT,MPI_COMM_WORLD);
     j=A->m+2;
     A->rowptr=CREATE(j,int);
     MPI_Unpack(sndvector,buf_count,&position,&A->rowptr[0],j,MPI_INT,MPI_COMM_WORLD);
     j=A->nonzeros+1;
     A->colind=CREATE(j,int);
     MPI_Unpack(sndvector,buf_count,&position,&A->colind[0],j,MPI_INT,MPI_COMM_WORLD);
     j=A->nonzeros+1;
     A->Valuecol=CREATE(j,double);
     MPI_Unpack(sndvector,buf_count,&position,&A->Valuecol[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
     A->Valuerow=CREATE(j,double);
     MPI_Unpack(sndvector,buf_count,&position,&A->Valuerow[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
     j=A->m+1;
     zeroptr_par=CREATE(j,int);
     MPI_Unpack(sndvector,buf_count,&position,&zeroptr_par[0],j,MPI_INT,MPI_COMM_WORLD);
     free((double *)(sndvector));
     buf_count=3+A->m;
     sndvector=CREATE(buf_count,double);
     buf_count=buf_count*sizeof(double); 
     position=0;
     MPI_Bcast(sndvector,buf_count,MPI_PACKED,0,MPI_COMM_WORLD);
     MPI_Unpack(sndvector,buf_count,&position,&nc,1,MPI_INT,MPI_COMM_WORLD);
     MPI_Unpack(sndvector,buf_count,&position,&ngroup,1,MPI_INT,MPI_COMM_WORLD);
     MPI_Unpack(sndvector,buf_count,&position,&antall,1,MPI_INT,MPI_COMM_WORLD);
     residual=v_get(A->m);
     j=residual->dim;
     MPI_Unpack(sndvector,buf_count,&position,&residual->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
     Apptr=CREATE(ngroup,VECTOR *);
     free((double *)(sndvector));
     
/* broadcast the Apptr */

     buf_count=100000;
     sndvector=CREATE(buf_count,double);
     buf_count=buf_count*sizeof(double); 
     position=0;
     j=A->m;
     MPI_Bcast(sndvector,buf_count,MPI_PACKED,0,MPI_COMM_WORLD);
     for (i=0;i<ngroup;i++) {
       Apptr[i]=v_get(j);
       MPI_Unpack(sndvector,buf_count,&position,&Apptr[i]->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
     }     
     free((double *)(sndvector));
     buf_count=A->m;
     indices=CREATE(buf_count,int);
     buf_count=buf_count*sizeof(double); 
     tag=12;
     MPI_Irecv(&indices[0],buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&request);

/* form tildeA */

     i=myid;
     tildeA=CREATE(1,MATRIX);
     tildeA=createtildeA(tildeA,A,Apptr,ngroup,i);

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
       if ( tildeA->rowptr[j] != 0 ) { 
         count=1;
         while ( tildeA->rowptr[j+count] == 0 )
           count++;
         NSUBS=tildeA->rowptr[j+count]-tildeA->rowptr[j];
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
     MPI_Wait(&request,&status);
    
/* begin iteration */  

     counter=0;
     sndvector=CREATE(100000,double);
     littles=v_get(A->n);
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

/* using the same tildeA matrix solve the system antall times with different */
/* residual values received from the master processor */
       
       for (k=1;k<=antall;k++) {  
         s=v_get(tildeA->n);
         s=multAtransr3(tildeA,residual,s);
         s=permdv(tildeA->n,s,invpv);
         s=gsslv(tildeA->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,s);
      
/* permute the solution to its original order of columns */

         s=permdv(tildeA->n,s,permv);          

/* take the d_i vector out of d_itilde  */

  /*       sndvector=CREATE(100000,double);
         littles=v_get(A->n); */
         for (j=0;j<A->n;j++) {
           verdi=s->ve[myid-1+j];
           sndvector[j]=verdi;
           littles->ve[j]=verdi;
         }                
    
/* compute A*d */

         sol=v_get(A->m);
         sol=m_fulvectormlt(A,littles,sol);    
       
/* send s and the nonzero entries of sol to master */
/* get the nonzeros of sol */
         
         eind=A->m-zeroptr_par[0];
         for (j=0;j<eind;j++) 
           sndvector[j+A->n]=sol->ve[indices[j]-1];
         tag=500;
         buf_count=A->n+eind;
         MPI_Isend(sndvector,buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&request);        
       
/* garbage collection */

         free((double *)(sol->ve));
         free((VECTOR *)(sol));
         free((double *)(s->ve));
         free((VECTOR *)(s));
  /*       free((double *)(littles->ve));
         free((VECTOR *)(littles)); */
         MPI_Wait(&request,&status);
         if (k < antall) {  /* receive the new residual only */
           buf_count=residual->dim;
           MPI_Recv(&residual->ve[0],buf_count,MPI_DOUBLE,0,600,MPI_COMM_WORLD,&status);
         }
         else {
           ++ddiag;
           ++lnzv;
           free((double *)(ddiag));
           free((double *)(lnzv));
         }
       }
       position=0;
       buf_count=100000;
       buf_count=buf_count*sizeof(double); 
/*       free((double *)(sndvector));
       sndvector=CREATE(100000,double); */
       MPI_Bcast(sndvector,buf_count,MPI_PACKED,0,MPI_COMM_WORLD);
       j=residual->dim;
       MPI_Unpack(sndvector,buf_count,&position,&residual->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
       if (residual->ve[0]!=10000.0) {
         j=A->m;
         for (i=0;i<ngroup;i++)
           MPI_Unpack(sndvector,buf_count,&position,&Apptr[i]->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
         i=myid;
         putvaluesintildeA3(tildeA,Apptr,A,ngroup,i); 
    /*     free((double *)(sndvector)); */
       }
       else { /* final */ 
         counter=1;                  
       }
     }
   }
   MPI_Finalize(); 
}    
