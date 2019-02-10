/* last revised 26.07.1999                                         */
/* jacobi_p0.c = routine for solving the LSQ problem synchronously    */
/* and using Jacobi updates and planar search                      */
/*******************************************************************/
#include "mpi.h"
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include "data_struct.h"

/* #define XERROR */  /* to check convergence with error on x */ 

/* #define LINESEARCH */

#ifdef XERROR
  VECTOR *xnormv;
#endif

MATRIX *C;
VECTOR **D;
VECTOR **Ad;
int **indicesptr;
  
/* variables used for parallel operations */

int myid, numprocs;                         /* to identify the processors */
int dest, tag;                   /* node to send/receive msg, and msg tag */
int buf_count;                       /* the number of elements to be sent */
double *sndvector,**sndptr;
int position;                                      
char *buffer; 

/**************************************************************************/
int main(argc,argv)
int argc;
char *argv[];
{  

/* local variables */
   
   static int i,j,k;
   int rnum, count, counter;
   MATRIX *AA;
   int converged;
   double *ip,number,num;
   int nc,sind,eind,choice1,choice2;
   double starttime, endtime, difftime;
   int *indices;
   int indicessize;
   double verdi;

#ifdef LINESEARCH
   double check,dotnorm;
#endif
      
   MPI_Request *req;
   MPI_Status status;
   
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
     sscanf(argv[1],"%d",&choice1); 
     sscanf(argv[2],"%d",&choice2); /* get the number of groups */

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
  
     x=v_get(A->n);
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
       
/* broadcast the residual and the number of columns in each block */     

     buf_count=residual->dim+1;
     buf_count=buf_count*sizeof(double); 
     buffer=CREATE(buf_count,char);
     position=0;
     MPI_Pack(&nc,1,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
     j=residual->dim;
     MPI_Pack(&residual->ve[0],j,MPI_DOUBLE,buffer,buf_count,&position,MPI_COMM_WORLD);
     MPI_Bcast(buffer,position,MPI_PACKED,0,MPI_COMM_WORLD);
     free((char *)(buffer));     

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
     
/* create hatA */

     C=CREATE(1,MATRIX);
     C->m=residual->dim;     
     C->n=ngroup;
     indicesptr=CREATE(ngroup,int *);     
     createC1(C,ngroup,zeroptr,indicesptr);

/* send the indices of nonzeros in A_id_i to the slaves */
     
     Ad=CREATE(ngroup,VECTOR *);
     --Ad;
     for ( i=1;i<=ngroup;i++ ) {
       count=count=residual->dim-zeroptr[i-1][0];
       Ad[i]=v_get(count);
       buf_count=1+count;
       buf_count=buf_count*sizeof(int);
       buffer=CREATE(buf_count,char);
       position=0;
       MPI_Pack(&count,1,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
       MPI_Pack(&indicesptr[i-1][0],count,MPI_INT,buffer,buf_count,&position,MPI_COMM_WORLD);
       MPI_Send(buffer,position,MPI_PACKED,i,20,MPI_COMM_WORLD);
     }     

/* garbage collection */

     for (i=0;i<ngroup;i++) 
       free((int *)(zeroptr[i]));
     free((int **)(zeroptr));
     
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
     spbusr_.MCOLS=C->n;
     spbusr_.MSEQNS=C->m;
     spbcon_.NCOLS=spbusr_.MCOLS; 
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
     
#ifdef XERROR
     xnormv=v_get(xexact->dim);
     for (i=0;i<xnormv->dim;i++)
       xnormv->ve[i]=x->ve[i]-xexact->ve[i];
     err=norm2(xnormv);
#else
     r_norm=r_zero_norm;
     err=r_norm/r_zero_norm;   
#endif         

     
printf("\n before loop...\n");
     counter=0;
     converged=0;
     buf_count=nc*residual->dim;
     sndptr=CREATE(ngroup,double *);
     --sndptr;
     for (i=1;i<=ngroup;i++) {
       sndvector=CREATE(buf_count,double);
       sndptr[i]=sndvector;
     }
     req=CREATE(count,MPI_Request);
     MPI_Barrier(MPI_COMM_WORLD);   
     starttime=MPI_Wtime();  
     while (converged==0) {
       buf_count=nc*residual->dim;
       for (i=1;i<=ngroup;i++ ) {
         MPI_Irecv(sndptr[i],buf_count,MPI_DOUBLE,i,500,MPI_COMM_WORLD,&req[i-1]);           
       }         
       for (i=1;i<=ngroup;i++ ) {
         MPI_Waitany(ngroup,req,&j,&status);
         dest=status.MPI_SOURCE;
         count=D[dest]->dim;
         for (j=0;j<count;j++) 
           D[dest]->ve[j]=sndptr[dest][j];
         for (j=0;j<Ad[dest]->dim;j++)
           Ad[dest]->ve[j]=sndptr[dest][j+count];
       }
#ifdef LINESEARCH
       s=v_get(ngroup);
       for (i=1;i<=ngroup;i++) { 
         check=0.0;
         sind=C->Jind[i];
         eind=C->Jind[i+1]-1;
         count=0;
         dotnorm=0.0;
         for (j=sind;j<=eind;j++) {
           verdi=Ad[i]->ve[count];
           check+=verdi*residual->ve[C->Iind[j]-1];
           dotnorm+=verdi*verdi;
           count++;
         }
         s->ve[i-1]=-check/dotnorm; 
       }
#else          
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
       
       createC222(C,Ad,ngroup);
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
#endif
     
/* update the residual vector */

       for (i=1;i<=ngroup;i++) {
         sind=C->Jind[i];
         eind=C->Jind[i+1]-1;
         count=0;
         for (j=sind;j<=eind;j++) {
           residual->ve[C->Iind[j]-1]+=Ad[i]->ve[count]*s->ve[i-1];
           count++;
         }
       }
        
/* find the error */

#ifdef XERROR
       for (i=0;i<xnormv->dim;i++)
         xnormv->ve[i]=x->ve[i]-xexact->ve[i];
       err=norm2(xnormv);
       r_norm=norm2(residual);

/* printf("\n error: %16.14e\n",err);  
printf("\n residual norm: %16.14e\n",r_norm);   */     

#else
       r_norm=norm2(residual);

/* printf("\n residual norm: %16.14e\n",r_norm);           */
       
       err=r_norm/r_zero_norm;
#endif            

       counter++;
       if ((err < epsilon)||(r_norm==0.0)||(counter>15000)) { /* send final message */
         converged=1;     
         residual->ve[0]=10000.0;
         buf_count=1; 
         MPI_Bcast(&residual->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD); 
         free((char *)(buffer));
       }
       else {
         buf_count=residual->dim;                
         MPI_Bcast(&residual->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD); 
       
/* garbage collection */

         free((double *)(s->ve));
         free((VECTOR *)(s)); 
#ifndef LINESEARCH                 
         ++ddiag;
         ++lnzv;
         free((double *)(ddiag));
         free((double *)(lnzv)); 
#endif         
       }  

/* update x vector */

       k=0;
       for (i=1;i<=ngroup;i++) {
         for (j=1;j<=D[i]->dim;j++) {
           x->ve[k]+=D[i]->ve[j-1]*s->ve[i-1];
           k++;
         }
       }                      
     }
     endtime=MPI_Wtime();     
     difftime=endtime-starttime;  
     fp=fopen("results_par_p0","a");    
     fprintf(fp,"\n Converged after %d loops.",counter);
#ifdef XERROR
     fprintf(fp,"\n The error is: %16.14e \n",err   );
#else   
     fprintf(fp,"\n Two norm of residual is: %16.14e \n",r_norm); 
#endif 
     fprintf(fp,"\n Elapsed time: %16.14e seconds \n",difftime); 
     fclose(fp);    
   }    
   else  {  /* slave */
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
     MPI_Unpack(buffer,buf_count,&position,&nc,1,MPI_INT,MPI_COMM_WORLD);
     residual=v_get(A->m);
     j=residual->dim;
     MPI_Unpack(buffer,buf_count,&position,&residual->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
     free((char *)(buffer));
     
/* get the size of A_id_i and the indices */

     buf_count=residual->dim+1;
     buf_count=buf_count*sizeof(int);
     buffer=CREATE(buf_count,char); 
     MPI_Recv(buffer,buf_count,MPI_PACKED,0,20,MPI_COMM_WORLD,&status);
     position=0;
     MPI_Unpack(buffer,buf_count,&position,&indicessize,1,MPI_INT,MPI_COMM_WORLD);
     indices=CREATE(indicessize,int);
     MPI_Unpack(buffer,buf_count,&position,&indices[0],indicessize,MPI_INT,MPI_COMM_WORLD);
     
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
     spbusr_.MSEQNS=indicessize;  
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

     
/* begin the solution */

     s=v_get(A->n);
     sol=v_get(A->m);
     tag=500; 
     
     converged=0;
     buf_count=2*residual->dim;
     sndvector=CREATE(buf_count,double);
     MPI_Barrier(MPI_COMM_WORLD); 

     while ( converged != 1 ) {
       s=multAtransr3(A,residual,s);
       s=permdv(A->n,s,invpv);
       s=gsslv(A->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,s);

/* permute the solution to its original order of columns */

       s=permdv(A->n,s,permv);   
       sol=m_fulvectormlt(A,s,sol);
   
/* send both s and sol to master */

       for (j=0;j<A->n;j++)
         sndvector[j]=s->ve[j];
       for (j=0;j<indicessize;j++) 
         sndvector[j+A->n]=sol->ve[indices[j]-1];      
       buf_count=indicessize+s->dim;
       MPI_Send(&sndvector[0],buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD); 
       free((double *)(s->ve));
       free((double *)(sol->ve));
       free((VECTOR *)(s));
       free((VECTOR *)(sol));  
       
/* probe for a final signal or new residual from master */

       buf_count=residual->dim;
       MPI_Bcast(&residual->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD);  
    
       if (residual->ve[0]!=10000.0) {
         s=v_get(A->n);
         sol=v_get(A->m);
       }
       else {    /* final message */         
         converged=1;
       }
     }
   }   /* end of slave */
   MPI_Finalize();
}    
       
