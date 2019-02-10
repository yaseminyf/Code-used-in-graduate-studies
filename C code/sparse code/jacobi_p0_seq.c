/* last revised 26.07.1999                                             */
/* jacobi_p0_seq.c = routine for solving the LSQ problem iteratively       */
/* and using jacobi updates and subspace corrections                   */
/* This is the sequential algorithm                                    */
/***********************************************************************/
#include "mpi.h"
#include <unistd.h>
#include <limits.h>
#include <float.h> 
#include <time.h>
#include "data_struct.h"

VECTOR **D;
MATRIX *C;
VECTOR **Ad;
int **indicesptr;

/* #define XERROR  */
/* #define LINESEARCH
*/ 
#ifdef XERROR
  VECTOR *xnormv;
#endif

MATRIX *AA;  
int myid, numprocs;                         /* to identify the processors */
                                   
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
   int nc,converged;
   int sind,eind;
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
  
   x=v_get(A->n); /* initial x vector, all zeros */
   residual=v_get(rhs->dim);
   for ( i=0;i<rhs->dim;i++ )
     residual->ve[i]=-1.0*-rhs->ve[i];

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
   printf("\n initial norm: %16.14e\n",r_zero_norm);

/* allocate the necessary storage to keep the vectors for each group */

   rnzptr=CREATE(ngroup,double *);
   ddiagptr=CREATE(ngroup,double *);
   xrnzptr=CREATE(ngroup,int *);
   nzsubptr=CREATE(ngroup,int *);
   xnzsubptr=CREATE(ngroup,int *);
   permptr=CREATE(ngroup,int *);
   invpptr=CREATE(ngroup,int *);
   --rnzptr;
   --ddiagptr;
   --xrnzptr;
   --nzsubptr;
   --xnzsubptr;
   --permptr;
   --invpptr;
   
/* initialize for SPARSPAK */  

   spbusr_.MAXSB=12000;
   TOL=0.0000000001;
   TYPTOL=1;
   spksys_.MAXINT=SHRT_MAX;
   spbusr_.IERRB=0;
   spbusr_.MSGLVB=1; 
   
/* for each group do the factorization and store the necessary vectors */
    
   for ( l=0; l<ngroup; l++ ) {
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
     AA=sub_ptr[l];
     spbusr_.MCOLS=AA->n;
     spbusr_.MSEQNS=AA->m-zeroptr[l][0];
     spbcon_.NCOLS=spbusr_.MCOLS;  
        
     T=CREATE(spbusr_.MAXSB,int);
     for (i=0;i<spbusr_.MAXSB;i++) {
       T[i]=-(i+1);
     }
     --T;

/* read in the matrix into SPARSPAK */

     TYPE=1;   
     rnum=1;
    
     for ( j=1;j<=AA->m;j++ ) {
       if ( AA->rowptr[j] != 0 ) { 
         count=1;
         while ( AA->rowptr[j+count] == 0 )
           count++;
         NSUBS=AA->rowptr[j+count]-AA->rowptr[j];
         spbcon_.NOFNZ+=NSUBS;
         ROWNUM=rnum;
         rnum++;
         SUBS=CREATE(NSUBS,int);
         for ( k=0;k<NSUBS;k++ ) {
           *(SUBS+k)=AA->colind[AA->rowptr[j]+k];
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
     
     sym_mult3(AA,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values);

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
   
     ddiagptr[l+1]=ddiag;
     xrnzptr[l+1]=xrnzv;
     nzsubptr[l+1]=nzsubv;
     xnzsubptr[l+1]=xnzsubv;
     permptr[l+1]=permv;
     invpptr[l+1]=invpv;
     rnzptr[l+1]=lnzv;
     
   } /* end of for loop for all groups */  

/* do the symmetric factorization of the C matrix */

   C=CREATE(1,MATRIX);
   C->m=residual->dim;     
   C->n=ngroup;
   indicesptr=CREATE(ngroup,int *);     
   createC1(C,ngroup,zeroptr,indicesptr);

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

   converged=0;
   counter=0;

   D=CREATE(ngroup,VECTOR *);
   --D;

   Ad=CREATE(ngroup,VECTOR *);
   --Ad;

printf("\n before loop..\n");       

   starttime=MPI_Wtime();  
   while (converged==0) {
     for (i=1;i<=ngroup;i++) {
       AA=sub_ptr[i-1];
       s=v_get(AA->n);
       sol=v_get(AA->m);
       s=multAtransr3(AA,residual,s);
       s=permdv(AA->n,s,invpptr[i]);
       s=gsslv(AA->n,xrnzptr[i],rnzptr[i],xnzsubptr[i],nzsubptr[i],ddiagptr[i],s);
       
/* permute the solution to its original order of columns */

       s=permdv(AA->n,s,permptr[i]);      
       sol=m_fulvectormlt(AA,s,sol);
             
       Ad[i]=sol;
       D[i]=s;    
     }
#ifdef LINESEARCH
       s=v_get(ngroup);
       for (i=1;i<=ngroup;i++) { 
         check=0.0;
         sind=C->Jind[i];
         eind=C->Jind[i+1]-1;
         dotnorm=0.0;
         for (j=sind;j<=eind;j++) {
           verdi=Ad[i]->ve[C->Iind[j]-1];
           check+=verdi*residual->ve[C->Iind[j]-1];
           dotnorm+=verdi*verdi;
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
     createC22(C,Ad,ngroup);
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
     
/* update x vector */

     k=0;
     for (i=1;i<=ngroup;i++) {
       for (j=1;j<=D[i]->dim;j++) {
         x->ve[k]+=D[i]->ve[j-1]*s->ve[i-1];
         k++;
       }
     }
     
/* update the residual vector */

     for (i=1;i<=ngroup;i++) {
       sind=C->Jind[i];
       eind=C->Jind[i+1]-1;
       for (j=sind;j<=eind;j++) 
         residual->ve[C->Iind[j]-1]+=Ad[i]->ve[C->Iind[j]-1]*s->ve[i-1];
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
/* printf("\n residual norm: %16.14e\n",r_norm);         */
     
     err=r_norm/r_zero_norm;
#endif                 

     counter++;
     if ((err < epsilon)||(r_norm==0.0)||(counter>15000)) {
       converged=1;
     }
     else {

/* garbage collection */

       free((double *)(s->ve));
       free((VECTOR *)(s));
#ifndef LINESEARCH
       ++ddiag;
       ++lnzv;
       free((double *)(ddiag));
       free((double *)(lnzv));
#endif           
       for (i=1;i<=ngroup;i++) {
         free((double *)(D[i]->ve));
         free((double *)(Ad[i]->ve));
       }
     }     
   }
   endtime=MPI_Wtime();     
   difftime=endtime-starttime;     
   
   fp=fopen("results_seq_p0","a");
   fprintf(fp,"\n Converged after %d loops.",counter);
#ifdef XERROR
   fprintf(fp,"\n The error is: %16.14e \n",err   );
#else   
   fprintf(fp,"\n Two norm of residual is: %16.14e \n",r_norm); 
#endif    
   fprintf(fp,"\n Elapsed time: %16.14e seconds \n",difftime);  
   fclose(fp);  
   MPI_Finalize();
}    
       
