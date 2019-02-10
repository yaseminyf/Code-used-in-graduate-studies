/* last revised 22.07.1999                                              */
/* jacobi_predict_new_seq.c = routine for solving the LSQ problem iteratively  */
/* and using jacobi updates                                            */
/* this version is to test the predictions with asynchronous communication */
/* and synchronization after predictions. in between two updates of the */
/* residual and x vector a temporary copy of the residual is used in the */
/* computations.                                                        */
/* this is the sequential version of above algo. */
/***********************************************************************/

#include "mpi.h" 
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include "data_struct.h"

VECTOR **D,**Ad,**AD;

/* Ad is used for multiplying A-i with d-i */

MATRIX *hatA,*tildeA,**tildeAptr;
VECTOR **Apptr;
VECTOR *pvec,**pvecptr;
MATRIX *AA;
VECTOR *littles;
int *indices,**indicesptr;

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
int totalapsize;
int **indicesptr,*nofnzptr;
/* variables used for parallel operations */

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
   int nc,t,converged;
   double starttime, endtime, difftime;
   int sind,eind;

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

/* create hatA */

   hatA=CREATE(1,MATRIX);
   hatA->m=residual->dim;     
   hatA->n=ngroup;
   indicesptr=CREATE(ngroup,int *);
   createC1(hatA,ngroup,zeroptr,indicesptr);
printf("created C\n");

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
printf("created Apptr\n");
     
   D=CREATE(ngroup,VECTOR *);
   --D;
   for (i=1;i<=ngroup;i++ ) 
     D[i]=v_get(sub_ptr[i-1]->n);
       
/* allocate the necessary storage to keep the vectors for each group */

   ddiagptr=CREATE(ngroup,double *);
   --ddiagptr;
   rnzptr=CREATE(ngroup,double *);
   --rnzptr;
   nofnzptr=CREATE(ngroup,int);
   xrnzptr=CREATE(ngroup,int *);
   nzsubptr=CREATE(ngroup,int *);
   xnzsubptr=CREATE(ngroup,int *);
   permptr=CREATE(ngroup,int *);
   invpptr=CREATE(ngroup,int *);
   --nofnzptr;
   --xrnzptr;
   --nzsubptr;
   --xnzsubptr;
   --permptr;
   --invpptr;
   tildeAptr=CREATE(ngroup,MATRIX *);
   
/* for each group do the factorization and store the necessary vectors */
    
   for ( l=0; l<ngroup; l++ ) { 
     AA=sub_ptr[l];
     i=l+1; 
     tildeA=CREATE(1,MATRIX);
     tildeA=createtildeA(tildeA,AA,Apptr,ngroup,i,indicesptr);
printf("created tildeA\n");

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
printf("did orcolb\n");
     
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
     
     nofnzptr[l+1]=spbcon_.NOFNZ;
     tildeAptr[l]=tildeA;
     xrnzptr[l+1]=xrnzv;
     nzsubptr[l+1]=nzsubv;
     xnzsubptr[l+1]=xnzsubv;
     permptr[l+1]=permv;
     invpptr[l+1]=invpv;
     
   } /* end of for loop for all groups */       
     
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
printf("did orcolb on hatA\n");
   
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
   counter=0;
   AD=CREATE(ngroup,VECTOR *);
   --AD;
    
printf("before loop...\n");   
   converged=0;   

/* initialize the temporary residual */

   vvec=v_get(residual->dim);
   for (i=0;i<vvec->dim;i++)
     vvec->ve[i]=residual->ve[i];    
   starttime=MPI_Wtime();  
   while (converged==0) {
     for (i=1;i<=ngroup;i++) {
       tildeA=tildeAptr[i-1];
       AA=sub_ptr[i-1];
       ddiag=CREATE(tildeA->n,double);
       lnzv=CREATE(nofnzptr[i],double);
       --ddiag;
       --lnzv;
       j=nofnzptr[i]+tildeA->n;
       isub=CREATE(j,int);
       jsub=CREATE(j,int);
       values=CREATE(j,double);
       --isub;
       --jsub;
       --values;  
       sym_mult3(tildeA,lnzv,ddiag,invpptr[i],xnzsubptr[i],nzsubptr[i],xrnzptr[i],isub,jsub,values);

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
       gsfct(tildeA->n,xrnzptr[i],lnzv,xnzsubptr[i],nzsubptr[i],ddiag,llink,first,temp,flag__); 

/* garbage collection */
      
       ++llink;
       ++temp;
       ++first;
       free((int *)(llink));
       free((double *)(temp));
       free((int *)(first));
       ddiagptr[i]=ddiag;
       rnzptr[i]=lnzv;
     }
     for (t=1;t<=antall;t++) {       
       for (i=1;i<=ngroup;i++) {
         tildeA=tildeAptr[i-1];
         AA=sub_ptr[i-1];
         s=v_get(tildeA->n);
         s=multAtransr3(tildeA,vvec,s);
         s=permdv(tildeA->n,s,invpptr[i]);
         s=gsslv(tildeA->n,xrnzptr[i],rnzptr[i],xnzsubptr[i],nzsubptr[i],ddiagptr[i],s);

/* permute the solution to its original order of columns */

         s=permdv(tildeA->n,s,permptr[i]); 

/* take the littles out of s */

         littles=v_get(D[i]->dim);
         for (j=0;j<D[i]->dim;j++) {
           littles->ve[j]=s->ve[i-1+j];
           D[i]->ve[j]+=s->ve[i-1+j];
         }                
       
/* compute A*D */

         sol=v_get(vvec->dim);
         sol=m_fulvectormlt(AA,littles,sol);
         AD[i]=sol; 
         free((double *)(s->ve));
         free((VECTOR *)(s));        
       } 
       for (i=1;i<=ngroup;i++) {
         eind=Ad[i]->dim;
         for (j=0;j<eind;j++) {
           verdi=AD[i]->ve[indicesptr[i-1][j]-1];
           if (t!=antall) 
             vvec->ve[indicesptr[i-1][j]-1]+=verdi;
           Ad[i]->ve[j]+=verdi;
         }
         free((double *)(AD[i]->ve));
         free((VECTOR *)(AD[i]));
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
         Apptr[i-1]->ve[count]=verdi;
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
     }
     else {
       for (i=1;i<=ngroup;i++) {
         tildeA=tildeAptr[i-1];
         AA=sub_ptr[i-1];
         putvaluesintildeA3(tildeA,Apptr,AA,ngroup,i);  
       }  
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
   printf("\n The error is: %16.14e \n",err);
#else   
   printf("\n Two norm of residual is: %16.14e \n",r_norm); 
#endif
   printf("\n Elapsed time: %16.14e seconds \n",difftime);    
   
   MPI_Finalize(); 
}    
