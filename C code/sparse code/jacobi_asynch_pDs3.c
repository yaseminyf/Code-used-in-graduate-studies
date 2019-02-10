/* last revised 25.01.1999                                             */
/* jacobi_asynch_pDs1.c = routine for solving the LSQ problem iteratively  */
/* and using jacobi updates                                            */
/* this version computes Jacobi p=DS                                   */
/* iterations on a parallel computer asynchronously                    */
/* p is gradually updated after receiving results from the slaves      */
/* planar search is done after each update and each slave receives new */
/* A_ip_i at the end of sending an update to the master                */
/***********************************************************************/

/* tag = 100 : used for sending/receiving computed data from slaves    */

#include "mpi.h" 
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include "data_struct.h"

#define PONTUS 

#define CHOL  

/* #define PERM    */
/* #define MMDP    */
/* #define COLP      */
/* #define NEWPERM    */
  
VECTOR **D,**Ad,**AD;

/* Ad is used for multiplying A-i with d-i */
/* AD is used for computing the nonzeros of hatA */

MATRIX *hatA,*tildeA;
VECTOR **Apptr;
VECTOR *pvec,**pvecptr,*pd;
MATRIX *AA;
double *deltav;
VECTOR *littles;

/* #define XERROR */   /* to check convergence with error on x */

/* #define FMSTART */ /* starting values of p are F_M values */

#ifdef FMSTART
  VECTOR *evec;
#endif
  
#ifdef XERROR
  VECTOR *xnormv;
#endif
  
#ifdef PERM
  int *p,*ivp;
  double *tmp;
#endif
 
#define DIFFGROUPSIZE       

/* variables used for parallel operations */

int myid, numprocs;                         /* to identify the processors */
int dest, tag;                   /* node to send/receive msg, and msg tag */
int buf_count;                       /* the number of elements to be sent */
char *sndvector,**sndptr;
int position;                                                                          
/**************************************************************************/
int main(argc,argv)
int argc;
char *argv[];
{  

/* local variables */
   
   static int i,j,k,l,m;
   int rnum, count, counter;
   char choice1,c1;
   double verdi;
   double number,num,*ip;
   int nc,tt,t,converged,planeflag;
   double starttime, endtime, difftime;
   
   MPI_Request request,*req;
   MPI_Status status;

#ifdef DIFFGROUPSIZE   
   int choice2;
#endif

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
     sscanf(argv[1],"%c",&choice1);  /* get the test matrix number */

#ifdef DIFFGROUPSIZE  
     sscanf(argv[2],"%d",&choice2);  /* get the number of groups */
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
       xexact->ve=tmp;
#endif
      
       rhs=m_fulvectormlt(A,xexact,rhs); 

/* garbage collection */

#ifndef XERROR
       free((double *)(xexact->ve));
       free((VECTOR *)(xexact));      
#endif
      
     }
     else if ( (choice1 >= '6') && (choice1<= '9') ) {
       epsilon=1.0e-03;
       fp1=fopen("data","w");
       readnewmatrix(fp,fp1);
     }

     x=v_get(A->n);  /* starting vector x, all zeros */

     residual=v_get(rhs->dim);  

     for ( i=0;i<rhs->dim;i++ )
       residual->ve[i]=-1.0*rhs->ve[i];

/* garbage collection */

     free((double *)(rhs->ve));
     free((VECTOR *)(rhs)); 

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

       free((int *)(A->Jind));
       free((int *)(A->Iind));
       free((int *)(A->colind));
       free((int *)(A->rowptr));
       free((double *)(A->Valuerow));
       free((double *)(A->Valuecol));
       free((MATRIX *)(A)); 
     
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
       buf_count=1;    
       MPI_Send(&AA->n,buf_count,MPI_INT,dest,tag,MPI_COMM_WORLD);
      
       tag=3;
       buf_count=1;    
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
       buf_count=AA->nonzeros+1;
       MPI_Send(&AA->Valuecol[0],buf_count,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);

       tag=9;
       buf_count=AA->nonzeros+1;
       MPI_Send(&AA->Valuerow[0],buf_count,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
     }
     
/* broadcast number of groups */
     
     buf_count=1;
     MPI_Bcast(&ngroup,buf_count,MPI_INT,0,MPI_COMM_WORLD);
     
/* broadcast the residual */     

     buf_count=residual->dim;
     MPI_Bcast(&residual->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD);
     
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
printf("\n before broadcast..\n");
     buf_count=200000;
     sndvector=CREATE(buf_count,char);
     position=0;
     j=residual->dim;
     for (i=0;i<ngroup;i++) {
       MPI_Pack(&Apptr[i]->ve[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
     } 
     MPI_Bcast(sndvector,position,MPI_PACKED,0,MPI_COMM_WORLD);
    
     free((char *)(sndvector));
     
/* do the symmetric factorization of hatA */

     hatA=CREATE(1,MATRIX);
     hatA->m=residual->dim;     
     hatA->n=ngroup;
     createC1(sub_ptr,hatA,ngroup);
printf("\n created hata..\n");

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
    
     spbmap_.LLFREE=spbusr_.MAXSB;

     spbcon_.NEDGES=0;
     spbcon_.NZMAX=0;
     spbcon_.NZEQNS=0;
     spbcon_.NZCONS=0;
     spbusr_.MCOLS=hatA->n;
     spbusr_.MSEQNS=hatA->m;
    
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
     for ( j=1;j<=hatA->m;j++ ) {
       if ( hatA->rowptr[j] != 0 ) { 
         count=1;
         while ( hatA->rowptr[j+count] == 0 )
           count++;
         NSUBS=hatA->rowptr[j+count]-hatA->rowptr[j];
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
     }  
      
     orcolb(T);
printf("\n after orcolb..\n");

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
   
     D=CREATE(ngroup,VECTOR *);
     --D;
     for (i=1;i<=ngroup;i++)
       D[i]=v_get(pvecptr[i-1]->dim);
     
     Ad=CREATE(ngroup,VECTOR *);
     --Ad;
     for (i=1;i<=ngroup;i++)
       Ad[i]=v_get(residual->dim);

     AD=CREATE(ngroup,VECTOR *);
     --AD;
     
     sndptr=CREATE(ngroup,char *);
     for (i=0;i<ngroup;i++) { 
       sndvector=CREATE(200000,char);
       sndptr[i]=sndvector;
     } 
     req=CREATE(ngroup+1,MPI_Request);

printf("\n before loop...\n");   
     
     converged=0; 
     MPI_Barrier(MPI_COMM_WORLD);   
     starttime=MPI_Wtime();  

     while (converged==0) {
         
/* do for each slave */
/* get the solution vectors from all processors */
/* divide the gotten vector into 2 and assign to D and Ad */

       buf_count=200000;
       for (i=1;i<=ngroup;i++) {
         position=0;
         MPI_Irecv(&sndptr[i-1][0],buf_count,MPI_PACKED,i,100,MPI_COMM_WORLD,&req[i]);
       }
       for (i=1;i<=ngroup;i++) {
         MPI_Waitany(ngroup+1,req,&j,&status);
         dest=status.MPI_SOURCE;
         buf_count=200000; 
         position=0;
         MPI_Unpack(&sndptr[dest-1][0],buf_count,&position,&verdi,1,MPI_DOUBLE,MPI_COMM_WORLD);
         sol=v_get(residual->dim);
         MPI_Unpack(&sndptr[dest-1][0],buf_count,&position,&sol->ve[0],residual->dim,MPI_DOUBLE,MPI_COMM_WORLD);
         s=v_get((int)(verdi));
         MPI_Unpack(&sndptr[dest-1][0],buf_count,&position,&s->ve[0],s->dim,MPI_DOUBLE,MPI_COMM_WORLD);  
         D[dest]=s;
         Ad[dest]=sol;     

/* do planar search */
 
/* compute hatA^T*hatA and initialize lnz and diag */

         ddiag=CREATE(spbcon_.NCOLS,double);
         lnzv=CREATE(spbcon_.NOFNZ,double);
         --ddiag;
         --lnzv;
     
/* get the necessary storage for hatA^ThatA */

         k=spbcon_.NOFNZ+spbcon_.NCOLS;
         isub=CREATE(k,int);
         jsub=CREATE(k,int);
         values=CREATE(k,double);
         --isub;
         --jsub;
         --values;

/* compute D matrix first */

         pd=v_get(x->dim);
         for (m=1;m<=ngroup;m++) {
           count=0;
           if (D[m]->ve[count]!=0.0) {
             k=0;
             for (j=1;j<=ngroup;j++) {
               if (m!=j) {
                 for (l=1;l<=pvecptr[j-1]->dim;l++) {
                   pd->ve[k]+=pvecptr[j-1]->ve[l-1]*D[m]->ve[count];
                   k++;
                 }
                 count++;
               }
               else {
                 t=D[m]->dim-ngroup+1;
                 for (l=1;l<=t;l++) {
                   pd->ve[k]+=D[m]->ve[count];
                   count++;
                   k++;
                 }
               }
             }
           }
         }

/* compute the delta values for each group */

         deltav=CREATE(ngroup,double);
         --deltav;
         for (k=1;k<=ngroup;k++) {
           deltav[k]=0.0;
           for (j=1;j<=ngroup;j++) {
             if (k!=j) {
               if (j>k) {
                 deltav[k]+=D[j]->ve[k-1];
               }
               if (j<k) {
                 deltav[k]+=D[j]->ve[j-1+nc+k-j-1];
               }
             }
           }
         }
     
/* compute the columns of A*D=hatA matrix */
   
         for (k=1;k<=ngroup;k++) {
           sol=v_get(residual->dim);
           for (j=0;j<Ad[k]->dim;j++) {
             sol->ve[j]+=Ad[k]->ve[j]+Apptr[k-1]->ve[j]*deltav[k];
           }
           AD[k]=sol;
         }     
                   
/* put in the values of hatA */

         createC22(hatA,AD,ngroup);  
         sym_mult(hatA,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values);

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

         s=v_get(hatA->n);
         s=multAtransr3(hatA,residual,s);
         s=permdv(hatA->n,s,invpv);
         s=gsslv(hatA->n,xrnzv,lnzv,xnzsubv,nzsubv,ddiag,s);   
         s=permdv(hatA->n,s,permv);  

/* update x vector and the pvecptr at the same time */
            
         for (j=0;j<ngroup;j++) {
           free((double *)(pvecptr[j]->ve));
           free((VECTOR *)(pvecptr[j])); 
         }

         k=0;
         for (l=1;l<=ngroup;l++) {
           t=D[l]->dim-ngroup+1;
           pvec=v_get(t); 
           for (j=1;j<=t;j++) {
             verdi=pd->ve[k]*s->ve[l-1];
             x->ve[k]+=verdi;
             pvec->ve[j-1]=verdi;
             k++;       
           }
           pvecptr[l-1]=pvec;
         }
     
/* update the residual vector and Apptr at the same time */


         if ( choice1 >= '6' )
           rhs=v_get(residual->dim);
         
         for (l=1;l<=ngroup;l++) {                              
           for ( j=0;j<residual->dim;j++ ) {
             verdi=s->ve[l-1]*AD[l]->ve[j];
             residual->ve[j]+=verdi;
             if ( choice1 >= '6' )
               rhs->ve[j]+=verdi;
             Apptr[l-1]->ve[j]=verdi; 
           }                                  
         }
          
/* garbage collection */

         free((double *)(s->ve));
         free((VECTOR *)(s));
         ++ddiag;
         ++lnzv;
         free((double *)(ddiag));
         free((double *)(lnzv));
         for (j=1;j<=ngroup;j++) {
           free((double *)(AD[j]->ve));
         }
         free((double *)(pd->ve)); 
         ++deltav;
         free((double *)(deltav));
         free((double *)(hatA->Valuerow));
         free((double *)(hatA->Valuecol));          
/* end of planar search */
               
/* find the error */

#ifdef XERROR
        for (j=0;j<xnormv->dim;j++)
          xnormv->ve[j]=x->ve[j]-xexact->ve[j];
        err=norm2(xnormv);
        r_norm=norm2(residual);

printf("\n error: %16.14e\n",err);  
printf("\n residual norm: %16.14e\n",r_norm);        

#else
        if ( choice1 <= '5' ) {
          r_norm=norm2(residual);

printf("\n residual norm: %16.14e\n",r_norm);       
       
        }
        else {
          rvec2=multAtransr2(A,rhs,rvec2);
          free((double *)(rhs->ve));
          free((VECTOR *)(rhs));
          r_norm=norm2(rvec2);

printf("\n residual norm: %16.14e\n",r_norm);       
         
        }
        err=r_norm/r_zero_norm;
#endif     

/* if not converged then broadcast residual */
              
        if ((err < epsilon)||(r_norm==0.0)||(counter>=10000)) { /* final */
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
          buf_count=200000;
          sndvector=CREATE(buf_count,char);
          position=0;
          j=residual->dim;
          for (k=0;k<ngroup;k++) {
            MPI_Pack(&Apptr[k]->ve[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
          } 
          MPI_Pack(&residual->ve[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
          tag=2000;
          MPI_Isend(sndvector,position,MPI_PACKED,dest,tag,MPI_COMM_WORLD,&request);        
        }
        counter++;
      }                   /* end of for loop */     
    }                     /* end of while loop */

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
    buf_count=1;
    MPI_Recv(&A->m,buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

    tag=2;
    buf_count=1;
    MPI_Recv(&A->n,buf_count,MPI_INT,0,tag,MPI_COMM_WORLD,&status);

    tag=3;
    buf_count=1;
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
    buf_count=A->nonzeros+1;
    A->Valuecol=CREATE(buf_count,double);
    MPI_Recv(&A->Valuecol[0],buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);

    tag=9;
    buf_count=A->nonzeros+1;
    A->Valuerow=CREATE(buf_count,double);
    MPI_Recv(&A->Valuerow[0],buf_count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);

    buf_count=1;
    MPI_Bcast(&ngroup,buf_count,MPI_INT,0,MPI_COMM_WORLD);

    residual=v_get(A->m);
    buf_count=residual->dim;
    MPI_Bcast(&residual->ve[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD); 
         
    Apptr=CREATE(ngroup,VECTOR *);

/* broadcast the Apptr */

    buf_count=200000;
    sndvector=CREATE(buf_count,char);
    position=0;
    j=A->m;
    MPI_Bcast(sndvector,buf_count,MPI_PACKED,0,MPI_COMM_WORLD);
    for (i=0;i<ngroup;i++) {
      Apptr[i]=v_get(j);
      MPI_Unpack(sndvector,buf_count,&position,&Apptr[i]->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
    }     
    
    free((char *)(sndvector));
 
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
     
/* begin iteration */  

    counter=0;
    verdi=(double)(tildeA->n);
    planeflag=1;
    sndvector=CREATE(200000,char);
          
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
      sym_mult(tildeA,lnzv,ddiag,invpv,xnzsubv,nzsubv,xrnzv,isub,jsub,values);

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

/* take the littles out of s */

      littles=v_get(A->n);
      for (j=0;j<A->n;j++) {
        littles->ve[j]=s->ve[myid-1+j];
      }                
       
/* compute A*d */

      sol=v_get(A->m);
      sol=m_fulvectormlt(A,littles,sol);    
       
/* send s and sol with verdi to master */

      tag=100;
      buf_count=200000;
      position=0; 
      MPI_Pack(&verdi,1,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
      j=sol->dim;
      MPI_Pack(&sol->ve[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
      j=s->dim;
      MPI_Pack(&s->ve[0],j,MPI_DOUBLE,sndvector,buf_count,&position,MPI_COMM_WORLD);
      MPI_Isend(sndvector,position,MPI_PACKED,0,tag,MPI_COMM_WORLD,&request);        
        
      MPI_Wait(&request,&status);
      MPI_Probe(0,MPI_ANY_TAG,MPI_COMM_WORLD,&status); 
      if (status.MPI_TAG==2000) { /* new Aipi */
        MPI_Recv(&sndvector[0],buf_count,MPI_PACKED,0,2000,MPI_COMM_WORLD,&status);
        position=0;
        j=A->m;
        for (i=0;i<ngroup;i++) {
          MPI_Unpack(&sndvector[0],buf_count,&position,&Apptr[i]->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);
        } 
        MPI_Unpack(&sndvector[0],buf_count,&position,&residual->ve[0],j,MPI_DOUBLE,MPI_COMM_WORLD);    
        i=myid;
        putvaluesintildeA3(tildeA,Apptr,A,ngroup,i);         
        ++ddiag;
        ++lnzv;
        free((double *)(ddiag));
        free((double *)(lnzv));        
        free((double *)(s->ve));
        free((double *)(sol->ve));
        free((VECTOR *)(s));
        free((VECTOR *)(sol)); 
        free((double *)(littles->ve));
        free((VECTOR *)(littles));  
      }
      else { /* final */ 
        MPI_Recv(&tag,1,MPI_INT,0,1000,MPI_COMM_WORLD,&status);
        counter=1;                  
      }
    }
  }
  MPI_Finalize(); 
}    
