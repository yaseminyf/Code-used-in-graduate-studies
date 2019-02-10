/* last revised 18.05.1999                                             */
/* multcheck.c = benchmark for checking the time for execution of      */
/* one double multiplication                                           */
/***********************************************************************/

#include "mpi.h"
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

int myid, numprocs;
                                      
/**************************************************************************/
int main(argc,argv)
int argc;
char *argv[];
{  

/* local variables */
   
   int i;
   double value1,value2,value3;
   double starttime, endtime, difftime;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

   if (myid==0) {
   
     value1=rand()/0.95012928514718;
     value2=rand()/0.23113851357429;
   
     starttime=MPI_Wtime();
     for (i=1;i<=1000000;i++)   
        value3=value1*value2;
     endtime=MPI_Wtime();     
     difftime=endtime-starttime;  

     printf("\n time: %16.14e seconds \n",difftime);    
   }
   MPI_Finalize();
}    
       
