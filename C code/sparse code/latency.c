/* last revised 08.11.1999 */
/* benchmark for evaluating the latency of a supercomputer */
/* to overhead for sending one byte */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "mpi.h"

double start_time,diff_time;
int myid, numprocs;                         /* to identify the processors */
int dest, tag;                   /* node to send/receive msg, and msg tag */
int buf_count;                       /* the number of elements to be sent */

void main(int argc, char *argv[])
{
  static int j,k;
  char c;
  MPI_Status status;
  
/* initialize MPI and identify yourself */

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  
  tag=1;
  buf_count=0;
  diff_time=0.0;
  c='a';
  start_time=MPI_Wtime();
  for ( j=0;j<10000;j++ ) {
    for ( k=0;k<2;k++ ) {
      if ( myid == k ) {
        dest=(k+1)%2;
/*        start_time=MPI_Wtime(); */
        MPI_Send(&c,0,MPI_CHAR,dest,tag,MPI_COMM_WORLD);
 /*       diff_time+=MPI_Wtime()-start_time; */
      }
      else 
        start_time=MPI_Wtime();
        MPI_Recv(&c,0,MPI_CHAR,MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,&status);
        diff_time+=MPI_Wtime()-start_time;
    }
  }
  if ( myid == 0 ) {
    diff_time=diff_time/10000;
    printf("Latency : %16.14e \n",diff_time);
  }
  MPI_Finalize();
} 
