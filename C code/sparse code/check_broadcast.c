/* last revised 06.05.1999 */
/* benchmark for evaluating the time for broadcast */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "mpi.h"

#define CREATE(Size,Type) (Type *)calloc(Size,sizeof(Type))

double start_time,diff_time;
int myid, numprocs;                         /* to identify the processors */
int dest, tag;                   /* node to send/receive msg, and msg tag */
int buf_count;                       /* the number of elements to be sent */

double bw;

void main(int argc, char *argv[])
{
  static int i,j,k;
  MPI_Status status;
  double *vector;
  
/* initialize MPI and identify yourself */

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  
  buf_count=atoi(argv[1]);
  vector=CREATE(buf_count,double);
  for (i=0;i<buf_count;i++ )
    vector[i]=1.0;
  tag=1;
  start_time=MPI_Wtime();
  for ( i=1;i<=100;i++) {
    MPI_Bcast(&vector[0],buf_count,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  diff_time=MPI_Wtime()-start_time; 
  if ( myid == 0 ) {
    printf("Time elapsed : %lf \n",diff_time);
  }
  MPI_Finalize();
} 
