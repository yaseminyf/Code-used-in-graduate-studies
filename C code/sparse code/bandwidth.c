/* last revised 18.05.1999 */
/* benchmark for evaluating the bandwidth of a supercomputer */

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
  char *vector;
  
/* initialize MPI and identify yourself */

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  
  buf_count=atoi(argv[1]);
  vector=CREATE(buf_count,char);
  if ( myid==0) {
    for (i=0;i<buf_count;i++ )
      vector[i]='a';
  }
  tag=1;
  start_time=MPI_Wtime();
  for ( j=0;j<10000;j++ ) {
    for ( k=0;k<2;k++ ) {
      if ( myid == k ) {
        dest=(k+1)%2;
        MPI_Send(&vector[0],buf_count,MPI_CHAR,dest,tag,MPI_COMM_WORLD); 
      }
      else 
        MPI_Recv(&vector[0],buf_count,MPI_CHAR,MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,&status);
    }
  }
  diff_time=MPI_Wtime()-start_time; 
  if ( myid == 0 ) {
    bw=(2*10000*buf_count/diff_time)/1000000;
    printf("Bandwidth : %lf Mbytes/s\n",bw);
    printf("Time elapsed : %lf \n",diff_time);
    printf("Buffer size : %d \n",buf_count);
  }
  MPI_Finalize();
} 
