/* last revised 18.05.1999                                             */
/* multcheck.c = benchmark for checking the time for execution of      */
/* one double multiplication                                           */
/***********************************************************************/

#include <unistd.h>
#include <limits.h>
#include <time.h>
#include "data_struct.h"   
                                     
/**************************************************************************/
int main(argc,argv)
int argc;
char *argv[];
{  

/* local variables */
   
   clock_t time1,time2;
   double dtime;
   int i;
   double value1,value2,value3;
   
   value1=rand()/0.95012928514718;
   value2=rand()/0.23113851357429;
   
   time1=clock();
   for (i=1;i<=1000000;i++)   
      value3=value1*value2;
   time2=clock();

   dtime=(double) (time2)- (double) (time1);
   dtime=dtime/CLOCKS_PER_SEC;
   printf("\n value1 %lf, value2 %lf",value1,value2);
   printf("\n time: %16.10e seconds \n",dtime);    
   exit(0); 
}    
       
