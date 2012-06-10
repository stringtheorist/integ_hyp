#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "integ_hyp.h"
#include "test_utils.h"

/*this function generates a file containing error convergence rates for the generalized quadrature method for an arbitrary
  funciton specified by f. */

int integral_generalized_test(gsl_function f, double a, double b, double y,int low_lim_n, int up_lim_n, int exact_n, double ratio, FILE *fp)
{
  double exact, integ2, error;
  int err_msg;
  int i;
  double *x_gauss,*w_gauss,*x,*w;

  

  x_gauss = (double *)malloc((exact_n+1)*sizeof(double));
  w_gauss = (double *)malloc((exact_n+1)*sizeof(double));
  x = (double *)malloc((exact_n+1)*sizeof(double));
  w = (double *)malloc((exact_n+1)*sizeof(double));

  gauss_legendre(-1.0,1.0,x_gauss,w_gauss,exact_n);

  exact = integral_generalized_sing(f,a,b,y,exact_n,(exact_n/4)+1,x_gauss,w_gauss,x,w);

  free(x_gauss);
  free(w_gauss);
  free(x);
  free(w);
  

  err_msg = 0;

  for(i=low_lim_n;i<up_lim_n;i++)
    {
      x_gauss = (double *)malloc((i+1)*sizeof(double));
      w_gauss = (double *)malloc((i+1)*sizeof(double));
      x = (double *)malloc((i+1)*sizeof(double));
      w = (double *)malloc((i+1)*sizeof(double));

      /*generate gauss-legendre weights*/
      
      gauss_legendre(-1.0,1.0,x_gauss,w_gauss,i);
      
      /*perform generalized integration*/

      integ2 = integral_generalized_sing(f,a,b,y,i,((int)(i/ratio)) + 1,x_gauss,w_gauss,x,w);
      if(!isfinite(integ2)) continue;
      /*fprintf(stdout,"%.18e \n",integ2);*/
      error = fabs((integ2 - exact)/exact);
      
      fprintf(fp,"%d, %.18e \n",i,error);

      
      
      free(x_gauss);
      free(w_gauss);
      free(x);
      free(w);
    }
  
  return err_msg;
}



int integral_generalized_benchmark(gsl_function f, double a, double b, double y, int low_lim_n, int up_lim_n, int reps, FILE *fp)
{
  double t1, t2, time_elapsed, avg_time;
  double *x_gauss, *w_gauss,*x,*w;
  int i,j;
  for(i=low_lim_n;i<up_lim_n;i++)
    {
      x_gauss = (double *)malloc((i+1)*sizeof(double));
      w_gauss = (double *)malloc((i+1)*sizeof(double));
      x = (double *)malloc((i+1)*sizeof(double));
      w = (double *)malloc((i+1)*sizeof(double));

      /*generate gauss-legendre weights*/
      
      gauss_legendre(-1.0,1.0,x_gauss,w_gauss,i);
      
      /*perform generalized integration*/

      t1 = (double)clock()/CLOCKS_PER_SEC;
      for(j=0;j<reps;j++){
      integral_generalized_sing(f,a,b,y,i,((int)(i/4)) + 1,x_gauss,w_gauss,x,w);
      }
      t2 = (double)clock()/CLOCKS_PER_SEC;

      time_elapsed = (t2-t1)*1000.0;
      avg_time = time_elapsed/(double)reps;

     
      


      fprintf(fp,"%d, %lf\n",i,avg_time);

      
      
      free(x_gauss);
      free(w_gauss);
      free(x);
      free(w);
    }
  
  return 0;
}


      

  

  
