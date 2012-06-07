#include <time.h>
#include <stdio.h>
#include "integ_hyp.h"
#include "test_utils.h"
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

typedef struct sin_par{
  int i;
}sin_par;

double my_sin(double, void *);

int main(int argc, char *argv[]) {
  
  double *x,*x1,*x2,*x3,*x4;
  double *w,*w1,*w2,*w3,*w4;
  double y, yscale, integ1, integ2, integ3, integ4,ratio;
  sin_par p = {0};
  int i,npts,low_lim_n,up_lim_n,exact_n;
  FILE *fp;
  gsl_function f;
  
  y = PI/4.0;
  npts = atoi(argv[1]);
  low_lim_n = atoi(argv[2]);
  up_lim_n = atoi(argv[3]);
  exact_n = atoi(argv[4]);
  ratio = atof(argv[5]);
  w = (double *)malloc((npts+1)*sizeof(double));
  x = (double *)malloc((npts+1)*sizeof(double));
  x1 = (double *)malloc((npts+1)*sizeof(double));
  w1 = (double *)malloc((npts+1)*sizeof(double));
  x2 = (double *)malloc((npts+1)*sizeof(double));
  w2 = (double *)malloc((npts+1)*sizeof(double));
  x3 = (double *)malloc((npts+1)*sizeof(double));
  w3 = (double *)malloc((npts+1)*sizeof(double));
  x4 = (double *)malloc((npts+1)*sizeof(double));
  w4 = (double *)malloc((npts+1)*sizeof(double));
  gauss_legendre(-1.0,1.0,x,w,npts);

  /*scale y*/
  yscale = 0.0;
  
  quad_log_singularity_half(yscale,npts,x,w,x3,w3);
  quad_x_singularity(yscale,npts,x,w,x2,w2);
  quad_x2_singularity(yscale,npts,x,w,x1,w1);

  f.function = &my_sin;
  f.params = &p;

  /* integ1 = integral_gaussleg(f, -0.0, PI/2.0, npts, x, w); */
  /* integ2 = integral_gauss_log_sing(f,-0.0,PI/2.0,y,npts,x,w,x3,w3); */
  /* integ3 = integral_gauss_x_sing(f,0.0,PI/2.0,y,npts,x,w,x2,w2); */
  /* integ4 = integral_gauss_x2_sing(f,0.0,PI/2.0,y,npts,x,w,x1,w1); */
  /* integ1 = integral_generalized_sing(f,0.0,PI/2.0,y,npts,(2*npts),x,w,x4,w4); */
  
  integ1 = integral_gaussleg(f,-1.0,1.0, npts, x, w);
  integ2 = integral_gauss_log_sing(f,-1.0,1.0,0.0,npts,x,w,x3,w3);
  integ3 = integral_gauss_x_sing(f,-1.0,1.0,0.0,npts,x,w,x2,w2);
  integ4 = integral_gauss_x2_sing(f,-1.0,1.0,0.0,npts,x,w,x1,w1);
  integ1 = integral_generalized_sing(f,-1.0,1.0,0.0,npts,(npts/4)+1,x,w,x4,w4); 
  
  fp = fopen("testlog","a");
  fprintf(fp,"weights w1n: \n");
  for(i=1;i<=npts;i++) {
    fprintf(fp,"%.18e \n",w3[i]);
    
  }
  fprintf(fp,"weights w2n: \n");
  for(i=1;i<=npts;i++) {
    fprintf(fp,"%.18e \n",w2[i]);
    
  }
  fprintf(fp,"weights w3n: \n");
  for(i=1;i<=npts;i++) {
    fprintf(fp,"%.18e \n",w1[i]);
    
  }

  fprintf(fp,"integrand from -1 to 1: %.18e \n",integ1);
  fprintf(fp,"log(x)*integrand from -1, 1: %.18e \n",integ2);
  fprintf(fp,"integrand/(x) from -1, 1: %.18e \n",integ3);
  fprintf(fp,"integrand/((x)^2) from -1,1: %.18e \n",integ4);
  fclose(fp);
  fp = fopen("testlog2","a");
  integral_generalized_test(f,-1.0,1.0,0.0,low_lim_n,up_lim_n,exact_n,ratio,fp);
  fclose(fp);
  /* fp = fopen("testlog3","a"); */
  /* integral_generalized_benchmark(f,-1.0,1.0,0.0,1,100,10,fp); */
  /* fclose(fp); */

  /* fp = fopen("testlog3","a"); */

  /* for(j = 0;j<100;j++){ */
  /*   t1 = (double)clock()/CLOCKS_PER_SEC; */
  /*   for(i = 0; i< 10; i++) { */

  /*     integral_generalized_sing(f,-1.0,1.0,0.0,npts,(npts/4)+1,x,w,x4,w4); */
  /*   } */

  /*   t2 = (double)clock()/CLOCKS_PER_SEC; */
  /*   time_elapsed = t2 - t1; */

  /*   time_elapsed = time_elapsed*1000.00; */
	
  /*   avg_time = time_elapsed/10.00; */
  /*   fprintf(fp,"%d  %lf\n",j, */
  /* } */



  
  free(x);
  free(x1);
  free(x2);
  free(x3);
  free(x4);
  free(w);
  free(w1);
  free(w2);
  free(w3);
  free(w4);
  return 0;
}

double my_sin(double x, void *params){

  /* if(fabs(x)>EPS){ */
  /*   return gsl_sf_bessel_Kn(2,fabs(x)); */
  /* } else { */
  /*   return 0.0; */
  /* } */

  return (cos(x)+sin(x+(PI/3.0)))/x;
}
  
  
  
