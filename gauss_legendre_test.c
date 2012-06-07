#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_legendre.h>
#include "integ_hyp.h"
typedef struct test_sin_params{
  int i;
}test_sin_params;



double test_sin(double x, void *params);


int main(int argc, char *argv[])
{
  int npts, i, j;
  double x_left, x_right, y;
  double *x;
  double *x1;
  double *x2;
  double *x3;
  double *x_test;
  double *x_test2;
  double *w;
  double *w1;
  double *w2;
  double *w3;
  double *w_test;
  double *w_test2;
  double z1,z2,z3;
  FILE *fp;
  gsl_function f;
  test_sin_params p = {0};


  
  x_left = -1.0;
  x_right = 1.0;
  y = -0.986284;

  npts = atoi(argv[1]);
  /*allocate memory for arrays*/
  x = (double *)malloc((npts+1)*sizeof(double));
  x1 = (double *)malloc((npts+1)*sizeof(double));
  x2 = (double *)malloc((npts+1)*sizeof(double));
  x3 = (double *)malloc((npts+1)*sizeof(double));
  w = (double *)malloc((npts+1)*sizeof(double));
  w1 = (double *)malloc((npts+1)*sizeof(double));
  w2 = (double *)malloc((npts+1)*sizeof(double));
  w3 = (double *)malloc((npts+1)*sizeof(double));
  w_test = (double *)malloc((npts+1)*sizeof(double));
  w_test2 = (double *)malloc((npts+1)*sizeof(double));
  /*call to gauss_legendre()*/
  
  gauss_legendre(x_left, x_right, x, w, npts);
  fprintf(stderr,"OK1");
  quad_log_singularity_half(y,npts,x,w,x1,w1);
  fprintf(stderr,"OK2");
  quad_x_singularity(y,npts,x,w,x2,w2);
  fprintf(stderr,"OK3");
  quad_x2_singularity(y,npts,x,w,x3,w3);
  
  quad_log_singularity_half(0.0,npts,x,w,x_test,w_test);
  z1 = compute_log_integral(npts,x_test,w_test);
  quad_x2_singularity(0.0,npts,x,w,x_test2,w_test2);
  z2 = compute_x2_integral(npts,x_test,w_test);
  f.function = &test_sin;
  f.params = &p;
  z3 = integral_gaussleg(f,0,PI/2.0,npts,x,w);
  /*open file for printing*/
  fp = fopen("gauleg_test", "a");
  if (fp == NULL) {
    fprintf(stderr, "ERROR: could not open file \n");
    exit(0);
  }
  
  fprintf(fp, "abcissas: \n");
  /*print abcissas*/
  for(i = 1; i <= npts; i++)
    {
      fprintf(fp, "%lf ; \n", x[i]);
    }
  
  fprintf(fp, "weights: \n");
  for(i = 1; i <= npts; i++)
    {
      fprintf(fp,"%lf ; \n", w[i]);
    }
  
  fprintf(fp,"\n");
  fprintf(fp, "weights 1: \n");
  for(i = 1; i <= npts; i++)
    {
      fprintf(fp,"%lf ; \n", w1[i]);
    }
  
  fprintf(fp,"\n");
  fprintf(fp,"weights 2: \n");
  for(i = 1; i <= npts; i++)
    {
      fprintf(fp,"%lf ; \n", w2[i]);
    }
  
  fprintf(fp,"\n");
  fprintf(fp,"weights 3: \n");
  for(i = 1; i <= npts; i++)
    {
      fprintf(fp,"%lf ; \n", w3[i]);
    }
  
  fprintf(fp,"\n");
  /* for(i = 1; i <= npts; i++) */
  /*   { */
  /*     fprintf(fp,"%lf ; \n", w_test[i]); */
  /*   } */
  

  fprintf(fp,"\n");
  fprintf(fp,"log integral: %lf",z1);
  fprintf(fp,"\n");
  fprintf(fp,"x2 integral: %lf",z2);
  fprintf(fp,"\n");
  fprintf(fp,"sin integral 0 pi/2: %lf",z3);
  
  free(x);
  free(w);
  free(w1);
  free(w2);
  free(w3);
  free(w_test);
  free(w_test2);
  

  fclose(fp);
  return 0;
}
  
  
double test_sin(double x, void *params)
{
  return sin(x);
}
