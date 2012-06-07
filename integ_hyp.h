#ifndef __INTEG_HYP__H__
#define __INTEG_HYP__H__
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>
#define EPS 3.0e-15 /* relative precision */
#define PI 3.141592654



typedef struct pl_params
{
  int deg;
}pl_params;

void gauss_legendre(double, double, double*, double*, int);
double legendre_poly(int, double);
double func_legendre_pl(double, void *);
int quad_log_singularity_half(double y, int n, double *x_gauss, double *w_gauss, double *x, double *w);
int quad_x_singularity(double y, int n, double *x_gauss, double *w_gauss, double *x, double *w);
int quad_x2_singularity(double y, int n, double *x_gauss, double *w_gauss, double *x, double *w);
double compute_log_integral(int n, double *x, double *w); 
double compute_x2_integral(int n, double *x, double *w);
double integral_gaussleg(gsl_function f, double a, double b, int n, double *x, double *w);
double integral_gauss_log_sing(gsl_function f, double a, double b, double y, int n, double *x_gauss, double *w_gauss, double*x, double *w);
double integral_gauss_x_sing(gsl_function f, double a, double b, double y, int n, double *x_gauss, double *w_gauss, double*x, double *w);
double integral_gauss_x2_sing(gsl_function f, double a, double b, double y, int n, double *x_gauss, double *w_gauss, double*x, double *w);

double integral_generalized_sing(gsl_function f, double a, double b, double y, int n, int m, double *x_gauss, double *w_gauss, double *x, double *w);

#endif /* __INTEG_HYP__H__*/

