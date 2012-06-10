
#include <stdio.h>
#include <stdlib.h>
#include "integ_hyp.h"
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_legendre.h>


/*Limits of integration: x1, x2
  Returns: x[1...n] and w[1...n], abcissas and weights for Gauss Legendre Quadrature
  Tested OK on 3-2-2012*/

void gauss_legendre(double x1, double x2, double* x, double* w, int n)
{
  int m, j, i;
  double z, z_temp, x_mid, x_length, p1, p2, p3, pd;
  
  m = (n+1)/2; /*roots are symmetric in interval, only half must be found*/
  x_mid = 0.5*(x1 + x2);
  x_length = 0.5*(x2 - x1);

  for(i = 1; i <= m; i++) /*loop over roots*/
    { 
      /*initial approximation to root*/
      
      z = cos(PI*(i - 0.25)/(n + 0.5));
      
      /* main loop of newton's method*/
      
      do {
	p1 = 1.0;
	p2 = 0.0;
	/*looping up the recurrence relation to get legendre poly evaluated at z*/
	for(j = 1; j<=n; j++) 
	  {
	    p3 = p2;
	    p2 = p1;
	    
	    p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j;
	  }
	/*p1 is now the desired value of the appropriate legendre poly.*/
	
	/*computing derivative pd*/

	pd = n*(z*p1 - p2)/(z*z - 1.0);
	z_temp = z;
	z = z_temp - (p1/pd);
      } while(fabs(z - z_temp) > PREC);
      
      x[i] = x_mid - (x_length*z);
      x[n+1 -i] = x_mid + (x_length*z);
      
      /*compute weights*/
      
      w[i] = 2.0*x_length/((1.0 - z*z)*pd*pd);
      w[n+1 - i] = w[i];
    }

}


double legendre_poly(int deg, double z)
{
  int i;
  double p1, p2, p3;
  
  p1 = 1.0;
  p2 = 0.0;
  /*looping up the recurrence relation to get legendre poly evaluated at z*/
  for(i = 1; i<=deg; i++) 
    {
      p3 = p2;
      p2 = p1;
	    
      p1 = ((2.0*i - 1.0)*z*p2 - (i - 1.0)*p3)/i;
    }
  /*p1 is now the desired value of the appropriate legendre poly.*/
  
  return p1;
}


double r_func(int i, double x)
{
  return (gsl_sf_legendre_Ql(i,x) + ((0.25)*log((x-1.0)*(x-1.0))));
}

int quad_log_singularity_half(double y, int n, double *x_gauss, double *w_gauss, double *x, double *w)
{
  int i, j;
  double z;
  for(i=1;i<=n;i++)
    {
      x[i] = x_gauss[i];
    }


  for(i=1;i<=n;i++)
    {
      z = 0.0;
      for(j=1;j<=n-2;j++)
	{
	  z = z + ((legendre_poly(j-1,x_gauss[i]) - legendre_poly(j+1,x_gauss[i]))*r_func(j,y));
	}
      
      z = z + ((legendre_poly(0,x_gauss[i])-legendre_poly(1,x_gauss[i]))*r_func(0,y));
      
      z = z + (legendre_poly(n-2,x_gauss[i])*r_func(n-1,y)) + (legendre_poly(n-1,x_gauss[i])*r_func(n,y));
      
      w[i] = w_gauss[i]*z;
    }
  return 0;
}

int quad_x_singularity(double y, int n, double *x_gauss, double *w_gauss, double *x, double *w)
{
  int i, j;
  double z;
  
  for(i=1;i<=n;i++)
    {
      x[i] = x_gauss[i];
    }

  for(i=1;i<=n;i++)
    {
      z = 0.0;
      for(j=0;j<=n-1;j++)
	{
	  z = z + (((2*j) + 1)*legendre_poly(j,x_gauss[i])*gsl_sf_legendre_Ql(j,y));
	}
      
      w[i] = w_gauss[i]*z;
    }
  return 0;
}


int quad_x2_singularity(double y, int n, double *x_gauss, double *w_gauss, double *x, double *w)
{
  int i, j, k, n1, n2, n3;
  double z, temp, z_temp;
  for(i=1;i<=n;i++)
    {
      x[i] = x_gauss[i];
    }


  for(i=1;i<=n;i++)
    {
      z = 0.0;
      
      z_temp = 0.0;
      
      for(j=0;j<=(n-2);j++)
	{
	  temp = floor((j-1)/2.0);
	  n3 = (int)temp;
	  for(k=0;k<=n3;k++)
	    {
	      z_temp = z_temp + (((2*j) - (4*k) -1)*((2*j) + 1)*gsl_sf_legendre_Ql(j-1-(2*k),y)*legendre_poly(j,x_gauss[i]));
	    }
	}
      z = z - z_temp;
      
      for(j=0;j<=n-1;j++)
	{
	  if((j%2)==0){
	    z = z + (((2.0*j + 1.0)/2.0)*legendre_poly(j,x_gauss[i])*((1.0/(y-1.0))-(1.0/(y+1.0))));
	  } else {
	    z = z + (((2.0*j + 1.0)/2.0)*legendre_poly(j,x_gauss[i])*((1.0/(y-1.0))+(1.0/(y+1.0))));
	  }
	}
      w[i] = w_gauss[i]*z;
    }
  return 0;
}


/*testing quad log integral function*/

double compute_log_integral(int n, double *x, double *w)
{
  /*compute the integral of .5*log((x-y)^2), y =0 on -1 to 1*/

  int i;
  double z;

  z = 0.0;
  for(i=1;i<=n;i++)
    {
      /*fprintf(stdout,"%lf",z);*/
      z = z + w[i];
    }
  return z;
}


/*test x2 singularity integral*/

double compute_x2_integral(int n, double *x, double *w)
{
  int i;
  double z;

  z = 0.0;
  for(i=1;i<=n;i++)
    {
      /*fprintf(stdout,"%lf",z);*/
      z = z + w[i];
    }
  return z;
}

double integral_gaussleg(gsl_function f, double a, double b, int n, double *x, double *w)
{
  int i;
  double z, m, c;
  double *x_t;
  
  x_t = (double *)malloc((n+1)*sizeof(double));

  /*transform to -1 1*/
  
  m = (b-a)/2.0;
  c = (b+a)/2.0;
  for(i=1;i<=n;i++)
    {
      x_t[i] = m*x[i] + c;
    }
  
  /*integrate*/
  z = 0.0;
  for(i=1;i<=n;i++)
    {
      z = z + w[i]*(*(f.function))(x_t[i],f.params);
    }
  z = m*z;
  
  free(x_t);

  return z;
}


double integral_gauss_log_sing(gsl_function f, double a, double b, double y, int n, double *x_gauss, double *w_gauss, double*x, double *w)
{
  int i;
  double z, m, c,yscale;
  double *x_t;

  m = (b-a)/2.0;
  c = (b+a)/2.0;
  yscale = (1/m)*y - (c/m);

  if(w==NULL&x==NULL){
    x = (double *)malloc((n+1)*sizeof(double));
    w = (double *)malloc((n+1)*sizeof(double));
    for(i=1;i<=n;i++)
      {
	x[i] = x_gauss[i];
      }


    quad_log_singularity_half(yscale, n, x_gauss, w_gauss, x, w);
  }
  x_t = (double *)malloc((n+1)*sizeof(double));

  
  
  /*transform to -1 1*/
  


  for(i=1;i<=n;i++)
    {
      x_t[i] = m*x[i] + c;
    }
  
  /*generate weights based on point y*/
  
  

  /*integrate remember integral returned is of f(x)log|x-y|*/
  z = 0.0;
  for(i=1;i<=n;i++)
    {
      z = z + w[i]*(*(f.function))(x_t[i],f.params);
    }
  z = m*z;
  
  free(x_t);

  return z;
}
		    
double integral_gauss_x_sing(gsl_function f, double a, double b, double y, int n, double *x_gauss, double *w_gauss, double*x, double *w)
{
  
  int i;
  double z, m, c,yscale;
  double *x_t;
  m = (b-a)/2.0;
  c = (b+a)/2.0;
  
  yscale = (1/m)*y - (c/m);

  if((x==NULL)&&(w==NULL)){
    x = (double *)malloc((n+1)*sizeof(double));
    w = (double *)malloc((n+1)*sizeof(double));
    /*copy x_gauss to x*/
    for(i=1;i<=n;i++)
      {
	x[i] = x_gauss[i];
      }
    quad_x_singularity(yscale, n, x_gauss, w_gauss, x, w);

  }
  x_t = (double *)malloc((n+1)*sizeof(double));
  
  
  /*transform to -1 1*/
  
  for(i=1;i<=n;i++)
    {
      x_t[i] = m*x[i] + c;
    }
  
  /*generate weights based on point y*/
  


  /*integrate remember integral returned is of f(x)log|x-y|*/
  z = 0.0;
  for(i=1;i<=n;i++)
    {
      z = z + w[i]*(*(f.function))(x_t[i],f.params);
    }
  z = m*z;
  
  free(x_t);

  return z;
}
  
double integral_gauss_x2_sing(gsl_function f, double a, double b, double y, int n, double *x_gauss, double *w_gauss, double*x, double *w)
{
  int i;
  double z, m, c,yscale;
  double *x_t;
  c = (b+a)/2.0;
  m = (b-a)/2.0;

  yscale = (1/m)*y - (c/m);

  if((x==NULL)&&(w==NULL)){
    x = (double *)malloc((n+1)*sizeof(double));
    w = (double *)malloc((n+1)*sizeof(double));
    /*copy x_gauss to x*/
    for(i=1;i<=n;i++)
      {
	x[i] = x_gauss[i];
      }
    quad_x2_singularity(yscale, n, x_gauss, w_gauss, x, w);
  }
  x_t = (double *)malloc((n+1)*sizeof(double));
  
  
  /*transform to -1 1*/
  
  
  
  for(i=1;i<=n;i++)
    {
      x_t[i] = m*x[i] + c;
    }
  
  /*generate weights based on point y*/
  
  

  
  z = 0.0;
  for(i=1;i<=n;i++)
    {
      z = z + w[i]*(*(f.function))(x_t[i],f.params);
    }
  z = m*z;
  
  free(x_t);

  return z;
}

double func_legendre_pl(double z, void *params)
{
  pl_params *p;
  int i, deg;
  double p1, p2, p3;
 
  p = (pl_params *)params;

  deg = p->deg;
  return legendre_poly(deg,z);
}
  
  
  
double integral_generalized_sing(gsl_function f, double a, double b, double y, int n, int m, double *x_gauss, double *w_gauss, double*x, double *w)
{
  gsl_vector *rhs,*soln,*tau,*res;
  gsl_matrix *A;
  double *w1,*w2,*w3,*x1,*x2,*x3,*x_t;
  gsl_function f_temp;
  pl_params p;
  double y_scaled,x_mid,x_halflength,integral;
  int i,j;

  
  x_mid = (b+a)/2.0;
  x_halflength = (b-a)/2.0;

  
  /*scale y to -1 to 1*/
  y_scaled = (1/x_halflength)*y - (x_mid/x_halflength);
  
  /*allocate memory for aux. quadrature weight calculations*/

  w1 = (double *)malloc((n+1)*sizeof(double));
  w2 = (double *)malloc((n+1)*sizeof(double));
  w3 = (double *)malloc((n+1)*sizeof(double));
  x1 = (double *)malloc((n+1)*sizeof(double));
  x2 = (double *)malloc((n+1)*sizeof(double));
  x3 = (double *)malloc((n+1)*sizeof(double));
  x_t = (double *)malloc((n+1)*sizeof(double));
  /*allocate memory for system matrix rhs vector and solution vector along with
    vector of householder coeffs*/

  rhs = gsl_vector_calloc(4*m);
  res = gsl_vector_calloc(4*m);
  soln = gsl_vector_calloc(n);

  /*Note that here we assume that 4*m > n*/
  tau = gsl_vector_calloc(n);
  A = gsl_matrix_calloc(4*m,n);


  /*fill in the entries of the matrix*/

  for(i=0;i<n;i++)
    {
      for(j=0;j<m;j++)
	{
	  gsl_matrix_set(A,j,i,legendre_poly(j,x_gauss[i+1]));
	  gsl_matrix_set(A,j+m,i,legendre_poly(j,x_gauss[i+1])*log(fabs(y_scaled - x_gauss[i+1])));
	  gsl_matrix_set(A,j+(2*m),i,legendre_poly(j,x_gauss[i+1])*(1/(y_scaled - x_gauss[i+1])));
	  gsl_matrix_set(A,j+(3*m),i,legendre_poly(j,x_gauss[i+1])*(1/((y_scaled - x_gauss[i+1])*(y_scaled - x_gauss[i+1]))));
	  
	}
    }

  /*calculate quadrature points for the exact evaluation of the integrals of the various
    phi functions */
  

  quad_log_singularity_half(y_scaled, n, x_gauss, w_gauss, x1, w1);
  quad_x_singularity(y_scaled, n, x_gauss, w_gauss, x2, w2);
  quad_x2_singularity(y_scaled, n, x_gauss, w_gauss, x3, w3);

  /*fill in the rhs vector*/

  for(i=0;i<m;i++)
    {
      
      p.deg = i;
      f_temp.function = &func_legendre_pl;
      f_temp.params = &p;
      
      gsl_vector_set(rhs,i,integral_gaussleg(f_temp,-1.0,1.0,n,x_gauss,w_gauss));
      gsl_vector_set(rhs,i+m,integral_gauss_log_sing(f_temp,-1.0,1.0,y_scaled,n,x_gauss,w_gauss,x1,w1));
      gsl_vector_set(rhs,i+(2*m),integral_gauss_x_sing(f_temp,-1.0,1.0,y_scaled,n,x_gauss,w_gauss,x2,w2));
      gsl_vector_set(rhs,i+(3*m),integral_gauss_x2_sing(f_temp,-1.0,1.0,y_scaled,n,x_gauss,w_gauss,x3,w3));
    }
  
  /*solve the linear system in the least squares sense*/
  
  gsl_linalg_QR_decomp(A,tau);
  gsl_linalg_QR_lssolve(A,tau,rhs,soln,res);
  
  for(i=0;i<n;i++)
    {
      w[i+1] = gsl_vector_get(soln,i);
      x[i+1] = x_gauss[i+1];
    }

  /*integrate*/
  
  for(i=1;i<=n;i++)
    {
      x_t[i] = x_halflength*x[i] + x_mid;
    }
  
  integral = 0;
  for(i=1;i<=n;i++)
    {
      integral = integral + (w[i]*(*(f.function))(x_t[i],f.params));
    }
  integral = x_halflength*integral;

  /*free allocated memory*/
  free(x1);
  free(x2);
  free(x3);
  free(w1);
  free(w2);
  free(w3);
  free(x_t);
  gsl_vector_free(rhs);
  gsl_vector_free(soln);
  gsl_vector_free(tau);
  gsl_vector_free(res);
  gsl_matrix_free(A);
  
  
  
  return integral;
}

  
		     
			 
  
  
      
  
	
  
  








