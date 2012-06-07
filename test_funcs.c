#include "test_funcs.h"
#include <math.h>
#include <gsl/gsl_math.h>


double f1(double x, void *params)
{
  params_f1 *p2;
  p2 =(params_f1 *)params;
  
  return sin(x);
}

