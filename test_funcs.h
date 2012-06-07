#ifndef __TEST_FUNCS__H__
#define __TEST_FUNCS__H__
#include <gsl/gsl_math.h>

typedef struct params_f1{
  int i;
}params_f1;

gsl_function f;
double f1(double x, void *params);

/* params_f1 p = {0}; */



#endif /*__TEST_FUNCS__H__*/
