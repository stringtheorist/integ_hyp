#ifndef __TEST_UTILS__H__
#define __TEST_UTILS__H__
#include <gsl/gsl_math.h>
#include "integ_hyp.h"

int integral_generalized_test(gsl_function f, double a, double b, double y,int low_lim_n, int up_lim_n, int exact_n, double ratio, FILE *fp);

int integral_generalized_benchmark(gsl_function f, double a, double b, double y, int low_lim_n, int up_lim_n, int reps, FILE *fp);

/* int integral_log_test(gsl_function f, double a, double b, double y, int low_lim_n, int up_lim_n, FILE *fp); */

/* int integral_x_test(gsl_function f, double a, double b, double y, int low_lim_n, int up_lim_n, FILE *fp); */

/* int integral_x2_test(gsl_function f, double a, double b, double y, int low_lim_n, int up_lim_n, FILE *fp); */




#endif /*__TEST_UTILS__H__*/
