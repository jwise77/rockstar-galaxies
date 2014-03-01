#include <math.h>
#include <inttypes.h>
#include <stdio.h>
#include <assert.h>
#include "universal_constants.h"
#include "config_vars.h"
#include "universe_time.h"
#include "hubble.h"
#include "integrate.h"

#define STEPS 4096
#define MAXA  5.0

static double times[STEPS+1]={0};
static double H_CONV;

void init_time_table(void) {
  double a;
  int64_t i;
  
  H_CONV = HUBBLE_TIME_CONVERSION/h0;
      
  for (i=STEPS; i>=0; --i) {
    a = MAXA*((double)i)/((double)STEPS);
    times[i] = exact_scale_to_time(a);
  }
}

double inv_hubble_scaling(double a, void *extra_data) {
  return (1.0/(a*hubble_scaling(1.0/a - 1.0)));
}

double exact_scale_to_time(double scale)
{  
  if (scale == 1.0) return 0.0;
  if (scale < 1e-30) scale = 1e-30;
  double tol = fabs(scale-1.0)*1e-7;
  if (tol>1e-7) tol = 1e-7;
  return (adaptiveSimpsons(inv_hubble_scaling, NULL, 1.0, scale, tol, 20));  
}

// Linearly interpolate between calculated values.
double scale_to_time(double scale) {
  double s = scale;
  int64_t l = (int)(s/MAXA*STEPS);
  double f = s/MAXA*STEPS - l;
  if (scale > MAXA) return exact_scale_to_time(scale);
  if (scale < 0) return times[0];
  return (times[l]+f*(times[l+1]-times[l]));
}

double scale_to_years(double scale) {
  return (scale_to_time(scale)*H_CONV);
}

double _exact_time_to_scale(double t) {
  assert((W0 == -1) && (WA == 0));
  double m = sinh(1.5*t*sqrt(1-Om));
  return pow(Om*m*m/(1.0-Om), 1.0/3.0);
}

/*int main(void) {
  double fact = 10.0;
  double i;
  h0 = 0.7;
  Om = 0.7;
  Ol = 0.3;
  W0 = -1.0;
  WA = 0.0;
  init_time_table();
  
  double exact_t0_conv = -1.0*exact_scale_to_time(0.0);
  
  for (i=0; i<100; i++) {
    double a = _exact_time_to_scale(i/fact);
    printf("a = %10f t = %10f exact func = %10f exact2 func = %10f interp func = %10f (frac error in time = %g)\n", 
	   a, i/fact-exact_t0_conv, exact_scale_to_time(a), exact_scale_to_time(a), scale_to_time(a), scale_to_time(a)/(i/fact-exact_t0_conv)-1.0);
  }

  printf("%f\n", fact);
  
  fprintf(stderr,"age of universe = %f\n",exact_t0_conv);
  fprintf(stderr,"scale(age of universe) = %f\n",_exact_time_to_scale(exact_t0_conv));
}
*/
