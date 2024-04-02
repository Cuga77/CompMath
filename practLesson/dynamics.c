#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// dimension 2
// x[0] -- position (height)
// x[1] -- velocity = dx[0]/dt
void f_fall (double *x, double *fx, void *params)
{
  double g = *(double *)params;
  fx[0] = x[1];
  fx[1] = -g;
}

void euler_step (int n, double *x0, double *xh, double h,
		 void (*f)(double *x, double *fx, void *params), void *params)
{
  int i;
  double *fx = (double *)malloc(sizeof(double) * n);

  f(x0, fx, params);
  
  for (i = 0; i < n; i++)
    xh[i] = x0[i] + h * fx[i];
  
  free(fx);
}

void midpoint (int n, double *x0, double *xh, double h,
	       void (*f)(double *x, double *fx, void *params), void *params)
{
  double *fx = (double *)malloc(sizeof(double) * n);
  int i;
  
  f(x0, fx, params);
  
  for (i = 0; i < n; i++)
    xh[i] = x0[i] + 0.5 * h * fx[i]; // actually xmid
 
  f(xh, fx, params);
 
  for (i = 0; i < n; i++)
    xh[i] = x0[i] + h * fx[i];
  
  free(fx);
}

// x0 -- initial value (at t=0)
// T -- time (final)
// h -- time step
void integrate (int n, double *x0, double T, double h,
		void (*f)(double *x, double *fx, void *params),
		void (*step)(int n, double *x0, double *xh, double h,
			     void (*f)(double *x, double *fx, void *params), void *params),
		void *params)
{
  double t = 0;
  double *x = (double *)malloc(sizeof(double) * n);
  double *xh = (double *)malloc(sizeof(double) * n);
  int i;

  memcpy(x, x0, sizeof(double) * n);

  for (t = 0; t < T; t += h)
  {
    printf("%.3le", t);
    for (i = 0; i < n; i++)
      printf(" %.8le", x[i]);
    printf("\n");
    step(n, x, xh, h, f, params);
    memcpy(x, xh, sizeof(double) * n);
  }
  printf("%.3le", t);
  for (i = 0; i < n; i++)
    printf(" %.8le", x[i]);
  printf("\n");
  free(x);
  free(xh);
}
