#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void f_decay (double *x, double *fx, void *params)
{
  double k = *(double *)params;
  *fx = -k * *x;
}

// dimension 2
void f_merrygoround (double *x, double *fx, void *params)
{
  double k = *(double *)params;
  fx[0] = -x[1] * k;
  fx[1] =  x[0] * k;
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

int main (int argc, char **argv)
{
  double T, h;
  double k;
  double point[2] = {10, 0};
  
  if (argc < 4)
  {
    fprintf(stderr, "Usage: prog2 <K> <T> <h>\n");
    return -1;
  }

  if (sscanf(argv[1], "%lf", &k) < 1 ||
      sscanf(argv[2], "%lf", &T) < 1 ||
      sscanf(argv[3], "%lf", &h) < 1)
  {
    fprintf(stderr, "error reading values\n");
    return -1;
  }
  
  integrate(2, point, T, h, f_merrygoround, midpoint, &k);
  return 0;
}
