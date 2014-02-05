#include <cuda.h>

/*
 * Kernel for evaluation of f and g.
 *
 * This takes in vectors u & v of length N, and returns f & g.
 *
 */
__global__ void kernel_fg(double * f, double * g, const double * u, const double * v, const int N)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < N){
    f[idx] = 0.0;   
    f[idx] += -u[idx]*v[idx]*v[idx] + 0.054*(1-u[idx]);
    g[idx] = 0.0;
    g[idx] += u[idx]*v[idx]*v[idx] - (0.054+0.063)*v[idx];
  }
}


