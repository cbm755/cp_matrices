#include <cuda.h>

/*
 * Kernel for Sparse Matrix Vector Multiply
 *
 * This takes an N x N sparse matrix and multiplies it by an N-vector x.
 *
 * The matrix has exactly M entries per row (some could be zero) and is stored
 * in the 2D matrix As and Aj which contain the coefficients and column numbers
 * respectively.  Aj is assumed to count from 1 not zero (e.g., it comes from Matlab)
 * The result is returned in y.
 *
 * Colin Macdonald, 2013-01-13
 */
__global__ void kernel_spmatvec2(double * y, const double * As, const int * Aj,
                                 const double * x, const int N, const unsigned int M)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int k;
  unsigned int j;
  if (idx < N) {
    y[idx] = 0.0;   // todo: any advantage to adding into a temporary double?
    for (j = 0; j < M; j++) {
      k = Aj[idx + j*N] - 1;  // -1 here for matlab indexing
      y[idx] += As[idx + j*N] * x[k];
    }
  }
}

