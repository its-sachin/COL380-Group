#include <stdio.h>

__global__
void saxpy(int n, float a, float *x, float *y, int* res)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  printf("%d\n", i);
  res[i] = 1;
  if (i < n) y[i] = a*x[i] + y[i];
}

int main(void)
{
  int N = 100;
  float *x, *y, *d_x, *d_y;
  
  x = (float*)malloc(N*sizeof(float));
  y = (float*)malloc(N*sizeof(float));

  cudaMalloc(&d_x, N*sizeof(float)); 
  cudaMalloc(&d_y, N*sizeof(float));


  for (int i = 0; i < N; i++) {
    x[i] = 1.0f;
    y[i] = 2.0f;
  }

  cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, N*sizeof(float), cudaMemcpyHostToDevice);

  int *result = new int[N];
  memset(result, 0, N*sizeof(int));
  int *d_result;
  cudaMalloc(&d_result, N*sizeof(float));
  cudaMemcpy(d_result, result, N*sizeof(float), cudaMemcpyHostToDevice);

  // Perform SAXPY on 1M elements
  saxpy<<<(N+255)/256, 256>>>(N, 2.0f, d_x, d_y, d_result);

  cudaMemcpy(result, d_result, N*sizeof(int), cudaMemcpyDeviceToHost);

  cudaMemcpy(y, d_y, N*sizeof(float), cudaMemcpyDeviceToHost);

  for(int i=0; i<N; i++){
    printf("%d\n", result[i]);
  }

  cudaFree(d_x);
  cudaFree(d_y);
  free(x);
  free(y);
}