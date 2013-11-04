#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void printMatrix(double *a, int row, int col) {
  printf("*********** Matrix %d x %d **************\n", row, col);
  int i, j;
  for (i = 0; i < row; i++) {
    for (j = 0; j < col; j++) 
      printf("%15.10f ", a[i*col+j]);
    printf("\n");
  }
  return;
}
double* transpose(double *a, int row, int col) {
  double *result = (double*)malloc(row*col*sizeof(double));
  int i, j;
  for (i = 0; i < row; ++i) {
    for (j = 0; j < col; ++j) {
      result[j*row+i] = a[i*col+j];
    }
  }
  return result;
}
//p*q matrix times q*r matrix, result is p*r matrix
double* multiply(double *a, double *b, int p, int q, int r) {
  double *result = (double*)malloc(p*r*sizeof(double));
  int i, j, k;
  for (i = 0; i < p; i++) {
    for (j = 0; j < r; j++) {
      double sum = 0;
      for (k = 0; k < q; k++)
        sum += a[i*q+k] * b[k*r+j];
      result[i*r+j] = sum;
    }
  }
  return result;
}
double* addition(double *a, double *b, int row, int col) {
  int i, j;
  double *result = (double*)malloc(row*col*sizeof(double));
  for (i = 0; i < row; ++i)
    for (j = 0; j < col; ++j)
      result[i*col+j] = a[i*col+j] + b[i*col+j];
  return result;
}
double* subtraction(double *a, double *b, int row, int col) {
  int i, j;
  double *result = (double*)malloc(row*col*sizeof(double));
  for (i = 0; i < row; ++i)
    for (j = 0; j < col; ++j)
      result[i*col+j] = a[i*col+j] - b[i*col+j];
  return result;
}
double normal(double *a, int row, int col) {
  double sum = 0;
  int i, j;
  for (i = 0; i < row; i++)
    for (j = 0; j < col; j++)
      sum += a[i*col+j] * a[i*col+j];
  return sqrt(sum);
}
double* gradient(double *a, double *x, double *y, int n) {
  double *result = 
    multiply(transpose(a, n, n), subtraction(multiply(a, x, n, n, 1), y, n, 1), n, n, 1);
  int i;
  for (i = 0; i < n; ++i)
    result[i] *= 2;
  return result;
}
void dumb_solve(double *a, double *y, int n, double eps,
		int numit, double *x, int *niter, double *discreps) {
  int i, j;
  *niter = 0;
  //double r = -0.5 / normal(a, n, n);
  //double r = -0.5;
  double bbbb = 0;
  while (*niter < numit) {
    double *x_gradient = gradient(a, x, y, n);
    
    //calculate r to generate next x
    double r = normal(x_gradient, n, 1);
    r = r * r * (-0.5);
    double tmp = normal(multiply(a, x_gradient, n, n, 1), n, 1);
    r = r / (tmp * tmp);
    
    
    for (i = 0; i < n; i++)
      x_gradient[i] *= r;
    double *new_x = addition(x, x_gradient, n, 1);
    discreps[*niter] = normal(subtraction(multiply(a, x, n, n, 1), y, n, 1), n, 1);
    if (discreps[*niter] < eps)
      break;
    memcpy(x, new_x, n*sizeof(double));
    free(new_x);
    *niter += 1;
    //printf("Round %d: r = %lf\n", *niter, r);
    //printMatrix(x, n, 1);
    //getchar();
  }
  return;
}
void problem(int n) {
  double *a = (double *)malloc(n*n*sizeof(double));
  int i, j;
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j) {
      if (i == j)
        a[i*n+j] = 1.0/((i+1)*(i+1));
      else
        a[i*n+j] = 0;
    }

  double *y = (double *)malloc(n * sizeof(double));
  for (i = 0; i < n; ++i)
    y[i] = 1.0;

  int numit = 1000;
  double eps = 1E-6;

  //use to host the result
  double *x = (double *)malloc(n * sizeof(double));
  //initialize x as all zero
  for (i = 0; i < n; ++i)
    x[i] = 0.0;

  int niter = 0;
  double *discreps = (double *)malloc(numit * sizeof(double));
  //initialize discrep as all INT_MAX
  for (i = 0; i < numit; ++i)
      discreps[i] = -1234;

  printf("Matrix A:\n");
  printMatrix(a, n, n);
  printf("Matrix Y:\n");
  printMatrix(y, n, 1);
  //double *y_t = transpose(y, n, 1);
  //printMatrix(y_t, 1, n);
  //printMatrix(multiply(y, y_t, n, 1, n), n, n);
  
  dumb_solve(a, y, n, eps, numit, x, &niter, discreps);

  printf("Solution X: \n");
  printMatrix(x, n, 1);
  printf("Total Round: %d. Discreps: \n", niter);
  //printMatrix(discreps, niter, 1);
  printf("\n\n");


}
int main() {
  problem(3);
  //getchar();
  problem(4);
  //getchar();
  problem(6);
  //getchar();
  problem(10);
}
