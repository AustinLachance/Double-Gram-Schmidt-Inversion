/****************************************************************
Name: Austin Lachance
Email: austin.lachance@yale.edu
Date: 02/26/16

CPSC440 
Double Gram-Schmidt Iversion
Description: This program inverts a matrix A. It functions by
taking the identity matrix and applying to it the same procedures
as are applied to matrix A in the Gram-Schmidt process, creating
a matrix G. This matrix is multiplied by the transpose of matrix
U (where U is the matrix such that A*G = U) to give the inverse
of matrix A. 
****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void inv_double_gs(double* a, int n, double* u, double* b);
//Creates an identity matrix of size n
double *identity(int n);
//prints a square matrix of size n
void matrix_print(double* matrix, int n);
//transposes matrix of size n
double *transpose(double *matrix, int n);
//finds the dot product of two columns in a matrix of size n
double dot_product(double *matrix, int column1, int column2, int n);
//transforms matrix U into the gram schmidt of matrix A
void gram_schmidt(double *a, double *u, double *g, int n);
//Multiplies matrix1 and matrix2 of size n
double *matrix_multiply(double *matrix1, double *matrix2, int n);
//finds the magnitude of a vector "column" in a matrix of size n
double magnitude(double *matrix, int column, int n);



int main() {

  int n = 3;
  double *a = (double*)malloc(n*n*sizeof(double));//input matrix to be inverted
  double *u = (double*)calloc(n*n,sizeof(double));//orthogonal matrix
  double *b = (double*)malloc(n*n*sizeof(double));//inverse of A

  a[0] = 1;
  a[1] = 2;
  a[2] = 2;
  a[3] = -1;
  a[4] = 0;
  a[5] = 2;
  a[6] = 0;
  a[7] = 0;
  a[8] = 1;

  
  inv_double_gs(a, n, u, b);

   printf("A:\n");
  matrix_print(a,n);

  printf("B:\n");
  matrix_print(b,n);

  printf("U:\n");
  matrix_print(u,n);

  return 0;
}

void inv_double_gs(double* a, int n, double* u, double* b) {
  double *g = identity(n);
  gram_schmidt(a, u, g, n);

  double *u_tran = transpose(u, n);
  double *test = matrix_multiply(u, u_tran, n);
  matrix_print(test, n);
  
  double *output = matrix_multiply(g, u_tran, n);
  for(int i = 0; i < n*n; i++) {
    b[i] = output[i];
  }
  free(output);
  free(u_tran);
  free(g);
}

void gram_schmidt(double *a, double *u, double *g, int n) {
  int col = 0;
  for(int row=0; row < n*n; row++) {
    u[row] = a[row];
  }
  for(col=0; col < n; col++) {
    int iteration =0;
    int column = 0;
    while(column < col) 
    {
      double scalar = dot_product(u, col, column, n);
      int row=0;
      for(row = 0; row < n; row++) 
      {
        u[row*n + col] -= scalar * u[row*n + column];
        g[row*n+col] -= scalar* g[row*n + column];
      }   
      column++;
    }
    column = col + 1;
    while(column < n) {
      double scalar = dot_product(u, col, column, n);
      int row=0;
      for(row = 0; row < n; row++) {
        u[row*n + column] -= scalar * u[row*n + col];
        g[row*n+column] -= scalar* g[row*n + col];
      }   
      column++;
    }
    double u_magnitude = magnitude(u, col, n);
    int row=0;
    for(row=0; row < n; row++) 
    {
      u[row*n + col] = u[row*n + col] / u_magnitude;
      g[row*n + col] = g[row*n + col] / u_magnitude;
    }
  } 
}


double *identity(int n) {
  double *identity = (double*)calloc(n*n, sizeof(double));
  int i = 0;
  for(i = 0; i < n; i++) {
    identity[i*n + i] = (double)1;
  }
  return identity;
}

void matrix_print(double *matrix, int n) {
  int i = 0;
  for(i = 0; i < n; i++) {
    int j = 0;
    for(j = 0; j < n; j++) {
      printf("%e\t", matrix[i*n + j]);
    }
    printf("\n");
  }
} 

double *transpose(double *matrix, int n) {
  double *transpose = (double*)malloc(n*n*sizeof(double));
  int i = 0;
  for(i = 0; i < n; i++) {
    int j = 0;
    for(j = 0; j < n; j++) {
      transpose[i*n + j] = matrix[j*n + i];
    }
  }
  return transpose;
}

double dot_product(double *matrix, int column1, int column2, int n) {
  double total = 0.0;
  int i =0;
  for(i = 0; i < n; i++) {
    total += matrix[i*n + column1] * matrix[i*n + column2];
  }
  return total;
}

double magnitude(double *matrix, int column, int n) {
  double total = 0.0;
  int i = 0;
  for(i=0; i < n; i++) {
    total += matrix[i*n + column] * matrix[i*n + column];
  }
  double magnitude = sqrt(total);
  return magnitude;
}

double *matrix_multiply(double *matrix1, double *matrix2, int n){
  double* output = (double*)calloc(n*n,sizeof(double));
  int row = 0;
  for(row = 0; row < n; row++) {
    int col = 0;
    for(col=0; col < n; col++) {
      int add_row = 0;
      for(add_row = 0; add_row < n; add_row++) {
        output[row*n + col] += matrix1[row*n + add_row] *\
         matrix2[add_row*n + col];
      }
    }
  }
  return output;
}


