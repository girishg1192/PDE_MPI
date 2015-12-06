#include <iostream>
#include "vector.h"

using namespace std;

extern "C" int dgetrf_(int* SIZE1, int* SIZE2, double* a, int* lda, int* ipiv, int* info);
extern "C" void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );
extern "C" void dgesv_(int *TRANS, int *N, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );

int krylov_solver(struct sparse *A, double *rK, int mat_Size);


int main()
{
  std::cout.precision(15);


  int N =10;
  cin>>N;
  int SIZE=N-2, mat_Size= SIZE*SIZE;
#if 0
  int LDA = mat_Size;
  int LDB = mat_Size;
  int ipiv[mat_Size];
  int info;
  int nrhs = 1;
  char trans = 'N';
#endif

  //double **u=(double **)malloc((N*N)*sizeof(double*));
  //init_Matrix(u, N*N);

  double **A = new double*[mat_Size];
  double *B = new double[mat_Size];
  struct sparse A_sparse;

  A_sparse.values= new double[mat_Size*mat_Size];
  A_sparse.col_ind = new int[mat_Size*mat_Size];
  A_sparse.row_ptr = new int[mat_Size];
  A_sparse.nnz = mat_Size;


  init_Matrix(A, mat_Size);
  init_Vector(B, mat_Size);

  generateMatrix(A, SIZE);
  generateVector(B, SIZE);
  generateSparse(&A_sparse);

  //printMatrix(A, mat_Size);
 // printVector(B, mat_Size);
  //printSparse(&A_sparse, mat_Size);

  //cout<<endl;
  krylov_solver(&A_sparse, B, mat_Size);
  //cout<<endl<<endl;
 // printVector(B,mat_Size);
 
  delete[] A;
  delete[] B;
  return(0);
}

int krylov_solver(struct sparse *A, double *rK, int mat_Size)   //b passed to rk, for first step residual (r0)= b 
{
  bool is_converge = false;
  double alpha_k, beta_k;
  double *p_k = new double[mat_Size];

  double *result_k = new double[mat_Size];    //Stores xk for kth step
  memset(result_k, 0, mat_Size*sizeof(double));
  double *a_norm = new double[mat_Size];
  memset(a_norm, 0, mat_Size*sizeof(double));

#ifndef PRECOND
  double **M_precondition = new double*[mat_Size];    //declare and initialize preconditioning matrix
  init_Matrix(M_precondition, mat_Size);
  generatePreconditioner(M_precondition, mat_Size);

  double *z_k = new double[mat_Size];
  memset(z_k, 0, mat_Size*sizeof(double));

  mat_vector_mult(M_precondition, rK, mat_Size, &z_k);

  memcpy(p_k, z_k, mat_Size*sizeof(double));   //Stores minimizer for kth step pk
  int k;
  for(k =0; k<1000 && !is_converge; k++)
  {
    memset(a_norm, 0, mat_Size*sizeof(double));
    //mat_vector_mult(A, p_k, mat_Size, &a_norm);   //Stores the result of the matrix vector multiplication
    sparse_mat_vector(A, p_k, mat_Size, &a_norm);

    double rK_dot = vectorDot(rK, z_k, mat_Size);    //Stores the value of the dot product of the residual

    alpha_k = rK_dot/(vectorDot(p_k, a_norm, mat_Size));
    for(int i=0; i<mat_Size; i++)
    {
      result_k[i] = result_k[i] + alpha_k*p_k[i];
      rK[i] = rK[i] - alpha_k*a_norm[i];
    }
    mat_vector_mult(M_precondition, rK, mat_Size, &z_k);   //Stores the result of the matrix vector multiplication
    beta_k = vectorDot(z_k, rK, mat_Size)/rK_dot;
    for(int i=0; i<mat_Size; i++)
      p_k[i] = z_k[i] + beta_k*p_k[i];
    is_converge = checkConvergence(rK, mat_Size);
  }
  cout<<k<<" iter\n";
  printVector(rK, mat_Size);
  delete[] M_precondition;
#else
  int k;

  memcpy(p_k, rK, mat_Size*sizeof(double));   //Stores minimizer for kth step pk

  for(k =0; k<100 && !is_converge; k++)
  {
    memset(a_norm, 0, mat_Size*sizeof(double));
    //mat_vector_mult(A, p_k, mat_Size, &a_norm);   //Stores the result of the matrix vector multiplication
    sparse_mat_vector(A, p_k, mat_Size, &a_norm);

    double rK_dot = vectorDot(rK, rK, mat_Size);    //Stores the value of the dot product of the residual

    alpha_k = rK_dot/(vectorDot(p_k, a_norm, mat_Size));
    for(int i=0; i<mat_Size; i++)
    {
      result_k[i] = result_k[i] + alpha_k*p_k[i];
      rK[i] = rK[i] - alpha_k*a_norm[i];
    }
    beta_k = vectorDot(rK, rK, mat_Size)/rK_dot;
    for(int i=0; i<mat_Size; i++)
      p_k[i] = rK[i] + beta_k*p_k[i];
    is_converge = checkConvergence(rK, mat_Size);
  }
  cout<<k<<" itern\n";
#endif

  //printVector(rK, mat_Size);
  printVectorMat(result_k, mat_Size);
  cout<<endl;
  delete[] p_k;
  delete[] result_k;
  delete[] a_norm;
  return -1;
}
