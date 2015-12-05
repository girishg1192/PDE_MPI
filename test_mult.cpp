#include <iostream>
#include "vector.h"

using namespace std;

extern "C" int dgetrf_(int* SIZE1, int* SIZE2, double* a, int* lda, int* ipiv, int* info);
extern "C" void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );
extern "C" void dgesv_(int *TRANS, int *N, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );

int krylov_solver(double **A, double *rK, int mat_Size);

int main()
{
  char trans = 'N';
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
#endif

  //double **u=(double **)malloc((N*N)*sizeof(double*));
  //init_Matrix(u, N*N);

  double **A = new double*[mat_Size];
  double *B = new double[mat_Size];


  init_Matrix(A, mat_Size);
  init_Vector(B, mat_Size);

  generateMatrix(A, SIZE);
  generateVector(B, SIZE);

 // printMatrix(A, mat_Size);
 // printVector(B, mat_Size);

#if 0
  //dgetrf_(&mat_Size, &mat_Size, (double *)A, &LDA, ipiv, &info);
//  std::cout << "Info = " << info << std::endl;
  //dgetrs_(&trans, &mat_Size, &nrhs, (double *)A, &LDA, ipiv, (double *)B, &LDB, &info);
//  dgesv_(&mat_Size, &mat_Size, (double *)A, &LDA, ipiv, (double *)B, &LDB, &info);
#endif
  //cout<<endl;
  krylov_solver(A, B, mat_Size);
  //cout<<endl<<endl;
 // printVector(B,mat_Size);

  for(int i=0; i<mat_Size; i++)
  {
  //  delete[] A[i];
  }
  //delete [] A;
  //free(u);
  return(0);
}

int krylov_solver(double **A, double *rK, int mat_Size)   //b passed to rk, for first step residual (r0)= b 
{
  double alpha_k, beta_k;
  double *p_k = new double[mat_Size];
  memcpy(p_k, rK, mat_Size*sizeof(double));   //Stores minimizer for kth step pk

  double *result_k = new double[mat_Size];    //Stores xk for kth step
  memset(result_k, 0, mat_Size*sizeof(double));
  double *a_norm = new double[mat_Size];
  memset(a_norm, 0, mat_Size*sizeof(double));


  for(int k =0; k<150; k++)
  {
    memset(a_norm, 0, mat_Size*sizeof(double));
    mat_vector_mult(A, p_k, mat_Size, &a_norm);   //Stores the result of the matrix vector multiplication
 //   printVector(a_norm, mat_Size);
    //cout<<endl;

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
  }

  printVector(result_k, mat_Size);
  //printVectorMat(result_k, mat_Size);
  return -1;
}
