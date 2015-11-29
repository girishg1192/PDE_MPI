#include <iostream>
#include "vector.h"

using namespace std;

extern "C" int dgetrf_(int* SIZE1, int* SIZE2, double* a, int* lda, int* ipiv, int* info);
extern "C" void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );
extern "C" void dgesv_(int *TRANS, int *N, double *A, int *LDA, int *IPIV, double *B, int *LDB, int *INFO );


int main()
{
  char trans = 'N';
  std::cout.precision(15);


  int N =10;
  cin>>N;
  int info;
  int nrhs = 1;
  int SIZE=N-2, mat_Size= SIZE*SIZE;
  int LDA = mat_Size;
  int LDB = mat_Size;
  int ipiv[mat_Size];

  //vector<double> a, b;

  //double **u=(double **)malloc((N*N)*sizeof(double*));
  //init_Matrix(u, N*N);

  //double **A=(double **)malloc(mat_Size*sizeof(double*));
  double **A = new double*[mat_Size];
  //double *B = (double *)malloc(mat_Size*sizeof(double));
  double *B = new double[mat_Size];


  init_Matrix(A, mat_Size);
  init_Vector(B, mat_Size);

  generateMatrix(A, SIZE);
  generateVector(B, SIZE);

  printMatrix(A, mat_Size);
  printVector(B, mat_Size);
  cout<<endl;

  cout<<endl;

//  for(int i=0; i<mat_Size; i++)
//  {
//    double sum=0.0;
//    for(int j=0;j<mat_Size;j++)
//    {
//      sum+=A[i][j]*B[j];
//    }
//    cout<<sum<<endl;
//  }
//  cout<<endl;


  cout<<LDA<<" "<<LDB<<" "<<mat_Size<<endl;
//  dgetrf_(&mat_Size, &mat_Size, (double *)A, &LDA, ipiv, &info);
//  std::cout << "Info = " << info << std::endl;
//  dgetrs_(&trans, &mat_Size, &nrhs, (double *)A, &LDA, ipiv, (double *)B, &LDB, &info);
  dgesv_(&mat_Size, &mat_Size, (double *)A, &LDA, ipiv, (double *)B, &LDB, &info);
  std::cout << "solution is:";
  printVector(B,mat_Size);
  std::cout << "Info = " << info << std::endl;

  for(int i=0; i<mat_Size; i++)
  {
  //  delete[] A[i];
  }
  //delete [] A;
  //free(u);

  return(0);
}
