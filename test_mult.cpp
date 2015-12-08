#include <iostream>
#include "vector.h"
#include <cmath>

using namespace std;

int npes, rank;
int krylov_solver(struct sparse *A, double *rK, int mat_Size);


void Init_MPI(int *argc, char ***argv)
{
  int err;
  err = MPI_Init(argc, argv);
  err = MPI_Comm_size(MPI_COMM_WORLD, &npes);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}

void partition(unsigned long int* load_, unsigned long int* start_, int Ni)
{
  int load =0, start =0;
	int load_division = 1;
#ifndef SCALING
	load_division = npes;
#endif

	load = (Ni/load_division);
	start = (Ni/load_division)*rank;
	if(rank < (Ni)%load_division)
	{
		load++;	
		start = start + rank;
	}
	else if(rank == Ni%npes)
	{
		start = start + rank;
	}
	else
	{
		start = start + Ni%npes;
	}
  *start_ = start;
  *load_ = load;
}

int main(int argc, char **argv)
{
  std::cout.precision(15);


  int N =10;
  long unsigned int load, start;
  Init_MPI(&argc, &argv);

  int SIZE=N-2, mat_Size= SIZE*SIZE;
  partition(&load, &start, mat_Size);
  cout<<rank<<":"<<start<<" <--start point, load --> "<<load<<endl<<" PEs"<<npes<<endl;

  double **A = new double*[mat_Size];
  double *B = new double[load];
  struct sparse A_sparse;

  A_sparse.values= new double[load*mat_Size];
  A_sparse.col_ind = new int[load*mat_Size];
  A_sparse.row_ptr = new int[load+1];
  A_sparse.nnz = mat_Size;


  init_Matrix(A, mat_Size);
  init_Vector(B, mat_Size);

  generateMatrix(A, SIZE);
  generateVector(B, SIZE, start, start+load);
  generateSparse(&A_sparse, start, start+load);

  //printMatrix(A, mat_Size);
  //printVector(B, load);
 // printSparse(&A_sparse, load);

  //cout<<endl;
  krylov_solver(&A_sparse, B, mat_Size);
  //cout<<endl<<endl;
  // printVector(B,mat_Size);
  return 0;

  delete[] A_sparse.values;
  delete[] A_sparse.col_ind;
  delete[] A_sparse.row_ptr;
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
  if(p_k==NULL || result_k==NULL || a_norm==NULL)
    return -1;

#ifndef PRECOND
  //double **M_precondition = new double*[mat_Size];    //declare and initialize preconditioning matrix
  struct sparse M_precondition;
  M_precondition.values= new double[mat_Size*mat_Size];
  M_precondition.col_ind = new int[mat_Size*mat_Size];
  M_precondition.row_ptr = new int[mat_Size+1];
  M_precondition.nnz = mat_Size;

  //init_Matrix(M_precondition, mat_Size);
  //generatePreconditioner(M_precondition, mat_Size);
  generateSparsePreconditioner(&M_precondition);

  double *z_k = new double[mat_Size];
  memset(z_k, 0, mat_Size*sizeof(double));

  //mat_vector_mult(M_precondition, rK, mat_Size, &z_k);
  sparse_mat_vector(&M_precondition, rK, mat_Size, &z_k);
  printSparse(&M_precondition, mat_Size);

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
    // mat_vector_mult(M_precondition, rK, mat_Size, &z_k);   //Stores the result of the matrix vector multiplication
    sparse_mat_vector(&M_precondition, rK, mat_Size, &z_k);
    beta_k = vectorDot(z_k, rK, mat_Size)/rK_dot;
    for(int i=0; i<mat_Size; i++)
      p_k[i] = z_k[i] + beta_k*p_k[i];
    is_converge = checkConvergence(rK, mat_Size);
  }
  cout<<k<<" iter\n";
  delete[] M_precondition.values;
  delete[] M_precondition.col_ind;
  delete[] M_precondition.row_ptr;
#else
  int k;

  memcpy(p_k, rK, mat_Size*sizeof(double));   //Stores minimizer for kth step pk

  for(k =0; k<100 && !is_converge; k++)
  {
    memset(a_norm, 0, mat_Size*sizeof(double));
    //mat_vector_mult(A, p_k, mat_Size, &a_norm);   //Stores the result of the matrix vector multiplication
    sparse_mat_vector(A, p_k, mat_Size, &a_norm);

    double rK_dot= vectorDot(rK, rK, mat_Size);    //Stores the value of the dot product of the residual
    cout<<"result "<<rK_dot<<endl;
    return 0;

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

  printVectorMat(result_k, mat_Size);
  cout<<endl;
  delete[] p_k;
  delete[] result_k;
  delete[] a_norm;
  return -1;
}
