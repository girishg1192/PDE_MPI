#include <iostream>
#include "vector.h"
#include <cmath>

using namespace std;

#ifndef SCALING
int N =200;
#else
int N =10;
#endif

/* MPI variables:
 * rank -> rank of processor
 * npes -> number of processors
 * start -> offset in grid fro where processor handles data
 * load: -> stores load to be handled by processor
 */
long unsigned int load, start;
int npes, rank;
MPI_Comm comm = MPI_COMM_WORLD;

int krylov_solver();
void partition(unsigned long int* load_, unsigned long int* start_, int Ni);


void Init_MPI(int *argc, char ***argv)
{
  int err;
  err = MPI_Init(argc, argv);
  err = MPI_Comm_size(MPI_COMM_WORLD, &npes);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}


int main(int argc, char **argv)
{
  std::cout.precision(15);

  Init_MPI(&argc, &argv);

  double time_start, time_end;

  time_start = MPI_Wtime();

  krylov_solver();

  time_end = MPI_Wtime();
 
  if(rank==0)
  cout<<"Running Time: "<<time_end-time_start<<endl;

  MPI_Finalize();
  return(0);
}

int krylov_solver()   
{
  int SIZE=N-2, mat_Size= SIZE*SIZE;
#ifndef SCALING
  partition(&load, &start, mat_Size);
#else
  load = mat_Size;
  start = mat_Size*rank;
  mat_Size = mat_Size*npes;
#endif

  double **A_ = new double*[mat_Size];
  double *rK = new double[load];
  struct sparse A;

  A.values= new double[load*5];
  A.col_ind = new int[load*5];
  A.row_ptr = new int[load+1];
  A.nnz = mat_Size;


  init_Matrix(A_, mat_Size);

  generateVector(rK, SIZE, start, start+load);    //Generates Result vector B
  generateSparse(&A, start, start+load);        // Generates Coefficient vector A

  bool is_converge = false;
  double alpha_k, beta_k;
  double *p_k = new double[load];

  double result_k[load];
  memset(result_k, 0, load*sizeof(double));

  double *a_norm = new double[load];
  memset(a_norm, 0, load*sizeof(double));
  if(p_k==NULL || result_k==NULL || a_norm==NULL)
    return -1;

#ifndef PRECOND
  //declare and initialize preconditioning matrix
  struct sparse M_precondition;
  M_precondition.values= new double[load];
  M_precondition.col_ind = new int[load];
  M_precondition.row_ptr = new int[load+1];
  M_precondition.nnz = mat_Size;

  //init_Matrix(M_precondition, mat_Size);
  //generatePreconditioner(M_precondition, mat_Size);
  generateSparsePreconditioner(&M_precondition, start, start+load);

  double *z_k = new double[load];
  memset(z_k, 0, load*sizeof(double));

  //mat_vector_mult(M_precondition, rK, mat_Size, &z_k);
  sparse_mat_vector(&M_precondition, rK, &z_k);
//  if(rank == 0)
//    printVector(z_k, load);

  memcpy(p_k, z_k, load*sizeof(double));   //Stores minimizer for kth step pk
  int k;
  for(k =0; k<5000; k++)
  {
    memset(a_norm, 0, load*sizeof(double));
    //mat_vector_mult(A, p_k, mat_Size, &a_norm);   //Stores the result of the matrix vector multiplication
    sparse_mat_vector(&A, p_k, &a_norm);

    double rK_dot = vectorDot(rK, z_k, load);    //Stores the value of the dot product of the residual

    alpha_k = rK_dot/(vectorDot(p_k, a_norm, load));
    for(int i=0; i<load; i++)
    {
      result_k[i] = result_k[i] + alpha_k*p_k[i];
      rK[i] = rK[i] - alpha_k*a_norm[i];
    }
    //Stores the result of the matrix vector multiplication
    sparse_mat_vector(&M_precondition, rK, &z_k);
    beta_k = vectorDot(z_k, rK, load)/rK_dot;
    for(int i=0; i<load; i++)
      p_k[i] = z_k[i] + beta_k*p_k[i];
//    is_converge = checkConvergence(rK, mat_Size);
  }
  delete[] M_precondition.values;
  delete[] M_precondition.col_ind;
  delete[] M_precondition.row_ptr;
#else
  int k;

  memcpy(p_k, rK, load*sizeof(double));   //Stores minimizer for kth step pk

  for(k =0; k<5000; k++)
  {
    //Stores the result of the matrix vector multiplication
    memset(a_norm, 0, load*sizeof(double));
    sparse_mat_vector(&A, p_k, &a_norm);

    double rK_dot= vectorDot(rK, rK, load);    //Stores the value of the dot product of the residual
    alpha_k = rK_dot/(vectorDot(p_k, a_norm, load));
    for(int i=0; i<load; i++)
    {
      result_k[i] = result_k[i] + alpha_k*p_k[i];
      rK[i] = rK[i] - alpha_k*a_norm[i];
    }
    beta_k = vectorDot(rK, rK, load)/rK_dot;
    for(int i=0; i<load; i++)
      p_k[i] = rK[i] + beta_k*p_k[i];
//    is_converge = checkConvergence(rK, load);
  }
#endif
 
  {
    double B_vector[mat_Size];
    init_Vector(B_vector, mat_Size);
    int rcount[npes], displs[npes];
    for(int i=0, start=0; i<npes; i++)
    {
      displs[i] = start;
      rcount[i] = mat_Size/npes;
      if(i < mat_Size%npes)
        rcount[i]++;
      start+=rcount[i];
    }
//    MPI_Allgatherv(result_k, rcount[rank], MPI_DOUBLE, &B_vector[0], rcount, displs, MPI_DOUBLE, MPI_COMM_WORLD);
//  printVectorMat(B_vector, mat_Size);
  }
  delete[] p_k;
  delete[] a_norm;
  return 0;
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
