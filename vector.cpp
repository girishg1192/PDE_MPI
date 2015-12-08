#include<stdlib.h>
#include<cmath>
#include<iostream>
#include "vector.h"

float delta = 1.0;
int boundary=0;
double size = 0;
#define MIN (1e-15)

using namespace std;

double boundarySum(int X, int Y);

int init_Matrix(double **A ,int SIZE)
{
  boundary= (int)sqrt(SIZE)+1;
  size = 1/(float)SIZE;
  for(int i =0; i<SIZE; i++)
  {
    //A[i] = (double *)malloc(SIZE*sizeof(double));
    A[i] = (double *) new double[SIZE];
    if(A[i] == NULL)
      return -1;
  }

  for(int i=0; i<SIZE; i++)
  {
    for(int j=0; j<SIZE; A[i][j]=0.0,j++);
  }
  return 0;
}
void freeMatrix(double **A, int SIZE)
{
  for(int i =0; i<SIZE; i++)
  {
    double *free_ptr = A[i];
  }
}
void init_Vector(double *A, int SIZE)
{
  for(int i=0; i<SIZE; A[i]=0.0,i++);
}
bool generateMatrix(double **A, int SIZE)
{
  if(A == NULL)
    return false;
  for(int i=0; i<SIZE ; i++)
  {
    for(int j=0; j<SIZE; j++)
    {
      if(i == j)    // Filling Close diagonals
      {
        for(int k=0; k<SIZE; k++)
        {
          int off = i*SIZE;
          A[k + off][k + off] = -4;
          if(k+1 < SIZE)
            A[k + off][k+1 + off] = 1;
          if(k-1 >= 0)
            A[k + off][k-1 + off] = 1;
        }
      }
      else if(j == (i-1) || j == (i+1))    //Filing far diagonal with delta^2
      {
        int off_X = i*SIZE, off_Y = j*SIZE;
        for(int k=0 ; k<SIZE ; k++)
          A[k + off_X][k + off_Y] = 1;
      }
    }
  }
  return true;
}

void generateSparse(struct sparse *A)
{
  double *val = A->values;
  int *col_ind = A->col_ind;
  int *row_ptr = A->row_ptr;
  int SIZE = A->nnz;
  //int val[SIZE*SIZE], row_ptr[SIZE+1], col_ind[SIZE*SIZE];
  int N = sqrt(SIZE);
  unsigned int nnz=0;   //NonZero integers
  unsigned int vector_sizes = SIZE;
  row_ptr[0] = 0;
  for(int i=0; i<SIZE; i++)
  {
    if((i-N)>=0)
    {
      val[nnz] = 1.0;
      col_ind[nnz++] = i-N;
    }
    if((i-1)>=0 && i%N!=0)
    {
      val[nnz] = 1;
      col_ind[nnz++] = i-1;
    }
    val[nnz] = -4;
    col_ind[nnz++] = i;
    if((i+1)>=0 && (i+1)%N!=0)
    {
      val[nnz] = 1;
      col_ind[nnz++] = i+1;
    }
    if(i+N<(SIZE))
    {
      val[nnz] = 1;
      col_ind[nnz++] = i+N;
    }
    row_ptr[i+1] = nnz;
  }
  A->nnz = nnz;
}

bool generateVector(double *Res, int SIZE)
{
  int X,Y;
  double val;
  for(int k = 0; k< (SIZE*SIZE); k++)
  {
    X = k%SIZE+1;
    Y = k/SIZE + 1;
    Res[k] = boundarySum(X,Y);
  }
  return false;
}
bool generatePreconditioner(double **A, int SIZE)
{
  if(A==NULL)
    return false;
  for(int i=0; i<SIZE; i++)
    A[i][i] = -0.25;
  return true;
}
bool generateSparsePreconditioner(struct sparse *A)
{
  double *val = A->values;
  int *col_ind = A->col_ind;
  int *row_ptr = A->row_ptr;
  int SIZE = A->nnz;
  int nnz=0;
  row_ptr[0] = 0;
  for(int k=0; k<(SIZE*SIZE); k++)
  {
    val[k] = -0.25;
    col_ind[k] = k;
    row_ptr[k+1] = ++nnz;
  }
  return true;
}

bool isBoundary(int X, int Y)
{
  if(X == 0 || X == boundary || Y ==0 || Y == boundary )
    return true;
  return false;
}

double boundaryConditions(int X, int Y)
{
  double pt_X = ((double)X)/(boundary);
  double pt_Y = ((double)Y)/(boundary);
  if(X == 0 || X == boundary)
    return 0.0;
  if(Y == 0)
    return sin(M_PI*pt_X);
  if(Y == boundary)
    return sin(M_PI*pt_X)*exp(-M_PI);
  return 0.0;
}

double boundarySum(int X, int Y)
{
  double sum=0.0;
  for(int i=-1; i<=1; i=i+2)
    if(isBoundary(X+i, Y))
    {
      sum -= (delta*delta)*boundaryConditions(X+i, Y);
    }
  for(int i=-1; i<=1; i=i+2)
    if(isBoundary(X, Y+i))
    {
      sum -= (delta)*(delta)*boundaryConditions(X, Y+i);
    }
  return sum;
}
void printMatrix(double **A, int mat_Size)
{
  for(int i=0; i<mat_Size; i++)
  {
    for(int j=0; j<mat_Size; j++)
      cout<<A[i][j]<<" ";
    cout<<endl;
  }
}

void printVector(double *B, int vec_Size)
{
  for(int i=0; i<vec_Size; i++)
  {
    cout<<B[i]<<" ";
    cout<<endl;
  }
}
void printVectorMat(double *B, int vec_Size)
{
  for(int i=0; i<vec_Size; i++)
  {
    if(i % (int)sqrt(vec_Size) == 0)
      cout<<endl;
    cout<<B[i]<<" ";
  } 
  cout<<endl;
}
void printSparse(struct sparse *A, int mat_Size)
{
  cout<<A->nnz<<endl;
  cout<<"Values \n";
  for(int i=0; i<A->nnz; i++)
    cout<<A->values[i]<<" ";
  cout<<endl<<"Col: \n";
  for(int i=0; i<A->nnz; i++)
    cout<<A->col_ind[i]<<" ";
  cout<<endl<<"Row pointers\n";
  for(int i=0; i<(mat_Size+1); i++)
    cout<<A->row_ptr[i]<<" ";
  cout<<endl;
}
double vectorDot(double *r, double *rT, int vec_Size)
{
  double result = 0.0;
  for(int i=0; i<vec_Size; i++)
  {
    result +=r[i]*rT[i];
  }
  //cout<<"Dot result "<<result<<endl;
  return result;
}
void mat_vector_mult(double **mat, double *vec, int edge_Size, double **result)
{
  for(int i=0; i<edge_Size; i++)
    for(int j=0; j<edge_Size; j++)
      *(*result + i)+= mat[i][j] * vec[j];
}
void sparse_mat_vector(struct sparse *A, double *vec, int edge_Size, double **result)
{
  for(int i=0; i<edge_Size; i++)
    for(int j=A->row_ptr[i]; j<A->row_ptr[i+1]; j++)
      *(*result + i)+= A->values[j]*vec[A->col_ind[j]];
}
bool checkConvergence(double *vec, int vec_Size)
{
  int min_ctr=0;
  for(int i=0; i<vec_Size; i++)
    if(abs(vec[i]) > MIN) return false;  //if any value is greater than minimum
  return true;
}
