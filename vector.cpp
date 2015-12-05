#include<stdlib.h>
#include<cmath>
#include<iostream>
float delta = 1.0;
int boundary=0;
double size = 0;

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

void printVector(double *B, int mat_Size)
{
  for(int i=0; i<mat_Size; i++)
  {
    cout<<B[i]<<" ";
    cout<<endl;
  }
}
void printVectorMat(double *B, int mat_Size)
{
  for(int i=0; i<mat_Size; i++)
  {
    cout<<B[i]<<" ";
    if(i % (int)sqrt(mat_Size) == 0)
      cout<<endl;
  } 
}
