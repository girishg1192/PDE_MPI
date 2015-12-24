#include<cstring>
#include "mpi.h"

/* Routines to initialize memory for vectors/Matrices*/
int init_Matrix(double **A ,int SIZE);
void freeMatrix(double **A, int SIZE);
void init_Vector(double *A, int SIZE);

struct sparse;
/*Routine for generating initial Matrices
 * generateMatrix : Generates coefficient matrix which is SIZE*SIZE
 * generateVector: Generates result matrix
 * generatePreconditioner: Generates jacobi preconditioner
 */
bool generateMatrix(double **A, int SIZE);
void generateSparse(struct sparse *A, int start, int end);
bool generateVector(double *Res, int SIZE, int start, int end);
bool generatePreconditioner(double **A, int SIZE);
bool generateSparsePreconditioner(struct sparse *A, int start, int end);

/*Print routines for debugging*/
void printMatrix(double **A, int mat_Size);
void printVector(double *A, int mat_Size);
void printVectorMat(double *B, int mat_Size);
void printSparse(struct sparse *A, int mat_Size);

/*Matrix/Vector computation functions
 * vectorDot: returns dot product of vectors r and rT (r.rT)
 * mat_vector_mult: computes matrix vector product
 * mat*vec, Result value stored in vector result
 */
double vectorDot(double *r, double *rT, int vec_Size);
void mat_vector_mult(double **mat, double *vec, int edge_Size, double **result);
void sparse_mat_vector(struct sparse *A, double *vec, double **result);

/* checks convergence of the algorithm given the residual vector vec
 */
bool checkConvergence(double *vec, int mat_Size);

struct sparse
{
  double *values;
  int *col_ind;
  int *row_ptr;
  int nnz;
};
