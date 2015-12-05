int init_Matrix(double **A ,int SIZE);
void freeMatrix(double **A, int SIZE);
void init_Vector(double *A, int SIZE);

bool generateMatrix(double **A, int SIZE);
bool generateVector(double *Res, int SIZE);
bool generatePreconditioner(double **A, int SIZE);

void printMatrix(double **A, int mat_Size);
void printVector(double *A, int mat_Size);
void printVectorMat(double *B, int mat_Size);

double vectorDot(double *r, double *rT, int vec_Size);
void mat_vector_mult(double **mat, double *vec, int edge_Size, double **result);
bool checkConvergence(double *vec, int mat_Size);
