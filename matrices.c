/**
    Matrix Multiplication
    matrices.c
    Matrix data structure in C.

    @author Michael Asper
    @version 1.0 3/29/17
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

typedef struct Matrix {
    double              scalar;
    size_t const        numRows;
    size_t const        numCols;
    double**            matrix;
} Matrix;


/**
    Randomizes the elements of a matrix

    @param *m pointer to Matrix to randomize;
*/
void randomize(Matrix *m){
    for(size_t i = 0; i < m->numRows ; i++){
        for(size_t j = 0; j < m->numCols; j++){
             m->matrix[i][j] = rand() % 10;
        }
    }
}


/**
    Returns a r x c Matrix with all 0s.

    @param r The row size of the matrix
    @param c The column size of the matrix
    @return r x c Matrix
*/
Matrix createMatrix(size_t const r, size_t const c){
    Matrix temp = {1, r, c, calloc(r, sizeof(double))};

    if (temp.matrix == NULL) {
        fprintf(stderr, "Error: Memory allocation failed");
        exit(1);
    }

    for (size_t i = 0; i < r; i++) {
        temp.matrix[i] = calloc(c, sizeof temp.matrix[i][0]);

        if (temp.matrix[i] == NULL) {
             fprintf(stderr, "Error: Memory allocation failed");
             exit(1);
        }
    }

    return temp;
}

/**
    Copy Matrix

    @param m Matrix M
    @return copy of A
*/
Matrix copy(Matrix m){
    Matrix result = createMatrix(m.numRows, m.numCols);
    for(size_t i = 0; i < m.numRows ; i++){
        for(size_t j = 0; j < m.numCols; j++){
            result.matrix[i][j] = m.matrix[i][j];
        }
    }
    return result;
}

/**
    Returns a r x c Matrix with random numbers.

    @param r The row size of the matrix
    @param c The column size of the matrix
    @return r x c Matrix
*/
Matrix createRandMatrix(size_t const r, size_t const c){
    Matrix temp = createMatrix(r,c);
    randomize(&temp);
    return temp;
}

/**
    Prints matrix.

    @param *m Pointer to Matrix you want to print
*/
void printMatrix(Matrix m){
    
    for(size_t i = 0; i < m.numRows ; i++){
        for(size_t j = 0; j < m.numCols; j++){
            printf("%-5.2f ",  m.matrix[i][j]);
        }
        printf("\n");
    }
}

/**

    Multipy matrix by scalar;

    @param s scalar
    @param m matrix (M)
    @return s*M
*/
Matrix scalarMultiply(int s, Matrix m){
    Matrix result = createMatrix(m.numRows,m.numCols);
    for(size_t i = 0; i < m.numRows ; i++){
        for(size_t j = 0; j < m.numCols; j++){
            result.matrix[i][j]=m.matrix[i][j] * s;
        }
    }
    return result;
}

/**
    Adds two matrices together

    @param a first matrix (A);
    @param b second matrix (B);
    @return A+B
*/
Matrix add(Matrix a, Matrix b){
    //check if matrices are compatible
    if(a.numRows != b.numRows || a.numCols != b.numCols){
        fprintf(stderr, "Error: Incompatible sizes");
        exit(0);
    }
    //create result matrix
    size_t r = a.numRows;
    size_t c = a.numCols;
    Matrix result = createMatrix(r,c);
    //add matrices
    for(size_t i = 0; i < r ; i++){
        for(size_t j = 0; j < c; j++){
            result.matrix[i][j] = a.matrix[i][j]+b.matrix[i][j];
        }
    }
    return result;
}

/**
    Subtracts two matrices together

    @param a first matrix (A);
    @param b second matrix (B);
    @return A-B
*/
Matrix sub(Matrix a, Matrix b){
    return add(a,scalarMultiply(-1,b));
}

/**
    Multiplies two matrices together

    @param a first matrix (A);
    @param b second matrix (B);
    @return A*B
*/
Matrix multiply(Matrix a, Matrix b){
    //check if matrices are compatible
    if(a.numCols != b.numRows ){
        fprintf(stderr, "Error: Incompatible sizes");
        exit(0);
    }
    
    //initialize return matrix
    size_t r = a.numRows;
    size_t c = b.numCols;
    Matrix result = createMatrix(r,c);
    
    //multiply matrices
    for(size_t i = 0; i < r ; i++){
        for(size_t j = 0; j < c; j++){
            long int sum = 0;
            for(size_t k = 0; k < a.numCols; k++){
                sum += a.matrix[i][k] * b.matrix[k][j];
            }
            result.matrix[i][j] = sum;
        }
    }
    return result;
}


/**
    Swap rows in a Matrix

    @param r1 first row
    @param r2 second row
    @param m Matrix
*/
void swapRows(size_t r1, size_t r2, Matrix *m){
    for (size_t i=0; i<m->numCols; ++i){
        double temp = m->matrix[r2][i];
        m->matrix[r2][i] = m->matrix[r1][i];
        m->matrix[r1][i] = temp;
    }
}

/**

    Row reduces a matrix using Gaussian elimination with partial pivoting.

    @param m matrix (M)
*/
Matrix reduce(Matrix m){
    Matrix result = copy(m);
    for(size_t i = 0; (i < result.numRows && i < result.numCols); i++){

        int max = fabs(result.matrix[i][i]);
        int maxRow = i;
        
        for(size_t j = i; j < result.numRows; j++){
            if (max < fabs(result.matrix[j][i])){
                max = fabs(result.matrix[j][i]);
                maxRow = j;
            }
        }

        if(maxRow != i){
            swapRows(i, maxRow, &result);
            result.scalar *=-1; 
        }
        
        for(size_t j = i+1; j < result.numRows; j++){
            double factor = result.matrix[j][i]/result.matrix[i][i];
            for(size_t k = 0; k < result.numCols; k++){
                result.matrix[j][k] = result.matrix[j][k] - factor*result.matrix[i][k];
            }
        }
    }

    return result;
}

/**
    Computes determinant for a square Matrix

    @param *m pointer to matrix (A)
    @return double determinant of matrix (det(A))
*/
double calcDet(Matrix m){
    if(m.numCols != m.numRows){
        fprintf(stderr, "Error: only square matrices have determinants");
        exit(2);
    }
    
    Matrix temp = reduce(m);
    double result = 1;
    for(size_t i = 0; (i < temp.numCols && i < temp.numRows); ++i){
        result *= temp.matrix[i][i];
    }
    return temp.scalar * result;
}

int main(){
    // seed random with time
    time_t t;
    srand((unsigned) time(&t));

    //setup random matrices and multiply
    Matrix a = createRandMatrix(3,3);
    printMatrix(a);
    printf("\n\n");
    a = reduce(a);
    printf("Determinant: %.2f", calcDet(a));
    
    return 0;
}