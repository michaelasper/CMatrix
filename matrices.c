/**
    Matrix Multiplication
    matrices.c
    Matrix data structure in C.

    @author Michael Asper
    @version 1.0 2017-30-3
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>

typedef struct Matrix {
    double complex              scalar;
    size_t const                numRows;
    size_t const                numCols;
    double complex**            matrix;
} Matrix;


double complex sign(double complex a){
    return a / cabs(a);
}

/**
    Randomizes the elements of a matrix

    @param *m pointer to Matrix to randomize;
*/
void randomize(Matrix *m){
    for(size_t i = 0; i < m->numRows ; i++){
        for(size_t j = 0; j < m->numCols; j++){
             m->matrix[i][j] = ((double)rand()/RAND_MAX * 5) + (((double)rand()/RAND_MAX * 5) * I);
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
    Matrix temp = {1, r, c, calloc(r, sizeof(double[c]))};

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
    Returns a r x r Matrix

    @param r The row x row size of the matrix
    @return r x r Matrix
*/
Matrix createIdentity(size_t const r){
    Matrix temp = createMatrix(r,r);

    for(size_t i = 0; i < r; i++){
        temp.matrix[i][i] = 1;
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
    Delete a matrix
*/
void delete(Matrix m)
{
    for (size_t i = m.numRows; i > 0; i--) {
        free(m.matrix[i-1]);
    }
	free(m.matrix);
}
 

/**
    Prints matrix.

    @param *m Pointer to Matrix you want to print
*/
void printMatrix(Matrix m){
    
    for(size_t i = 0; i < m.numRows ; i++){
        for(size_t j = 0; j < m.numCols; j++){
            printf("(%5.3f + %5.3fi) ", crealf(m.matrix[i][j]), cimagf(m.matrix[i][j]));
            //printf("(%5.3f) ", crealf(m.matrix[i][j]));
        }
        printf("\n");
    }
            printf("\n");

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
            double complex sum = 0;
            for(size_t k = 0; k < a.numCols; k++){
                sum += a.matrix[i][k] * b.matrix[k][j];
            }
            result.matrix[i][j] = sum;
        }
    }
    return result;
}

/**
    Matrix exponentiation for n >= 0
    "Exponentiation by squaring"
    Note: This algorithm runs in O(logn) time 
    where n x n is the size of the square matrix
   
    Not sure what to do for n < 0
    
    @param a Matrix to be exponentiated
    @param n exponent
    @return A^n
*/
Matrix power(Matrix a, int n) {
	Matrix b = multiply(a, a);
	if(n == 0) {
	    return createIdentity(a.numRows);
	}
	else if(n == 1){
	    return a;
	}
	else if(n % 2 == 1) {
	    return multiply(a, power(b, (n-1)/2));
	}
	else {
	    return power(b, n/2);
	}
}


/**
    Swap rows in a Matrix

    @param r1 first row
    @param r2 second row
    @param m Matrix
*/
void swapRows(size_t r1, size_t r2, Matrix *m){
    for (size_t i=0; i<m->numCols; ++i){
        double complex temp = m->matrix[r2][i];
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

        double max = cabs(result.matrix[i][i]);
        int maxRow = i;
        
        for(size_t j = i; j < result.numRows; j++){
            if (max < cabs(result.matrix[j][i])){
                max = cabs(result.matrix[j][i]);
                maxRow = j;
            }
        }

        if(maxRow != i){
            swapRows(i, maxRow, &result);
            result.scalar *=-1; 
        }
        for(size_t j = i+1; j < result.numRows; j++){
            double complex factor = result.matrix[j][i]/result.matrix[i][i];
            for(size_t k = 0; k < result.numCols; k++){
                result.matrix[j][k] = result.matrix[j][k] - factor*result.matrix[i][k];
            }
        }
    }

    return result;
}

/**
    Transpose a Matrix

    @param m Matrix(A)
    @return A^t
*/
Matrix transpose(Matrix m){
    Matrix result = createMatrix(m.numCols,m.numRows);

    for(size_t i = 0; i < result.numRows; ++i){
        for(size_t j = 0; j < result.numCols; ++j){
            result.matrix[i][j] = conj(m.matrix[j][i]);
        }
    }

    return result;
}


/**
    Computes determinant for a square Matrix

    @param *m pointer to matrix (A)
    @return double determinant of matrix (det(A))
*/
double complex calcDet(Matrix m){
    if(m.numCols != m.numRows){
        fprintf(stderr, "Error: only square matrices have determinants");
        exit(2);
    }

    Matrix temp = reduce(m);
    double complex result = 1;
    for(size_t i = 0; (i < temp.numCols && i < temp.numRows); ++i){
        result *= temp.matrix[i][i];
    }
    return temp.scalar * result;
}




double norm(double complex a, double complex b){
    return sqrt(cabs(a)*cabs(a) + cabs(b)*cabs(b));
}

Matrix givensrotation(Matrix m, size_t i, size_t j){
    double complex a = m.matrix[i-1][j];
    double complex b = m.matrix[i][j];
    double c;
    double complex s;
    double abs_a = cabs(a);
    
    if (abs_a == 0){
        c = 0;
        s = 1;
    }else{
        double nrm = norm(a,b);
        c = abs_a/nrm;
        s = (a/abs_a) * (conj(b)/nrm);
    }
    
    Matrix givens = createIdentity(m.numRows);
    givens.matrix[i][i-1] = -conj(s);
    givens.matrix[i-1][i-1] = c;
    givens.matrix[i-1][i] = s;
    givens.matrix[i][i] = c;
    return givens;
}

/**
    QR factorization using Given's rotation

    Dear lord, please have mercy on me
*/
void eigenvalues(Matrix m){
    Matrix R = copy(m);
    Matrix givens;
    Matrix Q = createIdentity(m.numRows);
    for(int g = 0; g < 20; g++){
        givens = createIdentity(m.numRows);
        for(size_t j = 0; j < (m.numCols-1);j++){
            for(size_t i = m.numRows-1; i > (0+j); --i){
                Matrix temp = givensrotation(R, i, j);

                R = multiply(temp, R);
                givens = multiply(givens,transpose(temp));
                delete(temp);
            }
        }
        R = multiply(R, givens);
    }
    printf("Eigenvalues: ");
    for(size_t i = 0;i< m.numRows;i++){
        if(i != m.numRows-1)
            printf("%5.3f + %5.3fi, ", crealf(R.matrix[i][i]), cimagf(R.matrix[i][i]));
        else
            printf("%5.3f + %5.3fi", crealf(R.matrix[i][i]), cimagf(R.matrix[i][i]));

    }
}

int main(){
    // seed random with time
    time_t t;
    srand((unsigned) time(&t));

    //setup random matrices and multiply
    Matrix a = createRandMatrix(3,3);
    eigenvalues(a);
    return 0;
    
}
