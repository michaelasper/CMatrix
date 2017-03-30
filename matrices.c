#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct Matrix {
    int     rowSize;
    int     columnSize;
    long int*    matrix;
} Matrix;

Matrix createMatrix(int r, int c){
    Matrix temp = {r, c, calloc(r * c, sizeof(long int *))};
    return temp;
}

void printMatrix(Matrix *m){
    int i,j;
    for(i = 0; i < m->rowSize ; i++){
        for(j = 0; j < m->columnSize; j++){
            printf("%li ", *(m->matrix + i*m->rowSize  + j));
        }
        printf("\n");
    }
}

Matrix add(Matrix *a, Matrix *b){
    //check if matrices are compatible
    if(a->rowSize != b->rowSize || a->columnSize != b->columnSize){
        fprintf(stderr, "Error: Incompatible sizes");
        exit(0);
    }
    //create result matrix
    int r = a->rowSize;
    int c = a->columnSize;
    Matrix result = createMatrix(r,c);
    //add matrix
    int i,j;
    for(i = 0; i < r ; i++){
        for(j = 0; j < c; j++){
            *(result.matrix+ i*r  + j) = *(a->matrix + i*r  + j) + *(b->matrix + i*r  + j);
        }
    }
    return result;
}

Matrix sub(Matrix *a, Matrix *b){
    //check if matrices are compatible
    if(a->rowSize != b->rowSize || a->columnSize != b->columnSize){
        fprintf(stderr, "Error: Incompatible sizes");
        exit(0);
    }
    //create result matrix
    int r = a->rowSize;
    int c = a->columnSize;
    Matrix result = createMatrix(r,c);
    //add matrix
    int i,j;
    for(i = 0; i < r ; i++){
        for(j = 0; j < c; j++){
            *(result.matrix+ i*r  + j) = *(a->matrix + i*r  + j) - *(b->matrix + i*r  + j);
        }
    }
    return result;
}

Matrix multiply(Matrix *a, Matrix *b){
    if(a->columnSize != b->rowSize ){
        fprintf(stderr, "Error: Incompatible sizes");
        exit(0);
    }
    int r = a->rowSize;
    int c = b->columnSize;
    Matrix result = createMatrix(r,c);
    int i,j;
    for(i = 0; i < r ; i++){
        for(j = 0; j < c; j++){
            long int sum = 0;
            int k;
            for(k = 0; k < a->columnSize; k++){
                sum = sum + (*(a->matrix + i*a->rowSize  + k)**(b->matrix + k*b->rowSize  + j));
            }
            *(result.matrix+ i*r  + j) = sum;
        }
    }
    return result;
}

void randomize(Matrix *m){
    int i,j;
    for(i = 0; i < m->rowSize ; i++){
        for(j = 0; j < m->columnSize; j++){
            *(m->matrix + i*m->rowSize  + j)= rand() % 5000;
        }
    }
}

int main(){
    time_t t;
    srand((unsigned) time(&t));
    Matrix a = createMatrix(3,100);
    Matrix b = createMatrix(100,3);
    randomize(&a);
    randomize(&b);
    Matrix result = multiply(&a,&b);
    printMatrix(&result);
    return 0;
}