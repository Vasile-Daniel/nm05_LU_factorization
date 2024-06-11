#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// Author: Vasile-Daniel DAN 
// Start
// Project: LU Factorization - Doolittle (WITH pivoting)
//////////////////////////////////////////////////////////////////////////////////  
// Functions for printing matrices and arrays
void printMatrix(int n, double matrix[n][n]) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%10.4lf ", matrix[i][j]);
        }
        printf("\n");
    }
}

void printArray(int n, double arr[n]){
    int i;
    for (i = 0; i < n; i++) {
        printf("b[%d] = %10.4lf \n", i, arr[i]);
    }
}

// Function to swap two rows in a matrix
void swapRows(int n, double matrix[n][n], int row1, int row2) {
    double temp;
    for (int i = 0; i < n; i++) {
        temp = matrix[row1][i];
        matrix[row1][i] = matrix[row2][i];
        matrix[row2][i] = temp;
    }
}

// Function to swap two elements in an array
void swapElements(double arr[], int i, int j) {
    double temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;
}

// LU decomposition function (Doolittle variant) with pivoting
void luDecomposition(int n, double a[n][n], double l[n][n], double u[n][n], int pivot[n]) {
    int i, j, k;
    for (i = 0; i < n; i++) {
        pivot[i] = i;
    }

    for (i = 0; i < n; i++) {
        // Partial pivoting
        double max = fabs(a[i][i]);
        int maxRow = i;
        for (k = i + 1; k < n; k++) {
            if (fabs(a[k][i]) > max) {
                max = fabs(a[k][i]);
                maxRow = k;
            }
        }

        if (maxRow != i) {
            swapRows(n, a, i, maxRow);
            swapElements(pivot, i, maxRow);
        }

        // Initialize U
        for (j = i; j < n; j++) {
            u[i][j] = a[i][j];
            for (k = 0; k < i; k++) {
                u[i][j] -= l[i][k] * u[k][j];
            }
        }

        // Initialize L
        for (j = i; j < n; j++) {
            if (i == j) {
                l[i][i] = 1;
            } else {
                l[j][i] = a[j][i];
                for (k = 0; k < i; k++) {
                    l[j][i] -= l[j][k] * u[k][i];
                }
                l[j][i] /= u[i][i];
            }
        }
    }
}

// Forward substitution to solve Ly = Pb
void forwardSubstitution(int n, double l[n][n], double b[n], double y[n], int pivot[n]) {
    for (int i = 0; i < n; i++) {
        y[i] = b[pivot[i]];
        for (int j = 0; j < i; j++) {
            y[i] -= l[i][j] * y[j];
        }
    }
}

// Backward substitution to solve Ux = y
void backSubstitution(int n, double u[n][n], double y[n], double x[n]) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= u[i][j] * x[j];
        }
        x[i] /= u[i][i];
    }
}

int main() {
    int n = 3;
    double x[3], y[3];
    int pivot[3];

    double a[3][3] = {
        { 1, 2, -1},
        {-2, 3, 1},
        { 4, -1, -3}
    };

    printf("Initial matrix:\n");
    printMatrix(n, a);

    double b[3] = {-1, 0, -2};

    printf("Initial vector b:\n");
    printArray(n, b);

    double l[3][3] = {0}, u[3][3] = {0};

    luDecomposition(n, a, l, u, pivot);

    printf("L matrix:\n");
    printMatrix(n, l);
    printf("U matrix:\n");
    printMatrix(n, u);

    forwardSubstitution(n, l, b, y, pivot);

    printf("Vector y after forward substitution:\n");
    printArray(n, y);

    backSubstitution(n, u, y, x);

    printf("The solution to the system of equations is:\n");
    for(int i = 0; i < n; i++){
        printf("x[%d] = %lf \n", i, x[i]);
    }

    return 0;
}
