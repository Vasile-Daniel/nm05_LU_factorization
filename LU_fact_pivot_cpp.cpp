#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

class LUSolver {
public:
    int n;
    vector<vector<double>> A;
    vector<double> b;
    vector<vector<double>> L;
    vector<vector<double>> U;
    vector<int> pivot;

    LUSolver(int size) : n(size), A(size, vector<double>(size, 0)), b(size, 0),
                         L(size, vector<double>(size, 0)), U(size, vector<double>(size, 0)),
                         pivot(size) {}

    void setMatrix(const vector<vector<double>>& matrix) {
        A = matrix;
    }

    void setVector(const vector<double>& vector) {
        b = vector;
    }

    void printMatrix(const vector<vector<double>>& matrix) const {
        for (const auto& row : matrix) {
            for (double val : row) {
                cout << setw(10) << fixed << setprecision(4) << val << " ";
            }
            cout << endl;
        }
    }

    void printVector(const vector<double>& vector) const {
        for (int i = 0; i < vector.size(); i++) {
            cout << "b[" << i << "] = " << setw(10) << fixed << setprecision(4) << vector[i] << endl;
        }
    }

    void decompose() {
        for (int i = 0; i < n; i++) {
            pivot[i] = i;
        }

        for (int i = 0; i < n; i++) {
            int maxRow = i;
            double maxVal = fabs(A[i][i]);

            for (int k = i + 1; k < n; k++) {
                if (fabs(A[k][i]) > maxVal) {
                    maxVal = fabs(A[k][i]);
                    maxRow = k;
                }
            }

            if (maxRow != i) {
                swap(A[i], A[maxRow]);
                swap(pivot[i], pivot[maxRow]);
            }

            for (int j = i; j < n; j++) {
                U[i][j] = A[i][j] - dotProduct(L[i], U, i, j);
            }

            for (int j = i; j < n; j++) {
                if (i == j) {
                    L[i][i] = 1;
                } else {
                    L[j][i] = (A[j][i] - dotProduct(L[j], U, i, i)) / U[i][i];
                }
            }
        }
    }

    vector<double> forwardSubstitution() {
        vector<double> y(n, 0);
        for (int i = 0; i < n; i++) {
            y[i] = b[pivot[i]] - dotProduct(L[i], y, i);
        }
        return y;
    }

    vector<double> backSubstitution(const vector<double>& y) {
        vector<double> x(n, 0);
        for (int i = n - 1; i >= 0; i--) {
            x[i] = (y[i] - dotProduct(U[i], x, i + 1, n)) / U[i][i];
        }
        return x;
    }

    vector<double> solve() {
        decompose();

        cout << "L matrix:" << endl;
        printMatrix(L);

        cout << "U matrix:" << endl;
        printMatrix(U);

        vector<double> y = forwardSubstitution();

        cout << "Vector y after forward substitution:" << endl;
        printVector(y);

        vector<double> x = backSubstitution(y);

        cout << "The solution to the system of equations is:" << endl;
        for (int i = 0; i < n; i++) {
            cout << "x[" << i << "] = " << fixed << setprecision(4) << x[i] << endl;
        }
        return x;
    }

private:
    double dotProduct(const vector<double>& vec, const vector<vector<double>>& mat, int row, int col) const {
        double result = 0;
        for (int k = 0; k < row; k++) {
            result += vec[k] * mat[k][col];
        }
        return result;
    }

    double dotProduct(const vector<double>& vec1, const vector<double>& vec2, int start) const {
        double result = 0;
        for (int k = 0; k < start; k++) {
            result += vec1[k] * vec2[k];
        }
        return result;
    }
};

int main() {
    int n = 3;
    LUSolver solver(n);

    vector<vector<double>> a = {
        {1, 2, -1},
        {-2, 3, 1},
        {4, -1, -3}
    };
    vector<double> bValues = {-1, 0, -2};

    solver.setMatrix(a);
    solver.setVector(bValues);

    cout << "Initial matrix:" << endl;
    solver.printMatrix(solver.A);

    cout << "Initial vector b:" << endl;
    solver.printVector(solver.b);

    solver.solve();

    return 0;
}
