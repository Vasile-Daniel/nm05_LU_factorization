import numpy as np

class LUSolver:
    def __init__(self, size):
        self.n = size
        self.A = np.zeros((size, size))
        self.b = np.zeros(size)
        self.L = np.zeros((size, size))
        self.U = np.zeros((size, size))
        self.pivot = np.arange(size)

    def set_matrix(self, matrix):
        self.A = matrix

    def set_vector(self, vector):
        self.b = vector

    def print_matrix(self, mat):
        for row in mat:
            print(" ".join(f"{val:10.4f}" for val in row))

    def print_vector(self, vec):
        for i, val in enumerate(vec):
            print(f"b[{i}] = {val:10.4f}")

    def decompose(self):
        n = self.n
        for i in range(n):
            maxRow = np.argmax(np.abs(self.A[i:, i])) + i
            if maxRow != i:
                self.A[[i, maxRow]] = self.A[[maxRow, i]]
                self.pivot[[i, maxRow]] = self.pivot[[maxRow, i]]

            for j in range(i, n):
                self.U[i, j] = self.A[i, j] - np.dot(self.L[i, :i], self.U[:i, j])

            for j in range(i, n):
                if i == j:
                    self.L[i, i] = 1
                else:
                    self.L[j, i] = (self.A[j, i] - np.dot(self.L[j, :i], self.U[:i, i])) / self.U[i, i]

    def forward_substitution(self):
        n = self.n
        y = np.zeros(n)
        for i in range(n):
            y[i] = self.b[self.pivot[i]] - np.dot(self.L[i, :i], y[:i])
        return y

    def back_substitution(self, y):
        n = self.n
        x = np.zeros(n)
        for i in range(n - 1, -1, -1):
            x[i] = (y[i] - np.dot(self.U[i, i + 1:], x[i + 1:])) / self.U[i, i]
        return x

    def solve(self):
        self.decompose()

        print("L matrix:")
        self.print_matrix(self.L)

        print("U matrix:")
        self.print_matrix(self.U)

        y = self.forward_substitution()

        print("Vector y after forward substitution:")
        self.print_vector(y)

        x = self.back_substitution(y)

        print("The solution to the system of equations is:")
        for i in range(self.n):
            print(f"x[{i}] = {x[i]:.4f}")
        return x

def main():
    n = 3
    solver = LUSolver(n)

    a = np.array([[1, 2, -1],
                  [-2, 3, 1],
                  [4, -1, -3]])
    b_values = np.array([-1, 0, -2])

    solver.set_matrix(a)
    solver.set_vector(b_values)

    print("Initial matrix:")
    solver.print_matrix(solver.A)

    print("Initial vector b:")
    solver.print_vector(solver.b)

    solver.solve()

if __name__ == "__main__":
    main()


"""
Explanation:
LUSolver Class:

Manages the matrix 
𝐴
A, vector 
𝑏
b, L and U matrices, and the pivot array.
Methods to set the matrix and vector (set_matrix and set_vector).
Methods to print matrices and vectors (print_matrix and print_vector).
The decompose method performs LU decomposition with pivoting.
The forward_substitution and back_substitution methods solve the system of equations.
The solve method orchestrates the solving process, including decomposition and printing intermediate steps.
Main Function:

Initializes the LUSolver with matrix 
𝐴
A and vector 
𝑏
b.
Prints the initial matrix and vector.
Calls the solve method to perform LU decomposition, solve the system, and print the solution.
This implementation consolidates all functionality within the LUSolver class, making the code more organized and easier to manage.
"""