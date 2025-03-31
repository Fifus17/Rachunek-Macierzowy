import numpy as np
import time


def print_matrix(A):
    for row in A:
        print(" ".join(map(str, row)))


def check_solution(A, b, x, tol=1e-6):
    b_computed = np.dot(A, x)

    if np.allclose(b_computed, b, atol=tol):
        return True
    else:
        print("Solution is wrong")
        print("Good b:", b)
        print("Yours b:", b_computed)
        return False


def backward_substitution(A):
    n = len(A)
    result = [0] * n

    for i in range(n-1, -1, -1):
        result[i] = (A[i, -1] - np.dot(A[i, i + 1:n], result[i + 1:n])) / A[i, i]

    return result


def forward_substitution(A):
    n = len(A)
    for i in range(n):
        A[i, -1] -= np.dot(A[i, :i], A[:i, -1])
        A[i, -1] /= A[i, i]
    return A[:, -1]


def gauss_ones_no_pivot(A):
    n = len(A)

    for i in range(n):
        divide_by = A[i][i]
        for j in range(i, n+1):
            A[i][j] /= divide_by
        for j in range(i+1, n):
            factor = A[j][i]
            for k in range(i, n+1):
                A[j][k] -= factor * A[i][k]

    return backward_substitution(A)


def gauss_pivot(A):
    n = len(A)

    for i in range(n):
        max_row = np.argmax(abs(A[i:n, i])) + i
        A[[i, max_row]] = A[[max_row, i]]

        divide_by = A[i][i]
        for j in range(i, n+1):
            A[i][j] /= divide_by
        for j in range(i+1, n):
            factor = A[j][i]
            for k in range(i, n+1):
                A[j][k] -= factor * A[i][k]

    return backward_substitution(A)


def lu_decomposition_no_pivot(A, b):
    n = len(A)

    L = np.eye(n)
    U = A.copy()

    for i in range(n):
        for j in range(i + 1, n):
            factor = U[j, i] / U[i, i]
            L[j, i] = factor
            U[j, i:] -= factor * U[i, i:]

    L_expanded = np.hstack((L, b.reshape(-1, 1)))
    y = forward_substitution(L_expanded)

    U_expanded = np.hstack((U, y.reshape(-1, 1)))
    x = backward_substitution(U_expanded)

    return x


def lu_decomposition_pivot(A, b):
    n = len(A)

    L = np.eye(n)
    U = A.copy()

    for i in range(n):
        max_row = np.argmax(abs(U[i:n, i])) + i
        if i != max_row:
            U[[i, max_row], :] = U[[max_row, i], :]
            L[[i, max_row], :i] = L[[max_row, i], :i]
            b[i], b[max_row] = b[max_row], b[i]

        for j in range(i + 1, n):
            factor = U[j, i] / U[i, i]
            L[j, i] = factor
            U[j, i:] -= factor * U[i, i:]

    L_expanded = np.hstack((L, b.reshape(-1, 1)))
    y = forward_substitution(L_expanded)

    U_expanded = np.hstack((U, y.reshape(-1, 1)))
    x = backward_substitution(U_expanded)

    return x


if __name__ == "__main__":
    n = 42
    A = np.random.rand(n, n)
    b = np.random.rand(n)
    A_expanded = np.column_stack((A, b))

    start_time = time.time()
    result_gauss_no_pivot = gauss_ones_no_pivot(A_expanded)
    time_gauss_no_pivot = time.time() - start_time
    print(time_gauss_no_pivot)

    start_time = time.time()
    result_gauss_pivot = gauss_pivot(A_expanded)
    time_gauss_pivot = time.time() - start_time

    start_time = time.time()
    result_lu_no_pivot = lu_decomposition_no_pivot(A, b)
    time_lu_no_pivot = time.time() - start_time

    start_time = time.time()
    result_lu_pivot = lu_decomposition_no_pivot(A, b)
    time_lu_pivot = time.time() - start_time

    print("---CHECKING SOLUTIONS---")
    print(check_solution(A, b, result_gauss_no_pivot))
    print(check_solution(A, b, result_gauss_pivot))
    print(check_solution(A, b, result_lu_no_pivot))
    print(check_solution(A, b, result_lu_pivot))

    print("---TIME COMPARISON---")
    print(f"Gauss no pivot: {time_gauss_no_pivot:.4f}s")
    print(f"Gauss pivot: {time_gauss_pivot:.4f}s")
    print(f"Lu no pivot: {time_lu_no_pivot:.4f}s")
    print(f"Lu pivot: {time_lu_pivot:.4f}s")
