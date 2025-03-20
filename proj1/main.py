import numpy as np
import time
import matplotlib.pyplot as plt

add_count = sub_count = mul_count = 0


def split(matrix):
    row, col = matrix.shape
    row2, col2 = row // 2, col // 2
    return matrix[:row2, :col2], matrix[:row2, col2:], matrix[row2:, :col2], matrix[row2:, col2:]


def strassen_multiplication(A, B):
    global add_count, sub_count, mul_count

    if A.shape[0] == 1:
        mul_count += 1
        return A * B

    A11, A12, A21, A22 = split(A)
    B11, B12, B21, B22 = split(B)

    M1 = strassen_multiplication(A11 + A22, B11 + B22)
    M2 = strassen_multiplication(A21 + A22, B11)
    M3 = strassen_multiplication(A11, B12 - B22)
    M4 = strassen_multiplication(A22, B21 - B11)
    M5 = strassen_multiplication(A11 + A12, B22)
    M6 = strassen_multiplication(A21 - A11, B11 + B12)
    M7 = strassen_multiplication(A12 - A22, B21 + B22)

    add_count += 6
    sub_count += 4

    C11 = M1 + M4 - M5 + M7
    C12 = M3 + M5
    C21 = M2 + M4
    C22 = M1 - M2 + M3 + M6

    add_count += 6
    sub_count += 2

    C = merge_matrices(C11, C12, C21, C22)
    return C

def merge_matrices(X_11, X_12, X_21, X_22):
    return np.vstack((np.hstack((X_11, X_12)), np.hstack((X_21, X_22))))


def matrix_multiplication(A, B, l=8):
    n = len(A)
    if n <= l:
        return strassen_multiplication(A, B)
    
    A11, A12, A21, A22 = split(A)
    B11, B12, B21, B22 = split(B)

    C11 = matrix_multiplication(A11, B11, l) + matrix_multiplication(A12, B21, l)
    C12 = matrix_multiplication(A11, B12, l) + matrix_multiplication(A12, B22, l)
    C21 = matrix_multiplication(A21, B11, l) + matrix_multiplication(A22, B21, l)
    C22 = matrix_multiplication(A21, B12, l) + matrix_multiplication(A22, B22, l)

    return merge_matrices(C11, C12, C21, C22)


# Funkcja do liczenia czasu i operacji dla różnych wartości l
def benchmark_matrix_multiplication(ks, l=8):
    global add_count, sub_count, mul_count
    results = {k: {'sizes': [], 'times': [], 'operations': []} for k in ks}
    
    for k in ks:
        print("k: ", k)
        size = 2 ** k
        np.random.seed(42)
        A = np.random.rand(size, size)
        np.random.seed(43)
        B = np.random.rand(size, size)
        
        # for l in ls:
        #     print("l: ", l)
        start_time = time.time()
        matrix_multiplication(A, B, l)
        duration = time.time() - start_time
        print(f"time: {duration}, add: {add_count}, sub: {sub_count}, mul: {mul_count}")
        results[k] = {'time': duration, 'add': add_count, 'sub': sub_count, 'mul': mul_count}
        add_count = sub_count = mul_count = 0
    
    return results

def plot_results(results):
    x = list(results.keys())
    times = [result['time'] for result in results.values()]
    adds = [result['add'] for result in results.values()]
    subs = [result['sub'] for result in results.values()]
    muls = [result['mul'] for result in results.values()]

    plt.figure(figsize=(14, 6))

    # Execution Times: Normal Scale
    plt.subplot(1, 2, 1)
    plt.plot(x, times, marker="o", linestyle="-", color="b", label="Time")
    plt.xlabel(r'Rozmiar macierzy: $2^k \times 2^k$')
    plt.ylabel("Długość obliczeń (s)")
    plt.title("Długość obliczeń - skala liniowa")
    plt.legend()
    plt.grid()

    # Execution Times: Log Scale
    plt.subplot(1, 2, 2)
    plt.yscale("log")
    plt.plot(x, times, marker="o", linestyle="-", color="b", label="Time")
    plt.xlabel(r'Rozmiar macierzy: $2^k \times 2^k$')
    plt.ylabel("Długość obliczeń (s)")
    plt.title("Długość obliczeń - skala logarytmiczna")
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)

    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(14, 6))

    # Operations Count: Normal Scale
    plt.subplot(1, 2, 1)
    plt.plot(x, adds, marker="o", linestyle="-", color="g", label="Additions")
    plt.plot(x, subs, marker="s", linestyle="--", color="r", label="Subtractions")
    plt.plot(x, muls, marker="^", linestyle="-.", color="m", label="Multiplications")
    plt.xlabel(r'Rozmiar macierzy: $2^k \times 2^k$')
    plt.ylabel("Liczba operacji")
    plt.title("Liczba operacji - skala liniowa")
    plt.legend()
    plt.grid(True, linestyle="--", linewidth=0.5)

    # Operations Count: Log Scale
    plt.subplot(1, 2, 2)
    plt.yscale("log")
    plt.plot(x, adds, marker="o", linestyle="-", color="g", label="Additions")
    plt.plot(x, subs, marker="s", linestyle="--", color="r", label="Subtractions")
    plt.plot(x, muls, marker="^", linestyle="-.", color="m", label="Multiplications")
    plt.xlabel(r'Rozmiar macierzy: $2^k \times 2^k$')
    plt.ylabel("Liczba operacji")
    plt.title("Liczba operacji - skala logarytmiczna")
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)

    plt.tight_layout()
    plt.show()


def plot_results_with_l(results, ls):
    x = list(results[ls[0]].keys())
    times = [[0 for _ in ls] for _ in results[ls[0]].items()]

    for i, result in enumerate(results.values()):
        for j, r in enumerate(result.values()):
            times[i][j] = r['time']

    plt.figure(figsize=(8, 6))

    # Execution Times: Normal Scale
    for i, l in enumerate(ls):
        plt.plot(x, times[i], marker="o", label=f"l={l}")
    plt.xlabel(r'Rozmiar macierzy: $2^k \times 2^k$')
    plt.ylabel("Długość obliczeń (s)")
    plt.title("Długość obliczeń - skala liniowa")
    plt.legend()
    plt.grid()
    plt.show()

ks = range(5, 9)  # Ile się uda wykonać
ls = [2, 3, 4, 5]  # Wybrane wartości l

# results = benchmark_matrix_multiplication(ks)

# plot_results(results)

results = {l: benchmark_matrix_multiplication(ks, l) for l in ls}

plot_results_with_l(results, ls)