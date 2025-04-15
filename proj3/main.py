import numpy as np
import matplotlib.pyplot as plt
import time

# === Manual implementations ===
def norm_1(M):
    return max(sum(abs(M[i][j]) for i in range(len(M))) for j in range(len(M[0])))

def norm_inf(M):
    return max(sum(abs(M[i][j]) for j in range(len(M[0]))) for i in range(len(M)))

def norm_frobenius(M):
    return (sum(M[i][j]**2 for i in range(len(M)) for j in range(len(M[0]))))**0.5

def norm_2(M):
    MtM = np.dot(np.transpose(M), M)
    eigenvalues = np.linalg.eigvals(MtM)
    max_eigenvalue = max(abs(eig) for eig in eigenvalues)
    return max_eigenvalue ** 0.5

def condition_number(norm_func, M, M_inv):
    norm_M = norm_func(M)
    norm_M_inv = norm_func(M_inv)
    return norm_M * norm_M_inv

# === Testing and Plotting===
def test_implementation(matrix_size=1000, seed=42, epsilon=1e-6):
    np.random.seed(seed)
    M = np.random.randn(matrix_size, matrix_size)

    M_list = M.tolist()
    M_inv = np.linalg.inv(M)
    M_inv_list = M_inv.tolist()

    # --- Our implementations ---
    for name, func, norm_builtin in [
        ("Norm 1", norm_1, lambda M: np.linalg.norm(M, 1)),
        ("Norm ∞", norm_inf, lambda M: np.linalg.norm(M, np.inf)),
        ("Frobenius", norm_frobenius, lambda M: np.linalg.norm(M, 'fro')),
        ("Norm 2", norm_2, lambda M: np.linalg.norm(M, 2)),
    ]:
        val_manual = func(M_list)

        val_np = norm_builtin(M)

        assert abs(val_manual - val_np) < epsilon, f"Error for implementation {name}: {abs(val_manual - val_np)} < {epsilon}"

    # --- Condition numbers ---
    for name, norm_func, norm_builtin in [
        ("Cond 1", norm_1, lambda M: np.linalg.cond(M, 1)),
        ("Cond ∞", norm_inf, lambda M: np.linalg.cond(M, np.inf)),
        ("Cond Frobenius", norm_frobenius, lambda M: np.linalg.norm(M, 'fro') * np.linalg.norm(np.linalg.inv(M), 'fro')),
        ("Cond 2", norm_2, lambda M: np.linalg.cond(M, 2)),
    ]:
        cond_manual = condition_number(norm_func, M_list, M_inv_list)

        cond_np = norm_builtin(M)

        assert abs(cond_manual - cond_np) < epsilon, f"Error for condition in {name}: {abs(cond_manual - cond_np)} < {epsilon}"

def benchmark_and_plot(matrix_size=1000, seed=42):
    np.random.seed(seed)
    M = np.random.randn(matrix_size, matrix_size)
    M_list = M.tolist()

    norm_tests = [
        ("Norm 1", norm_1, lambda A: np.linalg.norm(A, 1)),
        ("Norm ∞", norm_inf, lambda A: np.linalg.norm(A, np.inf)),
        ("Norm p (Frobenius)", norm_frobenius, lambda A: np.linalg.norm(A, 'fro')),
        ("Norm 2", norm_2, lambda A: np.linalg.norm(A, 2)),
    ]

    names, time_manuals, time_np = [], [], []

    for name, manual_func, np_func in norm_tests:
        names.append(name)

        t0 = time.time()
        manual_func(M_list)
        time_manuals.append(time.time() - t0)

        t0 = time.time()
        np_func(M)
        time_np.append(time.time() - t0)

    # === Plotting ===
    x = np.arange(len(names))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x - width/2, time_manuals, width, label='Nasza implementacja')
    ax.bar(x + width/2, time_np, width, label='NumPy')

    ax.set_ylabel('Czas wykonania (s)')
    ax.set_title('Porównanie czasu wykonania: Nasza implementacja vs NumPy')
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.legend()
    plt.tight_layout()
    plt.show()

test_implementation(1000)
benchmark_and_plot(1000)