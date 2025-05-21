import numpy as np
import matplotlib.pyplot as plt

def power_method(A, p=2, epsilon=1e-4, max_iter=15, z=None):
    if z is None:
        z = np.random.rand(3)
    z = z / np.linalg.norm(z, ord=p)

    errors = []
    for _ in range(max_iter):
        w = A @ z
        λ = np.max(w)
        z = w / np.linalg.norm(w, ord=p)
        error = np.linalg.norm(A @ z - λ * z, ord=p)
        errors.append(error)
        if error < epsilon:
            break
    return errors


# Losowa macierz A
A = np.random.rand(3, 3)
print(A)

# SVD ręczne
AAT = A @ A.T
eigvals, U = np.linalg.eigh(AAT)
D = np.diag(np.sqrt(eigvals[::-1]))
U = U[:, ::-1]
V = A.T @ U @ np.linalg.inv(D)
S = U @ D @ V.T

print(f"A: {A}")
print(f"U: {U}")
print(f"D: {D}")
print(f"V: {V}")

# Błędy dla ||UDV - A||_p (z poprawką spłaszczania macierzy)
errors_svd = {
    p: np.linalg.norm((S - A).flatten(), ord=p) for p in [1, 2, 3, 4, np.inf]
}
print("Błędy ||UDV - A||_p:", errors_svd)

initial_vectors = [np.random.rand(3) for _ in range(3)]

# Loop over norms
for p in [1, 2, 3, 4, np.inf]:
    for z in initial_vectors:
        errors = power_method(A, p=p, z=z)
        label = f'z={np.round(z, 2)}'
        plt.plot(errors, label=label)
    plt.yscale('log')
    plt.xlabel('Iteracja')
    plt.ylabel(f'Błąd (norma p={p})')
    plt.title(f'Zbieżność metody potęgowej dla normy p={p}')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'zbieznosc_norma_p{p}.png')
    plt.clf()