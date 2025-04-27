import numpy as np
import matplotlib
import os
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def plot_graph(matrix, title, filename):
    plt.figure(figsize=(6, 5))
    vmax = np.max(np.abs(matrix))
    plt.imshow(matrix, cmap='seismic', vmin=-vmax, vmax=vmax)
    plt.colorbar()
    plt.title(title)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

n, m = 5, 4
np.random.seed(42)
A = np.random.randn(n, m)

# 1. A
plot_graph(A, 'Macierz A', 'A')
print("Macierz A:\n", A)

# 2. AA^T
AAT = A @ A.T
plot_graph(AAT, 'Macierz AA^T', 'AAT')
print("Macierz AA^T:\n", AAT)

# 3. Eigenvalues & Eigenvectors AA^T
eigvals_U, eigvecs_U = np.linalg.eigh(AAT)
idx = np.argsort(eigvals_U)[::-1]
eigvals_U = eigvals_U[idx]
eigvecs_U = eigvecs_U[:, idx]

print("Wartości własne AA^T:\n", eigvals_U)
print("Wektory własne U:\n", eigvecs_U)

# 4. U i S
S = np.diag(np.sqrt(np.maximum(eigvals_U, 0)))
plot_graph(eigvecs_U, 'Macierz U', 'U')
plot_graph(S, 'Macierz S', 'S')
print("Macierz S:\n", S)

# 5. V
S_inv = np.zeros_like(S)
for i in range(min(n, m)):
    if S[i, i] > 1e-10:
        S_inv[i, i] = 1.0 / S[i, i]

V = A.T @ eigvecs_U @ S_inv
plot_graph(V.T, 'Macierz V^T', 'VT')
print("Macierz V^T:\n", V.T)

# 7. ATA
ATA = A.T @ A
plot_graph(ATA, 'Macierz A^T A', 'ATA')
print("Macierz A^T A:\n", ATA)

# 8. Eigenvalues & Eigenvectors A^T A
eigvals_V, eigvecs_V = np.linalg.eigh(ATA)
idx2 = np.argsort(eigvals_V)[::-1]
eigvals_V = eigvals_V[idx2]
eigvecs_V = eigvecs_V[:, idx2]

print("Wartości własne A^T A:\n", eigvals_V)
print("Wektory własne V:\n", eigvecs_V)

# 9. S (z ATA)
S2 = np.diag(np.sqrt(np.maximum(eigvals_V, 0)))
plot_graph(eigvecs_V.T, 'Macierz V^T (z A^T A)', 'V_z_ATA')
plot_graph(S2, 'Macierz S (z A^T A)', 'S_z_ATA')
print("Macierz S:\n", S2)

# 10. U with AVS^-1
S2_inv = np.zeros_like(S2)
for i in range(min(n, m)):
    if S2[i, i] > 1e-10:
        S2_inv[i, i] = 1.0 / S2[i, i]

U2 = A @ eigvecs_V @ S2_inv
plot_graph(U2, 'Macierz U (z AVS^{-1})', 'U_z_AVS')
print("Macierz U:\n", U2)

# 12. Comparison
reconstructed_A1 = eigvecs_U @ S @ V.T
reconstructed_A2 = U2 @ S2 @ eigvecs_V.T

print("Błąd rekonstrukcji metodą pierwszą:", np.linalg.norm(A - reconstructed_A1))
print("Błąd rekonstrukcji metodą drugą:", np.linalg.norm(A - reconstructed_A2))

# 13. Rank & Nullity
rank_A = np.linalg.matrix_rank(A)
nullity_A = m - rank_A

print(f"dim R(A) = {rank_A}")
print(f"dim N(A) = {nullity_A}")
