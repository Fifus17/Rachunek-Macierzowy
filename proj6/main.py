import numpy as np


def givens_rotation(a, b):
    r = np.hypot(a, b)
    if r == 0:
        c = 1
        s = 0
    else:
        c = a / r
        s = -b / r
    return c, s

def qr_givens(A):
    (m, n) = A.shape
    Q = np.eye(m)
    R = A.copy()

    for j in range(n):
        for i in range(m - 1, j, -1):
            if R[i, j] != 0:
                c, s = givens_rotation(R[i - 1, j], R[i, j])
                G = np.eye(m)
                G[[i - 1, i], [i - 1, i]] = c
                G[i, i - 1] = s
                G[i - 1, i] = -s

                R = G @ R
                Q = Q @ G.T 
    return Q, R

def test_qr_givens():
    np.random.seed(42)
    A = np.random.randn(32, 32)
    Q, R = qr_givens(A)

    reconstruction_error = np.linalg.norm(A - Q @ R)

    orthogonality_error = np.linalg.norm(Q.T @ Q - np.eye(32))

    Q_np, R_np = np.linalg.qr(A)
    equivalent_to_numpy = np.allclose(Q @ R, Q_np @ R_np)

    print("Błąd rekonstrukcji ||A - QR|| =", reconstruction_error)
    print("Błąd ortogonalności ||QᵀQ - I|| =", orthogonality_error)
    print("Zgodność z NumPy QR? ", equivalent_to_numpy)

    return reconstruction_error, orthogonality_error, equivalent_to_numpy

test_qr_givens()