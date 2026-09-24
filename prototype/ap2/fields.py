"""Vektorisierte Feldbausteine fuer die MFS-Experimente (Teil II).
Multivektoren als Arrays (..., 8) in der Bitmasken-Basis von cl3.py."""
import numpy as np
from cl3 import SIGN

# Multiplikation vektorisiert: C[..., a^b] += SIGN[a,b] A[...,a] B[...,b]
_IDX = [(a, b, a ^ b, SIGN[a, b]) for a in range(8) for b in range(8)]
def gp(A, B):
    A, B = np.broadcast_arrays(A, B)
    C = np.zeros(A.shape, dtype=complex)
    for a, b, c, s in _IDX:
        C[..., c] += s * A[..., a] * B[..., b]
    return C

def vec(v):                     # (...,3) -> (...,8)
    M = np.zeros(v.shape[:-1] + (8,), dtype=complex)
    M[..., 1], M[..., 2], M[..., 4] = v[..., 0], v[..., 1], v[..., 2]
    return M

def blade(b):
    M = np.zeros(8, dtype=complex); M[b] = 1; return M

PSEUDO = blade(7)

def dirac_kernel(x, xs, k):
    """Phi_k(x - xs) als Multivektor, x: (N,3), xs: (3,)"""
    z = x - xs; r = np.linalg.norm(z, axis=-1)[..., None]; zh = z / r
    Phi = np.exp(1j * k * r) / (4 * np.pi * r)
    M = vec(Phi * (1 / r - 1j * k) * zh)
    M[..., 0] += (-1j * k * Phi)[..., 0]
    return M

def maxwell_dipole(x, xs, p, k, eps):
    """F = sqrt(eps)E + I sqrt(mu)H fuer elektrischen Dipol p bei xs (Normierung beliebig).
    E = curl curl(Phi p),  sqrt(mu) H = -i k sqrt(eps) grad(Phi) x p."""
    z = x - xs; r = np.linalg.norm(z, axis=-1)[..., None]; zh = z / r
    Phi = np.exp(1j * k * r) / (4 * np.pi * r)
    pr = np.sum(zh * p, axis=-1)[..., None]
    E = Phi * ((k**2 + 1j * k / r - 1 / r**2) * p + (-k**2 - 3j * k / r + 3 / r**2) * pr * zh)
    gPhi = Phi * (1j * k - 1 / r) * zh
    sH = -1j * k * np.sqrt(eps) * np.cross(gPhi, p)
    return np.sqrt(eps) * vec(E) + gp(PSEUDO, vec(sH))

def dirac_residual(F, x, k, h=1e-4):
    """(nabla - ik)F per zentraler Differenz; F: callable (N,3)->(N,8)."""
    R = -1j * k * F(x)
    for j, b in enumerate([1, 2, 4]):
        dx = np.zeros(3); dx[j] = h
        R = R + gp(blade(b), (F(x + dx) - F(x - dx)) / (2 * h))
    return R
