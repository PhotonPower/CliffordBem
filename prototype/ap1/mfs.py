"""Subraumwinkel-Test fuer Eindeutigkeit des Transmissionsproblems auf der Kugel (Teil II)."""
import numpy as np
from fields import *
from cl3 import SIGN

def fib_sphere(n, R=1.0):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5**0.5) * i
    return R * np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], -1)

def sphere_quad(nt, nphi):
    t, wt = np.polynomial.legendre.leggauss(nt)
    ph = 2 * np.pi * np.arange(nphi) / nphi
    T, P = np.meshgrid(t, ph, indexing='ij'); W = np.repeat(wt[:, None], nphi, 1) * 2 * np.pi / nphi
    st = np.sqrt(1 - T**2)
    X = np.stack([st * np.cos(P), st * np.sin(P), T], -1).reshape(-1, 3)
    return X, W.reshape(-1)

class Sphere:
    def __init__(self, nt=36, nphi=72, nsrc=150, r_in=0.45, r_out=2.2):
        self.X, self.W = sphere_quad(nt, nphi)
        self.N = self.X.copy()                  # Normale = Ortsvektor
        self.src_in, self.src_out = fib_sphere(nsrc, r_in), fib_sphere(nsrc, r_out)
    def traces(self, kind, k, eps, where):
        """Spalten = Spuren (8*Q) der Basisloesungen. kind: 'dirac' | 'maxwell';
        where: 'int' (Quellen aussen, innere Loesungen) | 'ext'.  (vektorisiert)"""
        S = self.src_out if where == 'int' else self.src_in
        x = self.X[:, None, :]
        if kind == 'dirac':
            K = dirac_kernel(x, S[None, :, :], k)                    # (Q,S,8)
            cols = [K[..., None, :] * 0 for _ in range(0)]
            out = np.zeros(K.shape[:2] + (8, 8), dtype=complex)      # (Q,S,blade,8)
            for b in range(8):
                for a_ in range(8):
                    out[:, :, b, a_ ^ b] += SIGN[a_, b] * K[:, :, a_]
            return out.reshape(len(self.X), -1, 8).transpose(0, 2, 1)
        out = []
        for p in np.eye(3):
            F = maxwell_dipole(x, S[None, :, :], p, k, eps)          # (Q,S,8)
            out += [F, gp(PSEUDO, F)]
        return np.stack(out, 2).reshape(len(self.X), -1, 8).transpose(0, 2, 1)
    def J(self, H, re_, rm, a, b, inverse=False):
        """wendet J punktweise an; H: (Q,8,m)."""
        n = self.N
        se, sm = np.sqrt(re_), np.sqrt(rm)
        if inverse: se, sm, a, b = 1/se, 1/sm, 1/a, 1/b
        v = H[:, [1, 2, 4], :]                  # Vektor
        Ib = H[:, [6, 5, 3], :] * np.array([1, -1, 1])[None, :, None]  # I*b: e23=I e1, e13=-I e2, e12=I e3
        vn = np.einsum('qi,qim->qm', n, v); bn = np.einsum('qi,qim->qm', n, Ib)
        vT = v - n[:, :, None] * vn[:, None]; bT = Ib - n[:, :, None] * bn[:, None]
        v2 = se * vT + (1/se) * n[:, :, None] * vn[:, None]
        b2 = sm * bT + (1/sm) * n[:, :, None] * bn[:, None]
        out = np.zeros_like(H)
        out[:, 0] = a * H[:, 0]; out[:, 7] = b * H[:, 7]
        out[:, [1, 2, 4]] = v2
        out[:, [6, 5, 3]] = b2 * np.array([1, -1, 1])[None, :, None]
        return out
    def orth(self, H, tol=1e-11):
        A = (np.sqrt(self.W)[:, None, None] * H).reshape(-1, H.shape[-1])
        U, s, _ = np.linalg.svd(A, full_matrices=False)
        return U[:, s > tol * s[0]]

def min_angle(sph, kind, om, e1, m1, e2, m2, a=None, b=None, dual=False, W='J'):
    """sin des kleinsten Hauptwinkels zwischen J^{-1} ran E1^+ und ran E2^- (bzw. dual)."""
    k1, k2 = om*np.sqrt(e1)*np.sqrt(m1), om*np.sqrt(e2)*np.sqrt(m2)
    re_, rm = e1/e2, m1/m2
    if a == 'energy':      # Energiewahl (Satz 8.x): a = conj(r_eps) sqrt(r_mu), b = conj(r_mu) sqrt(r_eps)
        a, b = np.conj(re_)*np.sqrt(rm), np.conj(rm)*np.sqrt(re_)
    if a is None: a = np.sqrt(rm)
    if b is None: b = np.sqrt(re_)
    if not dual:
        A = sph.J(sph.traces(kind, k1, e1, 'int'), re_, rm, a, b, inverse=True)
        B = sph.traces(kind, k2, e2, 'ext')
    else:   # ran E2^+  cap  J^{-1} ran E1^-
        A = sph.traces(kind, k2, e2, 'int')
        B = sph.traces(kind, k1, e1, 'ext')
        if W == 'J': B = sph.J(B, re_, rm, a, b, inverse=True)
    QA, QB = sph.orth(A), sph.orth(B)
    c = np.linalg.svd(QA.conj().T @ QB, compute_uv=False)[0]
    return np.sqrt(max(0.0, 1 - min(c, 1.0)**2))
