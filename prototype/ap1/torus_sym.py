"""Hauptwinkel auf rotationssymmetrischen Geometrien (Torus) mit Zerlegung in
azimutale Moden der zyklischen Gruppe C_n (Rotation um e3, kombiniert mit der
Rotor-Wirkung R h(R^-1 x) R~ auf die Multivektorwerte)."""
import numpy as np
from fields import gp, blade
from torus import torus_pts
from mfs import Sphere

def rotor_matrix(D):
    R = np.cos(D/2)*blade(0) - np.sin(D/2)*blade(3)       # exp(-e12 D/2)
    Rr = np.cos(D/2)*blade(0) + np.sin(D/2)*blade(3)      # Reversion
    M = np.zeros((8, 8), complex)
    for b in range(8): M[:, b] = gp(gp(R, blade(b)), Rr)
    return M.real

class TorusSym(Sphere):
    def __init__(self, R=1.0, r=0.4, nt=32, L=2, n=32, mt=12, d_in=0.2, d_out=0.25):
        self.n, self.L, self.nt = n, L, nt
        self.X, self.N, self.W = torus_pts(R, r, nt, n*L)
        # Basisquellen: nur eine Orbit-Repraesentantin je Tubuswinkel (phi = 0)
        t = 2*np.pi*(np.arange(mt)+0.5)/mt
        mk = lambda rr: np.stack([R+rr*np.cos(t), 0*t, rr*np.sin(t)], -1)
        self.src_in, self.src_out = mk(r-d_in), mk(r+d_out)
        self.rho = rotor_matrix(2*np.pi/n)
    def mode_sectors(self, H):
        """H: (Q, 8, m) Spuren der Basisspalten auf dem vollen Gitter ->
        Liste ueber Moden m von Sektorvektoren (nt*L*8, m), gewichtet."""
        nt, L, n = self.nt, self.L, self.n
        Hg = H.reshape(nt, n*L, 8, -1)
        W = self.W.reshape(nt, n*L)[:, :L]
        out = []
        rhop = [np.linalg.matrix_power(self.rho, q) for q in range(n)]
        # Terme T_q[t, l] = rho^q v[t, l - L q]
        Tq = [np.einsum('ij,tljm->tlim', rhop[q], np.roll(Hg, L*q, axis=1)[:, :L]) for q in range(n)]
        Tq = np.stack(Tq, 0)                           # (n, nt, L, 8, m)
        F = np.fft.fft(Tq, axis=0)                     # sum_q e^{-2 pi i m q/n} T_q
        for m in range(n):
            S = F[m] * np.sqrt(W)[:, :, None, None]
            out.append(S.reshape(-1, S.shape[-1]))
        return out

def orthq(A, tol=1e-11):
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    return U[:, s > tol*s[0]]

def angles_sym(geo, kind, om, e1, m1, e2, m2, a=None, b=None, dual=False, W='J', nang=4):
    k1, k2 = om*np.sqrt(e1)*np.sqrt(m1), om*np.sqrt(e2)*np.sqrt(m2)
    re_, rm = e1/e2, m1/m2
    if a == 'energy': a, b = np.conj(re_)*np.sqrt(rm), np.conj(rm)*np.sqrt(re_)
    if a is None: a = np.sqrt(rm)
    if b is None: b = np.sqrt(re_)
    if not dual:
        A = geo.J(geo.traces(kind, k1, e1, 'int'), re_, rm, a, b, inverse=True); B = geo.traces(kind, k2, e2, 'ext')
    else:
        A = geo.traces(kind, k2, e2, 'int'); B = geo.traces(kind, k1, e1, 'ext')
        if W == 'J': B = geo.J(B, re_, rm, a, b, inverse=True)
    res = []
    for m, (Am, Bm) in enumerate(zip(geo.mode_sectors(A), geo.mode_sectors(B))):
        c = np.linalg.svd(orthq(Am).conj().T @ orthq(Bm), compute_uv=False)
        res += [(np.sqrt(max(0, 1-min(ci, 1)**2)), m) for ci in c[:nang]]
    return sorted(res)[:nang]

if __name__ == '__main__':
    import time
    from torus import Torus, angles
    # Validierung gegen die direkte Rechnung (gleiche Gitter)
    g = Torus(nt=16, nphi=32, mt=8, mphi=16); gs = TorusSym(nt=16, L=2, n=16, mt=8)
    for kind in ['dirac', 'maxwell']:
        print(kind, 'direkt:', np.round(angles(g, kind, 0.5, 2.25, 1, 1.0, 1), 5),
              ' Moden:', [(round(v, 5), m) for v, m in angles_sym(gs, kind, 0.5, 2.25, 1, 1.0, 1)], flush=True)
