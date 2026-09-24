"""Eck-Mellin-Symbol mit chiralem Innenmedium: Kopplung von TM- und TE-Block ueber die
Normalbedingungen; Bestimmung der Determinante als Polynom in S(a)^2."""
import numpy as np
from corner import M_E, geometry
from chiral_J import Jchiral
from corner_weights import S
def T1_chiral(alpha, tau, m1, m2=(1.0, 1.0, 0.0), c=0.5):
    M = M_E(alpha, tau, c); I = np.eye(16); _, n = geometry(alpha)
    J = np.zeros((16, 16), complex)
    for i in range(2): J[8*i:8*i+8, 8*i:8*i+8] = Jchiral(n[i], m1, m2)
    return 0.5*(I + M) + 0.5*(I - M)@J
TMTE = [1, 2, 5, 6, 9, 10, 13, 14]
def coupled_det(alpha, tau, m1):
    T = T1_chiral(alpha, tau, m1); return np.linalg.det(T[np.ix_(TMTE, TMTE)]), T
if __name__ == '__main__':
    alpha = np.pi/2
    for m1 in [(2.25, 1.0, 0.0), (-3+0.3j, 1.0, 0.4), (-2+0.1j, 1.3, 0.25+0.02j)]:
        e, u, ch = m1
        taus = np.array([0.0, 0.3, 0.7, 1.2, 2.0]); Ss = np.array([S(alpha, 0.5+1j*t) for t in taus])
        ds = []
        for t in taus:
            d, T = coupled_det(alpha, t, m1)
            off = T.copy(); off[np.ix_(TMTE, TMTE)] = 0
            aux = [0, 3, 4, 7, 8, 11, 12, 15]; off[np.ix_(aux, aux)] = 0
            ds.append(d)
        ds = np.array(ds)
        V = np.stack([np.ones_like(Ss), Ss**2, Ss**4], 1)
        coef, *_ = np.linalg.lstsq(V[:3], ds[:3], rcond=None)
        err = np.abs(V@coef - ds).max()/np.abs(ds).max()
        # Vermutung: (1/16 e mu') * [((1+e)(1+u) - ch^2) - ... ]  -> Koeffizienten normieren
        cn = coef/coef[0]
        print(f"m1={m1}: det_TM+TE = c0 + c1 S^2 + c2 S^4, Rest {err:.1e}; c1/c0 = {cn[1]:.6f}, c2/c0 = {cn[2]:.6f}; c0 = {coef[0]:.6f}")
        print(f"     Kopplung Hilfsbloecke/Physik: max Nebenblock {np.abs(off).max():.1e}")
