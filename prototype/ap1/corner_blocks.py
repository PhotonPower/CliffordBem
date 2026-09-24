"""Analytische Struktur des Eck-Mellin-Symbols: Zerlegung in vier 2x2-Bloecke je Strahl
({1,e12}, {e1,e2}, {e3,e123}, {e13,e23}) und Faktorisierung der Determinante in skalare
Faktoren g(c;a) = S(a)^2 - ((1+c)/(1-c))^2 mit blockspezifischem Kontrast c."""
import numpy as np
from corner import M_E, Jfull
from corner_weights import S
BL = {'TM-E {e1,e2}': [1, 2], 'Hilfs A {1,e12}': [0, 3], 'Hilfs B {e3,e123}': [4, 7], 'TE-H {e13,e23}': [5, 6]}
def T1(alpha, tau, re_, rm, a, b, c=0.5):
    M = M_E(alpha, tau, c); I = np.eye(16)
    return 0.5*(I + M) + 0.5*(I - M)@Jfull(alpha, re_, rm, a, b)
def idx(bl): return [bl[0], bl[1], 8 + bl[0], 8 + bl[1]]
def g(cc, alpha, am): q = (1 + cc)/(1 - cc); return S(alpha, am)**2 - q**2
if __name__ == '__main__':
    alpha = np.pi/2; rng = np.random.default_rng(1)
    re_, rm = -3+0.3j, 1.7+0.2j; a, b = 0.6-0.4j, 1.3+0.5j
    Tm = T1(alpha, 0.37, re_, rm, a, b)
    allidx = sum((idx(v) for v in BL.values()), []); off = Tm.copy()
    for v in BL.values():
        ii = idx(v); off[np.ix_(ii, ii)] = 0
    print("Nebenblock-Anteil (sollte 0 sein):", np.abs(off).max())
    # Kandidaten fuer den Kontrast jedes Blocks
    se, sm = np.sqrt(re_), np.sqrt(rm)
    cands = {'r_eps': re_, 'r_mu': rm, 'a/sqrt(r_mu)': a/sm, 'sqrt(r_mu)/a': sm/a, 'b/sqrt(r_eps)': b/se,
             'sqrt(r_eps)/b': se/b, '1/r_eps': 1/re_, '1/r_mu': 1/rm}
    for name, v in BL.items():
        ii = idx(v); ratios = {}
        for cn, cv in cands.items():
            vals = []
            for tau in [0.0, 0.37, 1.1, 2.3]:
                am = 0.5 + 1j*tau
                d = np.linalg.det(T1(alpha, tau, re_, rm, a, b)[np.ix_(ii, ii)])
                vals.append(d/g(cv, alpha, am))
            vals = np.array(vals); ratios[cn] = np.abs(vals/vals[0] - 1).max()
        best = min(ratios, key=ratios.get)
        print(f"{name:20s}: bester Kontrast {best:14s} (Abw. des Quotienten {ratios[best]:.1e}; naechstbester {sorted(ratios.values())[1]:.1e})")
