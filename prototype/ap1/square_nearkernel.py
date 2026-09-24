"""Struktur der kleinsten (gewichteten) Singulaervektoren der Kollokationsmatrix:
Gradverteilung, Lokalisierung an den Ecken, Konvergenz mit N."""
import numpy as np
from square_bem import square_mesh, E_matrix, T1_matrix
def analyse(N, kappa, beta, nvec=4, a=None, b=None):
    P0, P1, Nr, V = square_mesh(N); E, X, L = E_matrix(P0, P1, Nr); T = T1_matrix(E, Nr, kappa, a, b)
    r = np.min([np.linalg.norm(X - v, axis=1) for v in V], axis=0)
    w = np.repeat(r**beta*np.sqrt(L), 8)
    A = w[:, None]*T/w[None, :]
    U, s, Vh = np.linalg.svd(A)
    out = []
    for k in range(1, nvec+1):
        g = Vh[-k].conj().reshape(-1, 8)               # gewichtete Koordinaten (L^2-normiert)
        e = np.abs(g)**2; tot = e.sum()
        aux = e[:, [0, 7]].sum()/tot
        loc = e[r < 0.1].sum()/tot
        out.append((s[-k], aux, loc, g, X))
    return out
if __name__ == '__main__':
    for kappa, name in [(-3+0.3j, 'innen'), (-4+0j, 'verlustfrei -4'), (-11+1.2j, 'Gold')]:
        for beta in [0.0, 0.35]:
            for N in [32, 64]:
                res = analyse(N, kappa, beta)
                print(f"{name:14s} beta={beta} N={N}: " + "; ".join(f"s={s:.4f} Hilfsgrade {aux:.2f} Ecke(r<0,1) {loc:.2f}" for s, aux, loc, _, _ in res[:4]), flush=True)
