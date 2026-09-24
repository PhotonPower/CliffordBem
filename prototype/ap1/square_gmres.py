"""GMRES-Iterationszahlen fuer T_1 h = h_inc auf geometrisch gradierten Netzen:
ungewichtet (L^2) gegen gewichtetes Skalarprodukt (Aehnlichkeit mit r^beta).
h_inc: homogenes Feld E0 = e1 (quasistatisch, Spur einer inneren Loesung)."""
import numpy as np, json, time
from scipy.sparse.linalg import gmres
from square_geo import build_geo
def run(ng, kappa, beta_norm, tol=1e-8):
    A, X, rr = build_geo(ng, 0.0, kappa)           # A = D T D^{-1}, D = sqrt(L) (L^2-Koordinaten)
    L = build_geo.last['L']; d = np.repeat(np.sqrt(L), 8)
    h = np.zeros((len(X), 8), complex); h[:, 1] = 1.0; f = d*h.ravel()
    w = rr**beta_norm; Aw = w[:, None]*A/w[None, :]; fw = w*f
    it = [0]
    def cb(res): it[0] += 1
    x, info = gmres(Aw, fw, rtol=tol, restart=len(fw), maxiter=1, callback=cb, callback_type='pr_norm')
    res = np.linalg.norm(Aw@x - fw)/np.linalg.norm(fw)
    return it[0], res, info
if __name__ == '__main__':
    out = {}
    for ng in [8, 16, 24, 32]:
        row = []
        for name, k in [('innen', -3+0.3j), ('-4', -4+0j), ('Gold', -11+1.2j), ('Glas', 2.25+0j)]:
            for b in [0.0, 0.35]:
                n, r, info = run(ng, k, b); out[f'{name}|{b}|{ng}'] = n
                row.append(f'{name} b={b}: {n}')
        print(f'ng={ng}: ' + ', '.join(row), flush=True)
    json.dump(out, open('square_gmres.json', 'w'))
