"""Eckangepasste Vorkonditionierung: Block-Jacobi mit exakter Inversen der Eckbloecke
(alle Elemente mit Abstand < R_c zur jeweiligen Ecke), rechtsvorkonditioniertes GMRES."""
import numpy as np, json, time, sys
from scipy.sparse.linalg import gmres, LinearOperator
import scipy.linalg as sl
from square_geo import build_geo
V = np.array([(1, -1), (1, 1), (-1, 1), (-1, -1)], float)
def corner_blocks(X, rr, Rc):
    near = np.argmin([np.linalg.norm(X - v, axis=1) for v in V], axis=0)
    r = np.min([np.linalg.norm(X - v, axis=1) for v in V], axis=0)
    blocks = []
    for k in range(4):
        el = np.where((near == k) & (r < Rc))[0]
        blocks.append(np.concatenate([8*el + c for c in range(8)]))
    return blocks
def solve(ng, kappa, beta_norm=0.0, Rc=None, tol=1e-8):
    A, X, rr = build_geo(ng, 0.0, kappa)
    L = build_geo.last['L']; d = np.repeat(np.sqrt(L), 8)
    h = np.zeros((len(X), 8), complex); h[:, 1] = 1.0; f = d*h.ravel()
    w = rr**beta_norm; A = w[:, None]*A/w[None, :]; f = w*f
    n = len(f)
    if Rc is None:
        Pinv = lambda v: v
    else:
        blocks = corner_blocks(X, rr, Rc); lus = [sl.lu_factor(A[np.ix_(b, b)]) for b in blocks]
        def Pinv(v):
            y = v.copy()
            for b, lu in zip(blocks, lus): y[b] = sl.lu_solve(lu, v[b])
            return y
    Op = LinearOperator((n, n), matvec=lambda v: A@Pinv(v), dtype=complex)
    it = [0]
    x, info = gmres(Op, f, rtol=tol, restart=n, maxiter=1, callback=lambda r: it.__setitem__(0, it[0]+1), callback_type='pr_norm')
    res = np.linalg.norm(A@Pinv(x) - f)/np.linalg.norm(f)
    return it[0], res
if __name__ == '__main__':
    out = {}
    Rcs = [None, 0.1, 0.5]
    for ng in [8, 16, 24, 32]:
        row = []
        for name, k in [('innen', -3+0.3j), ('-4', -4+0j), ('Gold', -11+1.2j), ('Glas', 2.25+0j)]:
            vals = []
            for Rc in Rcs:
                for b in [0.0, 0.35]:
                    n_it, res = solve(ng, k, b, Rc); out[f'{name}|{Rc}|{b}|{ng}'] = [n_it, res]; vals.append(f'{n_it}')
            row.append(f"{name}: " + '/'.join(vals))
        print(f'ng={ng}: ' + '; '.join(row), flush=True)
    json.dump(out, open('square_precond.json', 'w'))
