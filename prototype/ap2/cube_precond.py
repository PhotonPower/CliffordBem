"""Schritt 2.5: Kanten-/Ecken-Block-Vorkonditionierung auf dem Wuerfel (je Symmetriesektor).
Bloecke im Fundamentaloktanten: Ecke (1,1,1) [Abstand < R], drei Kanten [Abstand < R],
sonst 8x8-Diagonalbloecke je Dreieck. Rechtsvorkonditioniertes GMRES."""
import numpy as np, sys, json, os, time
import scipy.linalg as sl
from scipy.sparse.linalg import gmres, LinearOperator
from cube_mesh import cube_mesh
from cube_sym import sector_matrices, fundamental, GROUP, A_of, perm_of
from assemble import project
from scatter3d import plane_wave
def groups(C, R):
    v = np.linalg.norm(C - 1, axis=1)
    ed = np.stack([np.hypot(1 - C[:, 1], 1 - C[:, 2]), np.hypot(1 - C[:, 0], 1 - C[:, 2]), np.hypot(1 - C[:, 0], 1 - C[:, 1])], 1)
    lab = np.full(len(C), -1)
    lab[np.min(ed, 1) < R] = np.argmin(ed, 1)[np.min(ed, 1) < R]
    lab[v < R] = 3
    G = [np.nonzero(lab == g)[0] for g in range(4)]
    G += [np.array([i]) for i in np.nonzero(lab == -1)[0]]
    return [g for g in G if len(g)]
def make_prec(T, G):
    lus = []
    for g in G:
        idx = np.concatenate([8*g + c for c in range(8)]); lus.append((idx, sl.lu_factor(T[np.ix_(idx, idx)])))
    def apply(v):
        y = v.copy()
        for idx, lu in lus: y[idx] = sl.lu_solve(lu, v[idx])
        return y
    return apply, max(len(idx) for idx, _ in lus)
def its(T, b, P=None):
    n = len(b); cnt = [0]
    op = LinearOperator((n, n), matvec=(lambda v: T@v) if P is None else (lambda v: T@P(v)), dtype=complex)
    x, info = gmres(op, b, rtol=1e-8, restart=n, maxiter=1, callback=lambda r: cnt.__setitem__(0, cnt[0]+1), callback_type='pr_norm')
    return cnt[0]
if __name__ == '__main__':
    L = int(sys.argv[1]); kap = complex(sys.argv[2]); om = 0.5
    RS = [float(r) for r in sys.argv[3].split(',')] if len(sys.argv) > 3 else [0.25, 0.5]
    fn = 'cube_precond.json'; res = json.load(open(fn)) if os.path.exists(fn) else {}
    key = f'{sys.argv[2]}|{L}' + ('|' + sys.argv[3] if len(sys.argv) > 3 else '')
    if key in res: print(key, res[key]); sys.exit()
    m = cube_mesh(L); F = fundamental(m); t0 = time.time()
    sec, ta, nn = sector_matrices(m, om, kap)
    b = project(m, lambda x: plane_wave(x, om, 1.0)).reshape(-1, 8); perms = {g: perm_of(m, g) for g in GROUP}
    out = {'none': 0, 'diag': 0}; out.update({f'R{R}': 0 for R in RS}); bmax = {}
    for chi, T in sec.items():
        bc = sum(np.prod([c if gi < 0 else 1 for c, gi in zip(chi, g)])*(b[perms[g][F]]@A_of(g)) for g in GROUP).ravel()
        if np.linalg.norm(bc) < 1e-12*np.linalg.norm(b): continue
        out['none'] = max(out['none'], its(T, bc))
        Pd, _ = make_prec(T, [np.array([i]) for i in range(len(F))]); out['diag'] = max(out['diag'], its(T, bc, Pd))
        for R in RS:
            P, bs = make_prec(T, groups(m.cent[F], R)); out[f'R{R}'] = max(out[f'R{R}'], its(T, bc, P)); bmax[f'R{R}'] = bs
    res[key] = dict(N=len(m.T), its=out, blockmax=bmax)
    json.dump(res, open(fn, 'w'), indent=1)
    print(f"L={L} N={len(m.T)} kappa={kap}: {out}  (groesster Block {bmax}; {time.time()-t0:.0f}s)", flush=True)
