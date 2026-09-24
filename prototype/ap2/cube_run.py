"""Stabilitaetsstudie Wuerfel: kleinster Singulaerwert (L^2 und gewichtet) und GMRES je Sektor."""
import numpy as np, sys, json, os, time
import scipy.linalg as sl
from scipy.sparse.linalg import gmres, LinearOperator
from cube_mesh import cube_mesh
from cube_sym import sector_matrices, fundamental, GROUP, A_of, perm_of
from assemble import project
from scatter3d import plane_wave
def edge_dist(X):
    a = np.abs(X); s = np.sort(a, axis=1)     # Abstand zur naechsten Kante: sqrt((1-a_i)^2 + (1-a_j)^2) der zwei groessten
    return np.sqrt((1 - s[:, 2])**2 + (1 - s[:, 1])**2)
def smin_lu(T, w, iters=40, seed=0):
    lu = sl.lu_factor(T); rng = np.random.default_rng(seed)
    x = rng.normal(size=len(T)) + 1j*rng.normal(size=len(T)); x /= np.linalg.norm(x); nrm = 0
    for _ in range(iters):
        y = w*sl.lu_solve(lu, x/w)                          # B x,   B = W T^-1 W^-1
        z = (1/w)*sl.lu_solve(lu, w*y, trans=2)             # B^H y  (B^H = W^-1 T^-H W)
        nrm = np.sqrt(np.linalg.norm(z)); x = z/np.linalg.norm(z)
    return 1/nrm
def gmres_its(T, b):
    cnt = [0]
    x, info = gmres(T, b, rtol=1e-8, restart=len(b), maxiter=1, callback=lambda r: cnt.__setitem__(0, cnt[0]+1), callback_type='pr_norm')
    return cnt[0]
if __name__ == '__main__':
    L = int(sys.argv[1]); kap = complex(sys.argv[2]); om = 0.5; beta = 0.35
    fn = 'cube_results.json'; res = json.load(open(fn)) if os.path.exists(fn) else {}
    key = f'{sys.argv[2]}|{L}'
    if key in res: print('schon vorhanden', res[key]); sys.exit()
    m = cube_mesh(L); F = fundamental(m); t0 = time.time()
    sec, ta, nn = sector_matrices(m, om, kap)
    w = np.repeat(edge_dist(m.cent[F])**beta, 8)
    # rechte Seite: ebene Welle, auf Sektoren projiziert
    b = project(m, lambda x: plane_wave(x, om, 1.0)).reshape(-1, 8); perms = {g: perm_of(m, g) for g in GROUP}
    out = {}
    for chi, T in sec.items():
        bc = sum(np.prod([c if gi < 0 else 1 for c, gi in zip(chi, g)])*(b[perms[g][F]]@A_of(g)) for g in GROUP).ravel()
        s0 = smin_lu(T, np.ones(len(T))); sw = smin_lu(T, w)
        it = gmres_its(T, bc) if np.linalg.norm(bc) > 1e-12*np.linalg.norm(b) else -1
        out[str(chi)] = [s0, sw, it]
    s0m = min(v[0] for v in out.values()); swm = min(v[1] for v in out.values()); itm = max(v[2] for v in out.values())
    res[key] = dict(N=len(m.T), smin=s0m, smin_w=swm, gmres_max=itm, sectors=out, t_asm=ta, near=int(nn))
    json.dump(res, open(fn, 'w'), indent=1)
    print(f"L={L} N={len(m.T)} kappa={kap}: sigma_min L2 {s0m:.4f}, gewichtet {swm:.4f}, GMRES max {itm}  (Assembl. {ta:.0f}s, gesamt {time.time()-t0:.0f}s, Nahpaare {nn})", flush=True)
