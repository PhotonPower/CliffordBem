"""Kombinierter Loeser am Wuerfel: H-Matrix-Operatoren (AP 3) + Kanten-/Eck-Blockvorkonditionierung
(AP 2.5) auf dem anisotropen, nichtkonformen Netz.  Aufruf: cube_h.py L kappa stage (b1|b2|solve)"""
import numpy as np, sys, pickle, json, os, time, itertools
import scipy.linalg as sl
from scipy.sparse.linalg import gmres, LinearOperator
from cube_aniso import cube_aniso
from hmat import HMatrix, Kernel
from assemble import project, mv_times_n, LB
from corner import Jblock
from scatter3d import plane_wave
from mesh import quad_points
from fields import vec, gp
L, kap, stage = int(sys.argv[1]), complex(sys.argv[2]), sys.argv[3]
om = 0.5; eps = 1e-4; k1, k2 = om*np.sqrt(kap), om
m = cube_aniso(L); N = len(m.T); tag = f'cube_{L}_{sys.argv[2]}'
if stage in ('b1', 'b2'):
    H = HMatrix(m, k1 if stage == 'b1' else k2, eps=eps, mode='joint'); del H.ker
    pickle.dump(H, open(f'{tag}_{stage}.pkl', 'wb')); sys.exit()
H1 = pickle.load(open(f'{tag}_b1.pkl', 'rb')); H2 = pickle.load(open(f'{tag}_b2.pkl', 'rb'))
J = np.array([Jblock(m.n[j], kap, 1.0, 1.0, np.sqrt(kap)) for j in range(N)])
def Tmv(x):
    Jx = np.einsum('nab,nb->na', J, x.reshape(N, 8)).ravel()
    return 0.5*(x + H2.matvec(x)) + 0.5*(Jx - H1.matvec(Jx))
def dense_T(B):
    """exakter T-Block auf Dreiecksmenge B (8|B| x 8|B|)."""
    out = []
    for k in (k1, k2):
        K = Kernel(m, k).dense(B, B)
        M = mv_times_n(K[..., 0], K[..., 1:], np.broadcast_to(m.n[B][None], K[..., 1:].shape))
        E = np.einsum('ijb,bxy->ixjy', M*(-2/np.sqrt(m.area[B][:, None]*m.area[B][None, :]))[..., None], LB).reshape(8*len(B), 8*len(B))
        out.append(E)
    E1, E2 = out; nb = len(B)
    Jb = np.zeros((8*nb, 8*nb), complex)
    for a in range(nb): Jb[8*a:8*a+8, 8*a:8*a+8] = J[B[a]]
    return 0.5*(np.eye(8*nb) + E2) + 0.5*(Jb - E1@Jb)
def groups(R):
    C = m.cent; V = np.array(list(itertools.product([-1, 1], repeat=3)), float)
    dv = np.linalg.norm(C[:, None, :] - V[None], axis=2); iv = np.argmin(dv, 1)
    E = []
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        for si in [-1, 1]:
            for sj in [-1, 1]: E.append(np.hypot(1 - si*C[:, i], 1 - sj*C[:, j]))
    E = np.stack(E, 1); ie = np.argmin(E, 1)
    lab = np.full(N, -1); lab[E.min(1) < R] = ie[E.min(1) < R]; lab[dv.min(1) < R] = 12 + iv[dv.min(1) < R]
    G = [np.nonzero(lab == g)[0] for g in range(20)]; G = [g for g in G if len(g)]
    return G, np.nonzero(lab == -1)[0]
def make_prec(R):
    G, rest = groups(R); lus = []
    for B in G:
        idx = np.concatenate([np.arange(8*b, 8*b+8) for b in B]); lus.append((idx, sl.lu_factor(dense_T(B))))
    Pd = {r: np.linalg.inv(0.5*(np.eye(8) + J[r])) for r in rest}                     # grobe Naeherung Diagonale
    def apply(v):
        y = v.copy(); V = v.reshape(N, 8); Y = y.reshape(N, 8)
        for r in rest: Y[r] = Pd[r]@V[r]
        for idx, lu in lus: y[idx] = sl.lu_solve(lu, v[idx])
        return y
    return apply, max(len(i) for i, _ in lus), sum(len(i) for i, _ in lus)
def run(P=None):
    cnt = [0]
    op = LinearOperator((8*N, 8*N), matvec=(lambda v: Tmv(v)) if P is None else (lambda v: Tmv(P(v))), dtype=complex)
    b = project(m, lambda x: plane_wave(x, k2, 1.0))
    x, info = gmres(op, b, rtol=1e-6, restart=300, maxiter=3, callback=lambda r: cnt.__setitem__(0, cnt[0]+1), callback_type='pr_norm')
    h = x if P is None else P(x)
    X, W = quad_points(m, 2); hs = (h - b).reshape(N, 8)/np.sqrt(m.area)[:, None]; xh = np.array([0, 0, 1.])
    Gm = gp(np.eye(8)[0] + vec(xh[None, :]), gp(vec(m.n), hs))
    Finf = -(1j*k2/(4*np.pi))*np.einsum('nq,nb->b', np.exp(-1j*k2*X@xh)*W, Gm)
    return cnt[0], 4*np.pi/k2*np.imag(Finf[1])
t0 = time.time()
Pdiag = lambda v: np.einsum('nab,nb->na', np.array([2*np.linalg.inv(np.eye(8) + Jj) for Jj in J]), v.reshape(N, 8)).ravel()
i0, s0 = run(); idg = run(Pdiag)[0] if N <= 1200 else -1
Pb, bmax, btot = make_prec(0.25); ib, sb = run(Pb)
fn = 'cube_h.json'; R = json.load(open(fn)) if os.path.exists(fn) else {}
R[f'{sys.argv[2]}|{L}'] = dict(N=N, its_none=i0, its_diag=idg, its_block=ib, sigma_ext=sb, sigma_ext_noprec=s0,
                               mem_MB=(H1.mem + H2.mem)*16/2**20, block_max=bmax, block_frac=btot/(8*N))
json.dump(R, open(fn, 'w'), indent=1)
print(f"L={L} N={N} kappa={kap}: GMRES ohne {i0}, 2(1+J)^-1 {idg}, Kanten/Ecken R=0.25 {ib}; sigma_ext = {sb:.5f} ({s0:.5f}); "
      f"H {(H1.mem+H2.mem)*16/2**20:.0f} MB; groesster Block {bmax}, Blockanteil {btot/(8*N):.2f}; {time.time()-t0:.0f}s", flush=True)
