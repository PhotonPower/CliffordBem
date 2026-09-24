"""Streuung an der Kugel mit H-Matrix-Operatoren (ACA) und GMRES, fuer Netze jenseits der dichten Grenze.
Aufruf: python3 scatter_h.py n omega eps1 stage   (stage: b1 | b2 | solve)"""
import numpy as np, sys, pickle, json, os, time
from scipy.sparse.linalg import gmres, LinearOperator
from mesh import icosphere, Mesh, quad_points
from hmat import HMatrix
from assemble import project
from corner import Jblock
from scatter3d import plane_wave
from fields import vec, gp
import mie
n, om, e1, stage = int(sys.argv[1]), float(sys.argv[2]), complex(sys.argv[3]), sys.argv[4]
eps = 1e-4; e2 = 1.0; k1, k2 = om*np.sqrt(e1), om*np.sqrt(e2)
m = Mesh(*icosphere(n)); N = len(m.T); tag = f'h_{n}_{om}_{sys.argv[3]}'
if stage in ('b1', 'b2'):
    k = k1 if stage == 'b1' else k2
    H = HMatrix(m, k, eps=eps, mode='joint'); del H.ker
    pickle.dump(H, open(f'{tag}_{stage}.pkl', 'wb')); sys.exit()
H1 = pickle.load(open(f'{tag}_b1.pkl', 'rb')); H2 = pickle.load(open(f'{tag}_b2.pkl', 'rb'))
re_ = e1/e2; J = np.array([Jblock(m.n[j], re_, 1.0, 1.0, np.sqrt(re_)) for j in range(N)])
P = np.array([2*np.linalg.inv(np.eye(8) + Jj) for Jj in J])
def Tmv(x):
    Jx = np.einsum('nab,nb->na', J, x.reshape(N, 8)).ravel()
    return 0.5*(x + H2.matvec(x)) + 0.5*(Jx - H1.matvec(Jx))
prec = lambda v: np.einsum('nab,nb->na', P, v.reshape(N, 8)).ravel()
b = project(m, lambda x: plane_wave(x, k2, e2)); cnt = [0]; t0 = time.time()
op = LinearOperator((8*N, 8*N), matvec=lambda v: prec(Tmv(v)), dtype=complex)
h, info = gmres(op, prec(b), rtol=1e-6, restart=200, maxiter=5, callback=lambda r: cnt.__setitem__(0, cnt[0]+1), callback_type='pr_norm')
res = np.linalg.norm(Tmv(h) - b)/np.linalg.norm(b)
X, W = quad_points(m, 2); hs = (h - b).reshape(N, 8)/np.sqrt(m.area)[:, None]; xh = np.array([0, 0, 1.])
G = gp(np.eye(8)[0] + vec(xh[None, :]), gp(vec(m.n), hs))
Finf = -(1j*k2/(4*np.pi))*np.einsum('nq,nb->b', np.exp(-1j*k2*X@xh)*W, G)
Q = (4*np.pi/k2*np.imag(Finf[1]/np.sqrt(e2)))/np.pi
fn = 'scatter_h.json'; R = json.load(open(fn)) if os.path.exists(fn) else {}
R[f'{om}|{sys.argv[3]}|{n}'] = dict(N=N, Q=Q, its=cnt[0], res=res, mem_MB=(H1.mem + H2.mem)*16/2**20)
json.dump(R, open(fn, 'w'), indent=1)
print(f"n={n} N={N} ({8*N} Unbekannte): Qext = {Q:.5f} (Mie {mie.qext(om, e1, e2):.5f}), GMRES {cnt[0]} It., Residuum {res:.1e}, H-Speicher {(H1.mem+H2.mem)*16/2**20:.0f} MB, Loesung {time.time()-t0:.0f}s", flush=True)
