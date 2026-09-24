import numpy as np, sys, time, json, os
from mesh import icosphere, Mesh
from hmat import HMatrix, Kernel
from assemble import mv_times_n, LB
n = int(sys.argv[1]); mode = sys.argv[2]; eps = float(sys.argv[3]); k = 1.5
m = Mesh(*icosphere(n)); N = len(m.T)
H = HMatrix(m, k, eps=eps, mode=mode)
rng = np.random.default_rng(0); x = rng.normal(size=8*N) + 1j*rng.normal(size=8*N)
t = time.time(); y = H.matvec(x); tmv = time.time() - t
# exakte Referenz fuer 40 zufaellige Zeilen
rows = rng.choice(N, 40, replace=False); ker = Kernel(m, k)
K = ker.dense(rows, np.arange(N))
M = mv_times_n(K[..., 0], K[..., 1:], np.broadcast_to(m.n[None], K[..., 1:].shape))
Eb = np.einsum('ijb,bxy->ixjy', M*(-2/np.sqrt(m.area[rows][:, None]*m.area[None, :]))[..., None], LB)
y0 = np.einsum('ixjy,jy->ix', Eb, x.reshape(N, 8))
err = np.linalg.norm(y.reshape(N, 8)[rows] - y0)/np.linalg.norm(y0)
fn = 'hmat_scaling.json'; res = json.load(open(fn)) if os.path.exists(fn) else {}
res[f'{n}|{mode}|{eps:g}'] = dict(N=N, mem_MB=H.mem*16/2**20, mem_lr_MB=H.mem_lr*16/2**20, dense_MB=N*N*64*16/2**20, t_build=H.t_build, t_mv=tmv, err=err,
                                  mean_rank=float(np.mean(H.ranks)) if H.ranks else 0, nblocks=len(H.blocks))
json.dump(res, open(fn, 'w'), indent=1)
print(f"n={n} N={N} {mode} eps={eps:g}: Speicher {H.mem*16/2**20:.1f} MB (davon Niedrigrang {H.mem_lr*16/2**20:.1f} MB) (dichte E-Matrix waere {N*N*64*16/2**20:.0f} MB), Aufbau {H.t_build:.0f}s, MatVec {tmv:.2f}s, rel. Fehler (40 Zeilen) {err:.1e}", flush=True)
