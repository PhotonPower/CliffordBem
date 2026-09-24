"""Schritt 2.3 (erster Teil): diskrete Stabilitaet auf der Kugel. GMRES-Iterationen ohne und mit
punktweiser Vorkonditionierung 2(1+J)^{-1}; kleinster Singulaerwert fuer kleine N."""
import numpy as np, time, sys, json
from scipy.sparse.linalg import gmres, LinearOperator
from mesh import icosphere, Mesh
from assemble import assemble_E, project
from corner import Jblock
from scatter3d import plane_wave
def build(n, om, e1, e2=1.0):
    m = Mesh(*icosphere(n)); N = len(m.T); k1, k2 = om*np.sqrt(e1), om*np.sqrt(e2); re_ = e1/e2
    E1 = assemble_E(m, k1); E2 = assemble_E(m, k2)
    J = np.array([Jblock(m.n[j], re_, 1.0, 1.0, np.sqrt(re_)) for j in range(N)])
    T = 0.5*(np.eye(8*N) + E2) - 0.5*np.einsum('ixjy,jyz->ixjz', E1.reshape(N, 8, N, 8), J).reshape(8*N, 8*N)
    del E1, E2
    for j in range(N): T[8*j:8*j+8, 8*j:8*j+8] += 0.5*J[j]
    P = np.array([2*np.linalg.inv(np.eye(8) + J[j]) for j in range(N)])      # 2(1+J)^{-1}
    return m, T, P, project(m, lambda x: plane_wave(x, k2, e2))
def its(A, b, M=None):
    n = len(b); cnt = [0]
    op = LinearOperator((n, n), matvec=(lambda v: A@v) if M is None else (lambda v: M(A@v)), dtype=complex)
    rhs = b if M is None else M(b)
    x, info = gmres(op, rhs, rtol=1e-8, restart=n, maxiter=1, callback=lambda r: cnt.__setitem__(0, cnt[0]+1), callback_type='pr_norm')
    return cnt[0]
if __name__ == '__main__':
    om, e1 = float(sys.argv[1]), complex(sys.argv[2]); ns = [int(v) for v in sys.argv[3].split(',')]
    for n in ns:
        t0 = time.time(); m, T, P, b = build(n, om, e1); N = len(m.T)
        Mp = lambda v, P=P, N=N: np.einsum('nab,nb->na', P, v.reshape(N, 8)).ravel()
        i0, i1 = its(T, b), its(T, b, Mp)
        smin = np.linalg.svd(T, compute_uv=False)[-1] if N <= 320 else float('nan')
        print(f"n={n} N={N}: GMRES ohne {i0}, mit 2(1+J)^-1 {i1}; sigma_min(T) = {smin:.4f}  ({time.time()-t0:.0f}s)", flush=True)
