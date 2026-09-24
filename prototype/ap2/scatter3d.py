"""3D-Galerkin-Loesung von T_1 h = h_inc auf Ikosaeder-Kugeln, Extinktion gegen Mie."""
import numpy as np, time, sys, json
from mesh import icosphere, Mesh, quad_points
from assemble import assemble_E, project, mv_times_n, LB
from corner import Jblock
from fields import vec, gp, PSEUDO
import mie
def plane_wave(X, k, eps, d=np.array([0, 0, 1.]), p=np.array([1, 0, 0.])):
    ph = np.exp(1j*k*X@d)[:, None]
    return vec(np.sqrt(eps)*p[None, :]*ph) + gp(PSEUDO, vec(np.sqrt(eps)*np.cross(d, p)[None, :]*ph))
def solve(n, om, e1, e2=1.0, verbose=True):
    m = Mesh(*icosphere(n)); N = len(m.T)
    k1, k2 = om*np.sqrt(e1), om*np.sqrt(e2); re_ = e1/e2
    t0 = time.time()
    E1 = assemble_E(m, k1); E2 = assemble_E(m, k2)
    J = np.array([Jblock(m.n[j], re_, 1.0, 1.0, np.sqrt(re_)) for j in range(N)])       # Standardwahl
    I = np.eye(8*N)
    EJ = np.einsum('ixjy,jyz->ixjz', E1.reshape(N, 8, N, 8), J).reshape(8*N, 8*N)
    T = 0.5*(I + E2) - 0.5*EJ; del EJ, E1, E2
    for j in range(N): T[8*j:8*j+8, 8*j:8*j+8] += 0.5*J[j]
    hin = project(m, lambda x: plane_wave(x, k2, e2))
    h = np.linalg.solve(T, hin)
    cond_est = None
    # Fernfeld: F_inf(xh) = -(i k2/4pi) int e^{-i k2 xh.y} (1 + xh) n h^s dS
    X, W = quad_points(m, 2); hs = (h - hin).reshape(N, 8)/np.sqrt(m.area)[:, None]
    xh = np.array([0, 0, 1.])
    ph = np.exp(-1j*k2*X@xh)*W                                  # (N,q)
    # (1 + xh) n h^s : erst n h^s, dann (1 + xh) von links
    nm = vec(m.n); nh = gp(nm, hs)                               # (N,8)
    G = gp(np.eye(8)[0] + vec(xh[None, :]), nh)                  # (N,8)
    Finf = -(1j*k2/(4*np.pi))*np.einsum('nq,nb->b', ph, G)
    Einf = Finf[[1, 2, 4]]/np.sqrt(e2)
    Qext = (4*np.pi/k2*np.imag(Einf[0]))/np.pi
    aux = np.linalg.norm(h.reshape(N, 8)[:, [0, 7]])/np.linalg.norm(h)
    if verbose: print(f"  n={n} N={N} ({time.time()-t0:.0f}s): Qext = {Qext:.5f}  (Mie {mie.qext(om, e1, e2):.5f}),  Hilfsanteil {aux:.1e}", flush=True)
    return Qext, aux, N
if __name__ == '__main__':
    om, e1 = float(sys.argv[1]), complex(sys.argv[2]); ns = [int(v) for v in sys.argv[3].split(',')]
    res = [solve(n, om, e1) for n in ns]
    json.dump({'om': om, 'e1': [e1.real, e1.imag], 'res': res, 'mie': mie.qext(om, e1)}, open(f'scatter_{om}_{sys.argv[2]}.json', 'w'))
