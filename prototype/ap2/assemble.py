"""Galerkin-Assemblierung des Cauchy-Randoperators E_k (orthonormale stueckweise konstante Basis
chi_i/sqrt|tau_i|, 8 Komponenten je Dreieck) auf ebenen Dreiecken.
Kern Phi_k = Phi_0 + (-ik/(4 pi r)) + Rem;  Nahpaare: Phi_0 und 1/r analytisch (innen), Rem mit Gauss."""
import numpy as np, time
np.seterr(divide="ignore", invalid="ignore")
from cl3 import MV
from mesh import quad_points
from analytic import tri_integrals
LB = np.array([MV.blade(b).left_matrix() for b in range(8)])        # (8,8,8): Linksmultiplikation mit Blade b
def kernel_full(Z, k):
    r = np.linalg.norm(Z, axis=-1); Phi = np.exp(1j*k*r)/(4*np.pi*r)
    return (Phi*(1/r - 1j*k)/r)[..., None]*Z, -1j*k*Phi              # Vektorteil, Skalarteil
def kernel_rem(Z, k):
    r = np.linalg.norm(Z, axis=-1); rs = np.maximum(r, 1e-14); e = np.exp(1j*k*rs)
    av = (e*(1 - 1j*k*rs) - 1)/(4*np.pi*rs**3)                      # Koeffizient vor Z (= a/r)
    bs = -1j*k*(e - 1)/(4*np.pi*rs)
    small = r < 1e-10
    av = np.where(small, 0, av); bs = np.where(small, k**2/(4*np.pi), bs)
    return av[..., None]*Z, bs
def mv_times_n(s, v, n):
    """Multivektor (s + v) n -> 8 Komponenten; s: (...), v: (...,3), n: (...,3)."""
    M = np.zeros(s.shape + (8,), complex)
    M[..., 0] = np.sum(v*n, -1); M[..., 1] = s*n[..., 0]; M[..., 2] = s*n[..., 1]; M[..., 4] = s*n[..., 2]
    w = np.cross(v, n); M[..., 6] = w[..., 0]; M[..., 5] = -w[..., 1]; M[..., 3] = w[..., 2]
    return M
def assemble_E(mesh, k, near_fac=2.5, sub_near=4, chunk=24, verbose=False):
    t0 = time.time(); N = len(mesh.T)
    X, W = quad_points(mesh)                                          # (N,7,3), (N,7)
    S = np.zeros((N, N), complex); Vv = np.zeros((N, N, 3), complex)
    for i0 in range(0, N, chunk):
        i1 = min(N, i0 + chunk)
        Z = X[i0:i1, :, None, None, :] - X[None, None, :, :, :]
        vv, ss = kernel_full(Z, k)
        ww = W[i0:i1, :, None, None]*W[None, None, :, :]
        S[i0:i1] = np.einsum('iqjp,iqjp->ij', ww, ss); Vv[i0:i1] = np.einsum('iqjp,iqjpd->ijd', ww, vv)
    # Nahpaare
    D = np.linalg.norm(mesh.cent[:, None, :] - mesh.cent[None, :, :], axis=-1)
    near = D < near_fac*np.maximum(mesh.h[:, None], mesh.h[None, :])
    Xo, Wo = quad_points(mesh, sub_near)
    for i, j in zip(*np.nonzero(near)):
        xo, wo = Xo[i], Wo[i]
        Ig, Ii = tri_integrals(xo, mesh.V[j], mesh.n[j])
        v = (wo@Ig)/(4*np.pi) if i != j else np.zeros(3)             # Phi_0: Selbstterm exakt 0 (Antisymmetrie)
        s = -1j*k*(wo@Ii)/(4*np.pi)
        Z = xo[:, None, :] - X[j][None, :, :]
        vr, sr = kernel_rem(Z, k); ww = wo[:, None]*W[j][None, :]
        v = v + np.einsum('qp,qpd->d', ww, vr); s = s + np.sum(ww*sr)
        S[i, j] = s; Vv[i, j] = v
    M = mv_times_n(S, Vv, np.broadcast_to(mesh.n[None, :, :], Vv.shape))     # (N,N,8)
    scale = -2/np.sqrt(mesh.area[:, None]*mesh.area[None, :])
    E = np.einsum('ijb,bxy->ixjy', M*scale[..., None], LB).reshape(8*N, 8*N)
    if verbose: print(f"  Assemblierung N={N}: {time.time()-t0:.1f}s, Nahpaare {near.sum()}", flush=True)
    return E
def project(mesh, f):
    """L^2-Projektion einer Funktion f: (M,3) -> (M,8) auf die orthonormale Basis: (8N,)."""
    X, W = quad_points(mesh, 2); N = len(mesh.T)
    vals = f(X.reshape(-1, 3)).reshape(N, -1, 8)
    return (np.einsum('nq,nqb->nb', W, vals)/np.sqrt(mesh.area)[:, None]).ravel()
