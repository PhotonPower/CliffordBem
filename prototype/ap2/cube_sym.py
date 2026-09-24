"""Symmetrieangepasste Assemblierung auf dem Wuerfel: drei Spiegelungen x_i -> -x_i,
Wirkung auf Multivektorfelder (U h)(x) = e_i h^(S x) e_i (Hauptinvolution). 8 Sektoren (Charaktere),
T_chi[s,t] = sum_g chi(g) T[s, g t] A_g  fuer s, t im Fundamentalbereich (Oktant x,y,z > 0)."""
import numpy as np, itertools, time
from cl3 import MV
from mesh import quad_points
from analytic import tri_integrals
from assemble import kernel_full, kernel_rem, mv_times_n, LB
from corner import Jblock
np.seterr(divide='ignore', invalid='ignore')
GROUP = list(itertools.product([1, -1], repeat=3))
def sandwich(i):
    """8x8-Matrix von c -> e_i c^ e_i."""
    e = MV.blade([1, 2, 4][i]); A = np.zeros((8, 8))
    for b in range(8):
        c = MV.blade(b); ch = c.involute(); A[:, b] = (e*ch*e).c.real
    return A
SW = [sandwich(i) for i in range(3)]
def A_of(g):
    A = np.eye(8)
    for i in range(3):
        if g[i] < 0: A = SW[i]@A
    return A
def perm_of(mesh, g):
    key = {tuple(np.round(c, 10)): i for i, c in enumerate(mesh.cent)}
    return np.array([key[tuple(np.round(c*np.array(g), 10))] for c in mesh.cent])
def fundamental(mesh):
    return np.nonzero(np.all(mesh.cent > 0, axis=1))[0]
def assemble_rows(mesh, k, rows, near_fac=2.5, sub_near=4, chunk=16):
    """E_k-Zeilenbloecke fuer die Dreiecke rows: (R, N, 8, 8)."""
    N = len(mesh.T); R = len(rows)
    X, W = quad_points(mesh); Xo, Wo = quad_points(mesh, sub_near)
    S = np.zeros((R, N), complex); Vv = np.zeros((R, N, 3), complex)
    for a0 in range(0, R, chunk):
        ii = rows[a0:a0+chunk]
        Z = X[ii][:, :, None, None, :] - X[None, None, :, :, :]
        vv, ss = kernel_full(Z, k); ww = W[ii][:, :, None, None]*W[None, None, :, :]
        S[a0:a0+len(ii)] = np.einsum('iqjp,iqjp->ij', ww, ss); Vv[a0:a0+len(ii)] = np.einsum('iqjp,iqjpd->ijd', ww, vv)
    D = np.linalg.norm(mesh.cent[rows][:, None, :] - mesh.cent[None, :, :], axis=-1)
    near = D < near_fac*np.maximum(mesh.h[rows][:, None], mesh.h[None, :])
    for a, j in zip(*np.nonzero(near)):
        i = rows[a]
        if i == j or mesh.h[i] <= 1.5*mesh.h[j]:
            xo, wo = Xo[i], Wo[i]
            Ig, Ii = tri_integrals(xo, mesh.V[j], mesh.n[j])
            v = (wo@Ig)/(4*np.pi) if i != j else np.zeros(3); s = -1j*k*(wo@Ii)/(4*np.pi)
            Z = xo[:, None, :] - X[j][None, :, :]; ww = wo[:, None]*W[j][None, :]
        else:   # aeusseres Integral ueber das kleinere Dreieck j, inneres analytisch ueber i
            yo, wo = Xo[j], Wo[j]
            Ig, Ii = tri_integrals(yo, mesh.V[i], mesh.n[i])
            v = -(wo@Ig)/(4*np.pi); s = -1j*k*(wo@Ii)/(4*np.pi)
            Z = X[i][:, None, :] - yo[None, :, :]; ww = W[i][:, None]*wo[None, :]
        vr, sr = kernel_rem(Z, k)
        S[a, j] = s + np.sum(ww*sr); Vv[a, j] = v + np.einsum('qp,qpd->d', ww, vr)
    M = mv_times_n(S, Vv, np.broadcast_to(mesh.n[None, :, :], Vv.shape))
    scale = -2/np.sqrt(mesh.area[rows][:, None]*mesh.area[None, :])
    return np.einsum('ijb,bxy->ijxy', M*scale[..., None], LB), near.sum()
def sector_matrices(mesh, om, e1, e2=1.0):
    """Liefert dict chi -> T_chi (dim 8|F|) fuer T_1 mit Standardwahl."""
    F = fundamental(mesh); nF = len(F); N = len(mesh.T)
    k1, k2 = om*np.sqrt(e1), om*np.sqrt(e2); re_ = e1/e2
    J = np.array([Jblock(mesh.n[j], re_, 1.0, 1.0, np.sqrt(re_)) for j in range(N)])
    t0 = time.time()
    perms = {g: perm_of(mesh, g)[F] for g in GROUP}
    R2, nn = assemble_rows(mesh, k2, F)
    Bg = {g: np.einsum('stxy,yz->stxz', R2[:, perms[g]], A_of(g)) for g in GROUP}
    del R2
    R1, _ = assemble_rows(mesh, k1, F)
    for g in GROUP:
        p = perms[g]
        Bg[g] -= np.einsum('stxy,tyz,zw->stxw', R1[:, p], J[p], A_of(g))
    del R1
    ta = time.time() - t0
    base = np.zeros((nF, nF, 8, 8), complex)
    for a in range(nF): base[a, a] = 0.5*(np.eye(8) + J[F[a]])
    out = {}
    for chi in GROUP:
        acc = base.copy()
        for g in GROUP:
            acc += 0.5*np.prod([c if gi < 0 else 1 for c, gi in zip(chi, g)])*Bg[g]
        out[chi] = acc.transpose(0, 2, 1, 3).reshape(8*nF, 8*nF)
    return out, ta, nn
if __name__ == '__main__':
    from cube_mesh import cube_mesh
    from assemble import assemble_E
    m = cube_mesh(1); N = len(m.T)
    # Kommutierungstest mit der vollen Matrix
    E = assemble_E(m, 1.3 + 0.2j)
    for g in GROUP[1:]:
        p = perm_of(m, g); A = A_of(g)
        U = np.zeros((8*N, 8*N))
        for t in range(N): U[8*p[t]:8*p[t]+8, 8*t:8*t+8] = A
        print(g, 'Kommutator |EU - UE|/|E| =', np.linalg.norm(E@U - U@E)/np.linalg.norm(E))
    # J kommutiert?
    J = [Jblock(m.n[j], -3+0.3j, 1.0, 1.0, np.sqrt(-3+0.3j)) for j in range(N)]
    g = (-1, 1, -1); p = perm_of(m, g); A = A_of(g)
    print('J-Kommutator:', max(np.abs(J[p[t]]@A - A@J[t]).max() for t in range(N)))
    # Sektorzerlegung gegen volle Matrix: sigma_min
    from cube_study import full_T
    Tf = full_T(m, 0.5, -3+0.3j)
    smin_full = np.linalg.svd(Tf, compute_uv=False)[-1]
    sec, _, _ = sector_matrices(m, 0.5, -3+0.3j)
    smin_sec = min(np.linalg.svd(Ts, compute_uv=False)[-1] for Ts in sec.values())
    print('sigma_min voll:', smin_full, ' min ueber Sektoren:', smin_sec)
