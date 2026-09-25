"""Transparente Eckbedingung durch Anreicherung (uebertragbar auf 3D-Kanten, da ohne komplexe Geometrie):
Galerkin wie pml2d_galerkin, Netz nur bis zur Tiefe r_min gradiert; das Eckelement [0, r_min] traegt als
Ansatzfunktion r^{a-1} (physikalischer Singulaerexponent a, alle 8 Komponenten), getestet mit 1
(Petrov-Galerkin, da r^{2(a-1)} fuer Re a < 1/2 nicht integrierbar ist)."""
import numpy as np, mpmath as mp
from cl3 import MV
from corner import Jblock
from pml2d import V, vec_mat, E1m, E2m
from pml2d_galerkin import GalerkinGeo, gauss01, graded_nodes
from pml2d_exact import xi
def exponent(kappa, alpha=np.pi/2):
    """physikalischer Exponent (Re a in (0,1), kleinster Realteil ausser 0) der Kante mit Winkel alpha."""
    c = (1 + kappa)/(1 - kappa); best = None
    for sg in (1, -1):
        f = lambda a: mp.sin((mp.pi - alpha)*a) - sg*c*mp.sin(mp.pi*a)
        for re in np.linspace(0.02, 0.98, 15):
            for im in np.linspace(-2, 2, 9):
                try: r = complex(mp.findroot(f, mp.mpc(re, im)))
                except Exception: continue
                if 1e-6 < r.real < 1 and (best is None or r.real < best.real): best = r
    return best
LEVELS = 300   # Eckelement: r^{a-1} ist fuer kleines Re a nur langsam integrierbar (Schwanz 2^{-levels Re a})
def build(ng, kappa, a, q=0.5, q0=0.5):
    g = GalerkinGeo(ng, 0.0, 0.25, q, q0); Ns = len(g.S0)
    corner = (g.S0 == 0.0) | (g.S1 == 2.0)
    Tn, Wn = gauss01(8); Tg, Wg = graded_nodes(8, 14, True)
    def pts(j, graded=False):
        if corner[j]:
            Y, w, xx = g.pts_corner(j, levels=LEVELS); off = g._off; r = np.linalg.norm(off.real, axis=1)
            return Y, w, xx, off, r
        Y, w, xx = g.pts(j, Tg, Wg) if graded else g.pts(j, Tn, Wn); return Y, w, xx, None, None
    P = [pts(j) for j in range(Ns)]; Pg = [pts(j, True) for j in range(Ns)]
    trial = lambda j, Q: Q[1]*(Q[4]**(a - 1) if corner[j] else 1.0)          # Gewicht * Ansatzfunktion
    def key(k, s):
        if s == 0.0: return ('V', k)
        if s == 2.0: return ('V', (k + 1) % 4)
        return ('E', k, round(s, 14))
    K0 = [key(g.SD[j], g.S0[j]) for j in range(Ns)]; K1 = [key(g.SD[j], g.S1[j]) for j in range(Ns)]
    Mass = np.array([np.sum(trial(j, P[j])) for j in range(Ns)])            # int 1 * phi_j
    Vall = np.zeros((Ns, Ns, 2), complex)
    for j in range(Ns):
        for i in range(Ns):
            adj = len({K0[i], K1[i]} & {K0[j], K1[j]}) > 0
            if i == j and not corner[j]: continue                            # Selbstterm 0 (antisymmetrisch)
            if not adj:
                X, wx = P[i][0], P[i][1]; Y, wy = P[j][0], trial(j, P[j])
                Z = X[:, None, :] - Y[None, :, :]
            elif i == j or g.SD[i] != g.SD[j] or corner[j]:
                # Eckelement mit sich selbst oder Paare an derselben Ecke: gradiert x gradiert (Versaetze)
                Qi, Qj = Pg[i], Pg[j]; wx, wy = Qi[1], trial(j, Qj)
                if corner[i] and corner[j] and (i == j or g.SD[i] != g.SD[j]):
                    Z = Qi[3][:, None, :] - Qj[3][None, :, :]
                else:
                    Z = Qi[0][:, None, :] - Qj[0][None, :, :]
            else:                                                            # kollinear, beide normal: analytisch
                X, wx, xx = Pg[i][0], Pg[i][1], Pg[i][2]
                x0, x1 = xi(np.array([g.S0[j], g.S1[j]]), 0.0, 0.25)
                Vall[i, j] = g.U[g.SD[j]]*np.sum(wx*np.log((xx - x0)/(xx - x1)))/(2*np.pi); continue
            zz = np.einsum('abd,abd->ab', Z, Z)
            with np.errstate(divide='ignore', invalid='ignore'):
                Kz = np.where(zz[..., None] != 0, Z/zz[..., None], 0.0)
            Vall[i, j] = np.einsum('a,b,abd->d', wx, wy, Kz)/(2*np.pi)
    Eg = np.zeros((8*Ns, 8*Ns), complex)
    for j in range(Ns):
        n = g.Nn[g.SD[j]]; M = -2*(Vall[:, j, 0, None, None]*E1m + Vall[:, j, 1, None, None]*E2m)@vec_mat(n)
        Eg[:, 8*j:8*j + 8] = M.reshape(8*Ns, 8)
    Jb = np.zeros_like(Eg)
    for j in range(Ns): n = g.Nn[g.SD[j]]; Jb[8*j:8*j + 8, 8*j:8*j + 8] = Jblock(np.array([n[0], n[1], 0.]), kappa, 1.0, 1.0, np.sqrt(kappa))
    Mm = np.diag(np.repeat(Mass, 8)); G = 0.5*(Mm + Eg) + 0.5*(Mm - Eg)@Jb
    Lr = np.array([np.sum(P[j][1]) for j in range(Ns)])                     # int 1 dl (rechte Seite)
    return G, Lr, g, P, trial
def far_field(ng, kappa, a, q=0.5, x=np.array([5.0, 3.0])):
    G, Lr, g, P, trial = build(ng, kappa, a, q); Ns = len(Lr)
    hin = np.zeros((Ns, 8), complex); hin[:, 1] = 1.0
    c = np.linalg.solve(G, (np.repeat(Lr, 8)*hin.ravel())).reshape(Ns, 8)
    F = np.zeros(8, complex)
    for j in range(Ns):
        Y, w = P[j][0], P[j][1]; Z = x[None, :] - Y; zz = np.einsum('gd,gd->g', Z, Z); n = g.Nn[g.SD[j]]
        Kt = np.einsum('g,gd->d', trial(j, P[j]), Z/zz[:, None])/(2*np.pi)    # mit Ansatzfunktion
        K1 = np.einsum('g,gd->d', w, Z/zz[:, None])/(2*np.pi)                  # fuer h_inc (konstant)
        nv = MV.vec([n[0], n[1], 0])
        F += (MV.vec([Kt[0], Kt[1], 0])*nv*MV(c[j])).c - (MV.vec([K1[0], K1[1], 0])*nv*MV(hin[j])).c
    return F
if __name__ == '__main__':
    import sys
    k = complex(sys.argv[1]); q = float(sys.argv[2]); depths = [float(v) for v in sys.argv[3].split(',')]
    a = exponent(k); print(f"kappa = {k}: Exponent a = {a:.4f}")
    for dep in depths:
        ng = int(np.ceil(np.log(dep/0.5)/np.log(q))); F = far_field(ng, k, a, q)
        print(f"  q={q}, Tiefe {dep:.0e} (ng={ng}): F1 = {F[1].real:+.7f}{F[1].imag:+.7f}i", flush=True)
