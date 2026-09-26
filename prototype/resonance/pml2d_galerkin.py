"""Galerkin-Variante (stueckweise konstante Ansatz- und Testfunktionen) des 2D-Prototyps mit komplexer
Streckung in log r. Testung mit dem Konturintegral int (.) dl~ (bilinear, ohne Konjugation), damit das
diskrete Problem die analytische Fortsetzung respektiert:
  G = 1/2 (M + Eg) + 1/2 (M - Eg) J,  M_ii = int_tau_i dl~,  Eg_ij = int_tau_i int_tau_j E-Kern dl~ dl~,
  G h = M h_inc.
Selbstterm des Doppelintegrals = 0 (antisymmetrischer Kern). Benachbarte Elemente auf derselben Kante:
inneres Integral analytisch (log), aeusseres gradiert zum gemeinsamen Endpunkt; uebrige Paare Gauss x Gauss,
nahe Paare unterteilt."""
import numpy as np
from cl3 import MV
from corner import Jblock
from pml2d import V, vec_mat, E1m, E2m
from pml2d_exact import Geo, xi, dxi
GL = {n: np.polynomial.legendre.leggauss(n) for n in (8, 16)}
def gauss01(n): t, w = GL[n]; return 0.5*(t + 1), 0.5*w
def graded_nodes(n=8, levels=12, both=True):
    """Knoten/Gewichte auf [0,1], geometrisch zu 0 (und 1) verfeinert."""
    t, w = gauss01(n)
    if both:
        edges = [0.0] + [0.5*0.5**m for m in range(levels, 0, -1)] + [0.5]
        e2 = [1 - e for e in edges[::-1]]; edges = edges + e2[1:]
    else:
        edges = [0.0] + [0.5**m for m in range(levels, 0, -1)] + [1.0]
    T = np.concatenate([a + (b - a)*t for a, b in zip(edges[:-1], edges[1:])])
    W = np.concatenate([(b - a)*w for a, b in zip(edges[:-1], edges[1:])])
    return T, W
class GalerkinGeo(Geo):
    def pts(self, j, T, W):
        s0, s1 = self.S0[j], self.S1[j]; ss = s0 + (s1 - s0)*T; k = self.SD[j]
        Y = V[k][None, :] + xi(ss, self.th, self.r0)[:, None]*self.U[k][None, :]
        return Y, W*(s1 - s0)*dxi(ss, self.th, self.r0), xi(ss, self.th, self.r0)
    def pts_corner(self, j, n=8, levels=40):
        """Eckelement: Teilintervalle geometrisch zur Ecke, direkt im Eckabstand r (ohne Ausloeschung in 2 - s)."""
        from pml2d_exact import rt, drt
        s0, s1 = self.S0[j], self.S1[j]; L = s1 - s0; t, w = gauss01(n)
        edges = [0.0] + [L*0.5**m for m in range(levels, -1, -1)]
        rr = np.concatenate([a + (b - a)*t for a, b in zip(edges[:-1], edges[1:])]); ww = np.concatenate([(b - a)*w for a, b in zip(edges[:-1], edges[1:])])
        k = self.SD[j]; r_t = rt(rr, self.th, self.r0)
        if s0 == 0.0: off = r_t[:, None]*self.U[k][None, :]; Y = V[k][None, :] + off; xx = r_t
        else: off = -r_t[:, None]*self.U[k][None, :]; Y = V[(k + 1) % 4][None, :] + off; xx = 2 - r_t
        self._off = off                                  # Versatz zur Ecke (fuer Paare an derselben Ecke)
        return Y, ww*drt(rr, self.th, self.r0), xx
def build(ng, kappa, th=0.0, r0=0.25, q=0.5, q0=0.5, parts_only=False):
    if r0 < q0: r0 = q0*q**round(np.log(r0/q0)/np.log(q))     # Beginn der Streckung auf einen Netzknoten legen (Knick der Abbildung)
    g = GalerkinGeo(ng, th, r0, q, q0); Ns = len(g.S0)
    corner = (g.S0 == 0.0) | (g.S1 == 2.0)
    Tn, Wn = gauss01(8); Tg, Wg = graded_nodes(8, 14, True)
    P = []                                                                           # Standardpunkte je Element
    for j in range(Ns): P.append(g.pts_corner(j) if corner[j] else g.pts(j, Tn, Wn))
    Pg = [g.pts_corner(j) if corner[j] else g.pts(j, Tg, Wg) for j in range(Ns)]    # gradiert (fuer Nachbarn)
    # Nachbarschaft: gemeinsamer Endpunkt (gleiche Kante oder Ecke)
    endpoints = [(g.SD[j], g.S0[j]) for j in range(Ns)]
    def key(k, s):   # Punkt identifizieren: Ecke oder Kantenpunkt
        if s == 0.0: return ('V', k)
        if s == 2.0: return ('V', (k + 1) % 4)
        return ('E', k, round(s, 14))
    K0 = [key(g.SD[j], g.S0[j]) for j in range(Ns)]; K1 = [key(g.SD[j], g.S1[j]) for j in range(Ns)]
    Mass = np.array([np.sum(P[j][1]) for j in range(Ns)])
    Eg = np.zeros((8*Ns, 8*Ns), complex)
    Vall = np.zeros((Ns, Ns, 2), complex)
    for j in range(Ns):
        Y, wy, _ = P[j]
        for i in range(Ns):
            if i == j: continue
            adj = len({K0[i], K1[i]} & {K0[j], K1[j]}) > 0
            if not adj:
                X, wx, _ = P[i]
                Z = X[:, None, :] - Y[None, :, :]; zz = np.einsum('abd,abd->ab', Z, Z)
                Vall[i, j] = np.einsum('a,b,abd->d', wx, wy, Z/zz[..., None])/(2*np.pi)
            else:
                X, wx, xx = Pg[i]
                if g.SD[i] == g.SD[j]:                                           # kollinear: inneres Integral analytisch
                    x0, x1 = xi(np.array([g.S0[j], g.S1[j]]), th, r0)
                    inner = np.log((xx - x0)/(xx - x1))
                    Vall[i, j] = g.U[g.SD[j]]*np.sum(wx*inner)/(2*np.pi)
                else:                                                            # Ecke: Gauss gradiert x gradiert, Versaetze zur Ecke
                    g.pts_corner(i); Oi = g._off; g.pts_corner(j); Oj = g._off; wyc = Pg[j][1]
                    Z = Oi[:, None, :] - Oj[None, :, :]; zz = np.einsum('abd,abd->ab', Z, Z)
                    Vall[i, j] = np.einsum('a,b,abd->d', wx, wyc, Z/zz[..., None])/(2*np.pi)
    for j in range(Ns):
        n = g.Nn[g.SD[j]]
        M = -2*(Vall[:, j, 0, None, None]*E1m + Vall[:, j, 1, None, None]*E2m)@vec_mat(n)
        Eg[:, 8*j:8*j + 8] = M.reshape(8*Ns, 8)
    if parts_only: return Eg, Mass, g, P
    return assemble(Eg, Mass, g, [kappa]*Ns), Mass, g, P
def assemble(Eg, Mass, g, kappas):
    """G = 1/2 (M + Eg) + 1/2 (M - Eg) J mit einem Kontrast je Element (quasistatisch geht kappa nur ueber J ein)."""
    Ns = len(Mass); Mm = np.diag(np.repeat(Mass, 8)); Jb = np.zeros_like(Eg)
    for j in range(Ns):
        n = g.Nn[g.SD[j]]; kj = kappas[j]
        Jb[8*j:8*j + 8, 8*j:8*j + 8] = Jblock(np.array([n[0], n[1], 0.]), kj, 1.0, 1.0, np.sqrt(kj))
    return 0.5*(Mm + Eg) + 0.5*(Mm - Eg)@Jb
def far_field(ng, kappa, th=0.0, r0=0.25, q=0.5, x=np.array([5.0, 3.0])):
    G, Mass, g, P = build(ng, kappa, th, r0, q); Ns = len(Mass)
    hin = np.zeros((Ns, 8), complex); hin[:, 1] = 1.0
    h = np.linalg.solve(G, (np.repeat(Mass, 8)*hin.ravel())).reshape(Ns, 8); hs = h - hin
    F = np.zeros(8, complex)
    for j in range(Ns):
        Y, wd, _ = P[j]; Z = x[None, :] - Y; zz = np.einsum('gd,gd->g', Z, Z)
        K = np.einsum('g,gd->d', wd, Z/zz[:, None])/(2*np.pi); n = g.Nn[g.SD[j]]
        F += (MV.vec([K[0], K[1], 0])*MV.vec([n[0], n[1], 0])*MV(hs[j])).c
    return F
if __name__ == '__main__':
    import sys
    k = complex(sys.argv[1]); ths = [float(v) for v in sys.argv[2].split(',')]; qs = [float(v) for v in sys.argv[3].split(',')]; depth = float(sys.argv[4])
    r0 = float(sys.argv[5]) if len(sys.argv) > 5 else 0.25
    for th in ths:
        line = []
        for q in qs:
            ng = int(np.ceil(np.log(depth/0.5)/np.log(q))); F = far_field(ng, k, th, r0=r0, q=q)
            line.append(f"q={q} (ng={ng}): {F[1].real:+.7f}{F[1].imag:+.7f}i")
        print(f"Galerkin kappa={k}, theta={th:+.1f}: " + '; '.join(line), flush=True)
