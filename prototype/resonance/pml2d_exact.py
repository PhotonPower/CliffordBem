"""Wie pml2d.py, aber die gestreckten Randstuecke werden exakt parametrisiert (keine Sehnen):
y(s) = V_k + xi(s) u,  xi(s) = r~(s) (s <= 1) bzw. 2 - r~(2 - s) (s > 1),  dl = xi'(s) ds,
Gauss in s; Hauptwert am eigenen Element analytisch: int Phi2 dl = u/(2 pi) log((xi_m - xi_0)/(xi_1 - xi_m))."""
import numpy as np
from cl3 import MV
from corner import Jblock
from pml2d import V, vec_mat, E1m, E2m, side_params
GLn = {n: np.polynomial.legendre.leggauss(n) for n in (4, 8, 16)}
def rt(r, th, r0):
    r = np.asarray(r, float); out = r.astype(complex); m = (r < r0) & (r > 0)
    out[m] = r[m]*np.exp(1j*th*np.log(r[m]/r0)); return out
def drt(r, th, r0):
    r = np.asarray(r, float); out = np.ones_like(r, dtype=complex); m = (r < r0) & (r > 0)
    out[m] = (1 + 1j*th)*np.exp(1j*th*np.log(r[m]/r0)); return out
def xi(s, th, r0):
    s = np.asarray(s, float); return np.where(s <= 1, rt(s, th, r0), 2 - rt(2 - s, th, r0))
def dxi(s, th, r0):
    s = np.asarray(s, float); return np.where(s <= 1, drt(s, th, r0), drt(2 - s, th, r0))
class Geo:
    def __init__(self, ng, th, r0, q, q0=0.5):
        s = side_params(ng, q, q0); self.s = s; self.th, self.r0 = th, r0
        S0, S1, SD = [], [], []
        for k in range(4):
            for i in range(len(s) - 1): S0.append(s[i]); S1.append(s[i + 1]); SD.append(k)
        self.S0, self.S1, self.SD = map(np.array, (S0, S1, SD))
        self.U = np.array([(V[(k + 1) % 4] - V[k])/2 for k in range(4)]); self.Nn = np.array([[u[1], -u[0]] for u in self.U])
        self.Sm = 0.5*(self.S0 + self.S1); self.X = self.point(self.Sm, self.SD)
    def point(self, s, side):
        return np.array([V[k] for k in np.atleast_1d(side)]) + xi(s, self.th, self.r0)[..., None]*self.U[np.atleast_1d(side)]
    def quad(self, j, n=16, sub=1):
        """Gausspunkte (Y, w*dl) auf Element j; fuer das Eckelement in s geometrisch zur Ecke verfeinert."""
        s0, s1 = self.S0[j], self.S1[j]; t, w = GLn[n]; t = 0.5*(t + 1); w = 0.5*w
        if s0 == 0.0 or s1 == 2.0:                                   # Eckelement: Teilintervalle geometrisch zur Ecke
            L = s1 - s0; edges = [0.0] + [L*0.5**m for m in range(30, -1, -1)]
            ss = np.concatenate([a + (b - a)*t for a, b in zip(edges[:-1], edges[1:])]); ww = np.concatenate([(b - a)*w for a, b in zip(edges[:-1], edges[1:])])
            ss = s0 + ss if s0 == 0.0 else s1 - ss
        else:
            ss = np.concatenate([(s0 + (s1 - s0)*(m + t)/sub) for m in range(sub)]); ww = np.tile(w*(s1 - s0)/sub, sub)
        side = self.SD[j]
        Y = V[side][None, :] + xi(ss, self.th, self.r0)[:, None]*self.U[side][None, :]
        return Y, ww*dxi(ss, self.th, self.r0)
def build(ng, kappa, th=0.0, r0=0.25, q=0.5, q0=0.5):
    g = Geo(ng, th, r0, q, q0); Ns = len(g.S0); X = g.X
    E = np.zeros((8*Ns, 8*Ns), complex)
    Lr = np.abs(g.S1 - g.S0)
    for j in range(Ns):
        # grobe Abstaende im Parameter (reelle Geometrie) fuer die Nahfeldentscheidung
        Xr = np.array([V[k] for k in g.SD]) + np.where(g.Sm <= 1, g.Sm, g.Sm)[:, None]*g.U[g.SD]
        yj = V[g.SD[j]] + g.Sm[j]*g.U[g.SD[j]]
        near = np.linalg.norm(Xr - yj, axis=1) < 3*Lr[j]
        Vint = np.zeros((Ns, 2), complex)
        for idx, sub in [(~near, 1), (near, 8)]:
            if not idx.any(): continue
            Y, wd = g.quad(j, 16, sub)
            Z = X[idx][:, None, :] - Y[None, :, :]; zz = np.einsum('ngd,ngd->ng', Z, Z)
            Vint[idx] = np.einsum('g,ngd->nd', wd, Z/zz[..., None])/(2*np.pi)
        # eigenes Element: Hauptwert analytisch
        k = g.SD[j]; x0, x1, xm = xi(np.array([g.S0[j], g.S1[j], g.Sm[j]]), th, r0)
        Vint[j] = g.U[k]*np.log((xm - x0)/(x1 - xm))/(2*np.pi)
        M = -2*(Vint[:, 0, None, None]*E1m + Vint[:, 1, None, None]*E2m)@vec_mat(g.Nn[k])
        E[:, 8*j:8*j + 8] = M.reshape(8*Ns, 8)
    Jb = np.zeros_like(E)
    for j in range(Ns): n = g.Nn[g.SD[j]]; Jb[8*j:8*j + 8, 8*j:8*j + 8] = Jblock(np.array([n[0], n[1], 0.]), kappa, 1.0, 1.0, np.sqrt(kappa))
    I = np.eye(8*Ns); return 0.5*(I + E) + 0.5*(I - E)@Jb, g
def far_field(ng, kappa, th=0.0, r0=0.25, q=0.5, x=np.array([5.0, 3.0])):
    T, g = build(ng, kappa, th, r0, q); Ns = len(g.S0)
    hin = np.zeros((Ns, 8), complex); hin[:, 1] = 1.0
    hs = np.linalg.solve(T, hin.ravel()).reshape(Ns, 8) - hin
    F = np.zeros(8, complex)
    for j in range(Ns):
        Y, wd = g.quad(j, 8); Z = x[None, :] - Y; zz = np.einsum('gd,gd->g', Z, Z)
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
        print(f"kappa={k}, theta={th:+.1f}: " + '; '.join(line), flush=True)
