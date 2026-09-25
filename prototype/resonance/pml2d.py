"""Resonanzfenster an Ecken, 2D-Quadrat, quasistatisch (k = 0): Kollokation von
T_1 = 1/2 (1 + E) + 1/2 (1 - E) J  mit stueckweise konstanten Multivektor-Dichten auf einem zu den Ecken
geometrisch gradierten Netz, und komplexe Streckung in log r nahe den Ecken ("Mellin-PML"):
    r -> r~ = r (r / r0)^{i theta}   fuer r < r0   (log r~ = log r + i theta (log r - log r0)).
Auf dem gestreckten Rand ist der Kern analytisch fortgesetzt: Phi2(z) = z / (2 pi z.z) mit bilinearem z.z,
dl -> dr~. Die Normalen bleiben reell (die Streckung wirkt entlang der Kante). Beobachtungsgroesse:
Streufeld F^s(x) = int Phi2(x - y) n h^s dl an x = (5, 3) bei homogener Anregung h_inc = e1."""
import numpy as np
from cl3 import MV
from corner import Jblock
V = [np.array(v, float) for v in [(1, -1), (1, 1), (-1, 1), (-1, -1)]]
def vec_mat(v): return MV.vec([v[0], v[1], 0]).left_matrix()
E1m, E2m = vec_mat([1, 0]), vec_mat([0, 1])
GL = {n: np.polynomial.legendre.leggauss(n) for n in (4, 16)}
def side_params(ng, q=0.5, q0=0.5, nmid=None):
    s = [0.0] + [q0*q**k for k in range(ng, -1, -1)]
    if nmid is None: nmid = max(2, int(np.ceil((2 - 2*q0)/(q0*(1 - q)))))   # Mittelteil so fein wie die groessten Eckelemente
    mid = list(np.linspace(q0, 2 - q0, nmid + 1)[1:-1])
    return np.array(s + mid + [2 - x for x in s[::-1]])
def rmap(r, theta, r0):
    r = np.asarray(r, float); out = r.astype(complex)
    m = (r < r0) & (r > 0)
    out[m] = r[m]*np.exp(1j*theta*np.log(r[m]/r0))
    return out
def build(ng, kappa, theta=0.0, r0=0.25, q=0.5, q0=0.5, nmid=None):
    s = side_params(ng, q, q0, nmid)
    P0, P1, Nr, U = [], [], [], []
    for k in range(4):
        A, B = V[k], V[(k + 1) % 4]; u = (B - A)/2; n = np.array([u[1], -u[0]])
        r = np.minimum(s, 2 - s); rt = rmap(r, theta, r0)
        pts = np.array([A + rt[i]*u if s[i] <= 1 else B - rt[i]*u for i in range(len(s))])
        for i in range(len(s) - 1):
            P0.append(pts[i]); P1.append(pts[i + 1]); Nr.append(n); U.append(u)
    P0, P1, Nr, U = map(np.array, (P0, P1, Nr, U))
    Ns = len(P0); X = 0.5*(P0 + P1); D = P1 - P0; dl = np.einsum('jd,jd->j', D, U)   # komplexe Laenge
    Lr = np.abs(dl)
    E = np.zeros((8*Ns, 8*Ns), complex)
    for j in range(Ns):
        # Gauss auf dem Segment, fuer nahe Punkte in 8 Teilsegmente unterteilt
        dist = np.abs(np.sqrt(np.einsum('id,id->i', X - X[j], X - X[j]) + 0j))
        near = dist < 3*Lr[j]
        Vint = np.zeros((Ns, 2), complex)
        for sub, idx in [(1, ~near), (8, near)]:
            if not idx.any(): continue
            t, w = GL[16]; t = 0.5*(t + 1); w = 0.5*w
            T = (np.arange(sub)[:, None] + t[None, :]).ravel()/sub; Wt = np.tile(w, sub)/sub
            Y = P0[j][None, :] + T[:, None]*D[j][None, :]                 # (G,2)
            Z = X[idx][:, None, :] - Y[None, :, :]                          # (n,G,2)
            zz = np.einsum('ngd,ngd->ng', Z, Z)
            Vint[idx] = np.einsum('g,ngd->nd', Wt, Z/zz[..., None])*dl[j]/(2*np.pi)
        Vint[j] = 0.0                                                       # Hauptwert am eigenen Segment
        M = -2*(Vint[:, 0, None, None]*E1m + Vint[:, 1, None, None]*E2m)@vec_mat(Nr[j])
        E[:, 8*j:8*j + 8] = M.reshape(8*Ns, 8)
    Jb = np.zeros_like(E)
    for j in range(Ns): Jb[8*j:8*j + 8, 8*j:8*j + 8] = Jblock(np.array([Nr[j][0], Nr[j][1], 0.]), kappa, 1.0, 1.0, np.sqrt(kappa))
    I = np.eye(8*Ns); T1 = 0.5*(I + E) + 0.5*(I - E)@Jb
    return T1, dict(P0=P0, P1=P1, X=X, dl=dl, Nr=Nr, Lr=Lr)
def far_field(ng, kappa, theta=0.0, r0=0.25, x=np.array([5.0, 3.0]), **kw):
    T, g = build(ng, kappa, theta, r0, **kw); Ns = len(g['X'])
    hin = np.zeros((Ns, 8), complex); hin[:, 1] = 1.0
    h = np.linalg.solve(T, hin.ravel()).reshape(Ns, 8); hs = h - hin
    t, w = GL[4]; t = 0.5*(t + 1); w = 0.5*w; F = np.zeros(8, complex)
    for j in range(Ns):
        Y = g['P0'][j][None, :] + t[:, None]*(g['P1'][j] - g['P0'][j])[None, :]
        Z = x[None, :] - Y; zz = np.einsum('gd,gd->g', Z, Z)
        K = np.einsum('g,gd->d', w, Z/zz[:, None])*g['dl'][j]/(2*np.pi)
        F += (MV.vec([K[0], K[1], 0])*MV.vec([g['Nr'][j][0], g['Nr'][j][1], 0])*MV(hs[j])).c
    s = np.linalg.svd(T, compute_uv=False)
    return F, s[-1], s[0]
if __name__ == '__main__':
    import sys
    for name, k in [('Gold', -11 + 1.2j), ('Fenster -2,8+0,3i', -2.8 + 0.3j)]:
        for th in [0.0, 0.5, -0.5]:
            line = []
            for ng in [8, 16, 24]:
                F, smin, smax = far_field(ng, k, th)
                line.append(f"ng={ng}: F1={F[1].real:+.6f}{F[1].imag:+.6f}i s_min={smin:.1e}")
            print(f"{name}, theta={th:+.1f} | " + '; '.join(line), flush=True)
