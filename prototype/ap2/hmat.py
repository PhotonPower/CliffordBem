"""AP 3: H-Matrix mit Adaptiver Kreuzapproximation (ACA) fuer den Cauchy-Randoperator E_k.
Komprimiert werden die Kernkomponenten K_c(i,j) = int int [s, v_x, v_y, v_z](x - y), c = 0..3
(Skalar- und Vektoranteil von Phi_k); daraus E_ij = -2/sqrt(A_i A_j) L((s + v) n_j).
Varianten: 'comp' = vier getrennte ACAs, 'joint' = eine ACA auf der gestapelten Matrix R x (4C)."""
import numpy as np, time
from mesh import quad_points
from analytic import tri_integrals
from assemble import kernel_full, kernel_rem, mv_times_n, LB
np.seterr(divide='ignore', invalid='ignore')
class Kernel:
    def __init__(self, mesh, k, near_fac=2.5, sub_near=4):
        self.m, self.k = mesh, k
        self.X, self.W = quad_points(mesh); self.Xo, self.Wo = quad_points(mesh, sub_near); self.near_fac = near_fac
    def far(self, rows, cols):
        """Kernkomponenten (R, C, 4) mit Gauss 7x7 (nur fuer getrennte Paare korrekt)."""
        Z = self.X[rows][:, :, None, None, :] - self.X[cols][None, None, :, :, :]
        vv, ss = kernel_full(Z, self.k); ww = self.W[rows][:, :, None, None]*self.W[cols][None, None, :, :]
        K = np.empty((len(rows), len(cols), 4), complex)
        K[..., 0] = np.einsum('iqjp,iqjp->ij', ww, ss); K[..., 1:] = np.einsum('iqjp,iqjpd->ijd', ww, vv)
        return K
    def dense(self, rows, cols):
        """wie far, aber mit Nahfeldbehandlung (fuer unzulaessige Bloecke)."""
        K = self.far(rows, cols); m = self.m; k = self.k; X, W, Xo, Wo = self.X, self.W, self.Xo, self.Wo
        D = np.linalg.norm(m.cent[rows][:, None, :] - m.cent[cols][None, :, :], axis=-1)
        near = D < self.near_fac*np.maximum(m.h[rows][:, None], m.h[cols][None, :])
        from analytic import tri_integrals_multi
        for a in np.unique(np.nonzero(near)[0]):
            i = rows[a]; bs = np.nonzero(near[a])[0]; js = cols[bs]
            normal = (js == i) | (m.h[i] <= 1.5*m.h[js])
            # Faelle 1: aeusseres Integral ueber i (vektorisiert ueber j)
            b1 = bs[normal]; j1 = cols[b1]
            if len(j1):
                Ig, Ii = tri_integrals_multi(Xo[i], m.V[j1], m.n[j1]); wo = Wo[i]
                v = np.einsum('q,tqd->td', wo, Ig)/(4*np.pi); v[j1 == i] = 0
                s = -1j*k*np.einsum('q,tq->t', wo, Ii)/(4*np.pi)
                Z = Xo[i][None, :, None, :] - X[j1][:, None, :, :]; ww = wo[None, :, None]*W[j1][:, None, :]
                vr, sr = kernel_rem(Z, k)
                K[a, b1, 0] = s + np.einsum('tqp,tqp->t', ww, sr); K[a, b1, 1:] = v + np.einsum('tqp,tqpd->td', ww, vr)
            for b in bs[~normal]:
                j = cols[b]
                Ig, Ii = tri_integrals(Xo[j], m.V[i], m.n[i]); wo = Wo[j]
                v = -(wo@Ig)/(4*np.pi); s = -1j*k*(wo@Ii)/(4*np.pi)
                Z = X[i][:, None, :] - Xo[j][None, :, :]; ww = W[i][:, None]*wo[None, :]
                vr, sr = kernel_rem(Z, k)
                K[a, b, 0] = s + np.sum(ww*sr); K[a, b, 1:] = v + np.einsum('qp,qpd->d', ww, vr)
        return K
# ---------------- Clusterbaum und Blockpartition ----------------
class Cluster:
    def __init__(self, idx, cent, leaf=32):
        self.idx = idx; P = cent[idx]; self.lo, self.hi = P.min(0), P.max(0)
        self.center = 0.5*(self.lo + self.hi); self.diam = np.linalg.norm(self.hi - self.lo)
        self.children = []
        if len(idx) > leaf:
            ax = np.argmax(self.hi - self.lo); med = np.median(P[:, ax])
            L, R = idx[P[:, ax] <= med], idx[P[:, ax] > med]
            if len(L) and len(R): self.children = [Cluster(L, cent, leaf), Cluster(R, cent, leaf)]
def box_dist(a, b):
    d = np.maximum(0, np.maximum(a.lo - b.hi, b.lo - a.hi)); return np.linalg.norm(d)
def partition(t, s, mesh, eta, kmax_diam, out):
    dist = box_dist(t, s); hmax = max(mesh.h[t.idx].max(), mesh.h[s.idx].max())
    adm = min(t.diam, s.diam) <= eta*dist and dist > 3*hmax and max(t.diam, s.diam)*kmax_diam <= 20
    if adm: out.append(('lr', t, s)); return
    if not t.children or not s.children: out.append(('dense', t, s)); return
    for a in t.children:
        for b in s.children: partition(a, b, mesh, eta, kmax_diam, out)
# ---------------- ACA mit partieller Pivotisierung ----------------
def aca(row_fn, col_fn, nr, nc, eps, rmax=200):
    """row_fn(i) -> (nc,), col_fn(j) -> (nr,). Liefert U (nr,r), V (nc,r) mit A ~ U V^T."""
    U, V = [], []; used = set(); i = 0; normS2 = 0.0
    for it in range(min(rmax, nr, nc)):
        r = row_fn(i) - (sum(u[i]*v for u, v in zip(U, V)) if U else 0)
        used.add(i); j = int(np.argmax(np.abs(r)))
        if abs(r[j]) < 1e-300:
            cand = [q for q in range(nr) if q not in used]
            if not cand: break
            i = cand[0]; continue
        v = r/r[j]; u = col_fn(j) - (sum(uu*vv[j] for uu, vv in zip(U, V)) if U else 0)
        nu, nv = np.linalg.norm(u), np.linalg.norm(v)
        normS2 += (nu*nv)**2 + 2*sum(np.real(np.vdot(uu, u)*np.vdot(vv, v)) for uu, vv in zip(U, V))
        U.append(u); V.append(v)
        if nu*nv <= eps*np.sqrt(abs(normS2)): break
        uu = np.abs(u); uu[list(used)] = -1; i = int(np.argmax(uu))
    if not U: return np.zeros((nr, 0), complex), np.zeros((nc, 0), complex)
    return np.array(U).T, np.array(V).T
def recompress(U, V, eps):
    if U.shape[1] == 0: return U, V
    qu, ru = np.linalg.qr(U); qv, rv = np.linalg.qr(V)
    w, s, zh = np.linalg.svd(ru@rv.T); r = max(1, int(np.sum(s > eps*s[0])))
    return qu@(w[:, :r]*s[:r]), qv@zh[:r].T          # U V^T = qu (w s zh) qv^T
# ---------------- H-Matrix fuer E_k ----------------
class HMatrix:
    def __init__(self, mesh, k, eps=1e-4, eta=1.0, leaf=32, mode='joint', verbose=True):
        t0 = time.time(); self.m, self.k, self.mode = mesh, k, mode; N = len(mesh.T)
        self.ker = Kernel(mesh, k); root = Cluster(np.arange(N), mesh.cent, leaf)
        blocks = []; partition(root, root, mesh, eta, abs(k), blocks)
        self.blocks = []; self.mem = 0; self.mem_lr = 0; ranks = []
        for typ, t, s in blocks:
            R, C = t.idx, s.idx
            if typ == 'dense':
                K = self.ker.dense(R, C); self.blocks.append(('dense', R, C, K)); self.mem += K.size
            else:
                if mode == 'joint':
                    rf = lambda i, R=R, C=C: self.ker.far(R[i:i+1], C)[0].T.ravel()          # (4C,)
                    cf = lambda j, R=R, C=C: self.ker.far(R, C[j % len(C):j % len(C)+1])[:, 0, j//len(C)]
                    U, V = aca(rf, cf, len(R), 4*len(C), eps); U, V = recompress(U, V, eps)
                    self.blocks.append(('joint', R, C, (U, V))); self.mem += U.size + V.size; self.mem_lr += U.size + V.size; ranks.append(U.shape[1])
                else:
                    facs = []
                    for c in range(4):
                        rf = lambda i, R=R, C=C, c=c: self.ker.far(R[i:i+1], C)[0, :, c]
                        cf = lambda j, R=R, C=C, c=c: self.ker.far(R, C[j:j+1])[:, 0, c]
                        U, V = aca(rf, cf, len(R), len(C), eps); U, V = recompress(U, V, eps)
                        facs.append((U, V)); self.mem += U.size + V.size; self.mem_lr += U.size + V.size; ranks.append(U.shape[1])
                    self.blocks.append(('comp', R, C, facs))
        self.ranks = ranks; self.t_build = time.time() - t0
        # Hilfsgroessen fuer das Matrix-Vektor-Produkt: L(e_c n_j) / sqrt(A_j)
        from cl3 import MV
        basis = [np.eye(8)[0], np.eye(8)[1], np.eye(8)[2], np.eye(8)[4]]
        self.Ln = np.zeros((N, 4, 8, 8), complex)
        for c in range(4):
            s = np.zeros(N); v = np.zeros((N, 3))
            if c == 0: s[:] = 1
            else: v[:, c-1] = 1
            Mc = mv_times_n(s.astype(complex), v.astype(complex), mesh.n)
            self.Ln[:, c] = np.einsum('jb,bxy->jxy', Mc, LB)/np.sqrt(mesh.area)[:, None, None]
        self.isqA = 1/np.sqrt(mesh.area)
        if verbose:
            nd = sum(1 for b in self.blocks if b[0] == 'dense'); nl = len(self.blocks) - nd
            print(f"  H-Matrix k={k:.3g} ({mode}, eps={eps:g}): {nd} dichte / {nl} Niedrigrangbloecke, "
                  f"Speicher {self.mem*16/2**20:.1f} MB (voll: {N*N*4*16/2**20:.1f} MB Kern), mittl. Rang {np.mean(ranks) if ranks else 0:.1f}, "
                  f"{self.t_build:.0f}s", flush=True)
    def matvec(self, x):
        """y = E x, x: (8N,)"""
        N = len(self.m.T); X = x.reshape(N, 8)
        Z = np.einsum('jcxy,jy->jcx', self.Ln, X)            # (N,4,8)
        Y = np.zeros((N, 8), complex)
        for typ, R, C, dat in self.blocks:
            if typ == 'dense': Y[R] += np.einsum('ijc,jcx->ix', dat, Z[C])
            elif typ == 'joint':
                U, V = dat; zs = Z[C].transpose(1, 0, 2).reshape(4*len(C), 8)   # Reihenfolge (c, j) wie V
                Y[R] += U@(V.T@zs)
            else:
                for c, (U, V) in enumerate(dat): Y[R] += U@(V.T@Z[C, c])
        return (-2*self.isqA[:, None]*Y).ravel()
