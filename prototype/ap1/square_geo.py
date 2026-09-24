"""Geometrisch gradiertes Netz zu den Ecken (Verhaeltnis q), gewichteter Ansatz h = r^{-beta} g.
Damit sind minimale Elementgroessen bis 1e-9 erreichbar: Unterscheidung von L^2 und L^2(r^{2 beta})."""
import numpy as np
from square_bem import E1m, E2m, vec_mat
from square_weighted import seg_integral
from corner import Jblock

def side_nodes(ng, q=0.5, nmid=4, q0=0.5):
    """q0: Beginn der geometrischen Zone (Abstand zur Ecke)."""
    s = [0.0] + [q0*q**k for k in range(ng, -1, -1)]          # 0, q0 q^ng, ..., q0
    mid = list(np.linspace(q0, 2 - q0, nmid+1)[1:-1])
    right = [2 - x for x in s[::-1]]
    return np.array(s + mid + right)

def build_geo(ng, beta, kappa, q=0.5, a=None, b=None, m=5, K=40, beta_corner=None, nmid=4, q0=0.5, kappa_fn=None):
    """beta_corner (komplex): Ansatz r^{-beta_corner} nur in den Eckelementen (Anreicherung, z.B. 1 - a_j)."""
    t = side_nodes(ng, q, nmid, q0); M = len(t) - 1
    V = [np.array(v, float) for v in [(1, -1), (1, 1), (-1, 1), (-1, -1)]]
    P0, P1, Nr, T0, T1, SD = [], [], [], [], [], []
    for k in range(4):
        A, B = V[k], V[(k+1) % 4]
        u = (B - A)/2; n = np.array([u[1], -u[0]])
        for i in range(M):
            Nr.append(n); T0.append(t[i]); T1.append(t[i+1]); SD.append(k)
    Nr, T0, T1, SD = map(np.array, (Nr, T0, T1, SD))
    Ns = len(T0); L = T1 - T0; tm = 0.5*(T0 + T1)
    U = np.array([(V[(k+1) % 4] - V[k])/2 for k in range(4)])
    # Kollokationspunkte exakt als (naechste Ecke, Versatz)
    Kc = np.array([V[SD[i]] if tm[i] <= 1 else V[(SD[i]+1) % 4] for i in range(Ns)])
    Oc = np.array([tm[i]*U[SD[i]] if tm[i] <= 1 else -(2 - tm[i])*U[SD[i]] for i in range(Ns)])
    rX = np.minimum(tm, 2 - tm); X = Kc + Oc
    def Fb(s, bt):
        s = np.asarray(s, complex)
        return np.where(s.real <= 1, s**(1-bt)/(1-bt), 2/(1-bt) - (2-s)**(1-bt)/(1-bt))
    corner_el = (T0 == 0.0) | (T1 == 2.0)
    bet = np.where(corner_el, beta if beta_corner is None else beta_corner, beta).astype(complex)
    Ew = np.zeros((8*Ns, 8*Ns), complex)
    for j in range(Ns):
        t0, t1 = T0[j], T1[j]; hj = t1 - t0; side = SD[j]; u = U[side]
        F = lambda s, bt=bet[j]: Fb(s, bt)
        if t0 == 0.0 or t1 == 2.0:
            d = 0.1*hj
            sl = sorted(set([0.0] + [(0.5*hj - d)*0.5**k for k in range(K)] + [0.5*hj - d, 0.5*hj + d, 0.75*hj, hj]))
            sl = np.array(sl); nodes = t0 + sl if t0 == 0.0 else t1 - sl[::-1]
        else:
            nodes = np.linspace(t0, t1, m+1)
        O = V[side] if t0 < 1 else V[(side+1) % 4]
        Xrel = (Kc - O) + Oc
        Vsum = np.zeros((Ns, 2), complex)
        for s0, s1 in zip(nodes[:-1], nodes[1:]):
            if not s1 > s0: continue
            y0 = s0*u if t0 < 1 else -(2 - s0)*u
            y1 = s1*u if t0 < 1 else -(2 - s1)*u
            if t0 == 0.0: y0, y1 = s0*u, s1*u
            if t1 == 2.0: y0, y1 = -(2 - s0)*u, -(2 - s1)*u
            Vsum += (F(s1) - F(s0))/(s1 - s0)*seg_integral(Xrel, y0, y1)
        Mb = -2*(Vsum[:, 0, None, None]*E1m + Vsum[:, 1, None, None]*E2m)@vec_mat(Nr[j])
        Ew[:, 8*j:8*j+8] = Mb.reshape(8*Ns, 8)
    Dr = np.repeat(rX**(-bet), 8)
    aa = 1.0 if a is None else a; bb = np.sqrt(kappa) if b is None else b
    Jb = np.zeros_like(Ew)
    for j in range(Ns):
        kj = kappa if kappa_fn is None else kappa_fn(rX[j], kappa)
        bj = (np.sqrt(kj) if b is None else b) if kappa_fn is not None else bb
        Jb[8*j:8*j+8, 8*j:8*j+8] = Jblock(np.array([Nr[j][0], Nr[j][1], 0.]), kj, 1.0, aa, bj)
    Bm = 0.5*(np.diag(Dr) + Ew) + 0.5*(np.diag(Dr) - Ew)@Jb
    rs = np.repeat(rX**bet*np.sqrt(L), 8); cs = np.repeat(1/np.sqrt(L), 8)
    build_geo.last = dict(Ew=Ew, X=X, L=L, rX=rX, bet=bet, T0=T0, T1=T1)
    return rs[:, None]*Bm*cs[None, :], X, np.repeat(rX, 8)

if __name__ == '__main__':
    import time, json, sys
    out = {}
    for ng in [4, 8, 16, 24, 32]:
        row = []
        for name, k in [('innen', -3+0.3j), ('-4', -4+0j), ('Gold', -11+1.2j)]:
            for beta in [0.0, 0.35, 0.45]:
                A, _, _ = build_geo(ng, beta, k)
                s = np.linalg.svd(A, compute_uv=False)[-1]; out[f'{name}|{beta}|{ng}'] = s
                row.append(f'{name} b={beta}: {s:.4f}')
        print(f'ng={ng:2d} (h_min={0.5*0.5**ng:.1e}): ' + ', '.join(row), flush=True)
    json.dump(out, open('square_geo.json', 'w'))
