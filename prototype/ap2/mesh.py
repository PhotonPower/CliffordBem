"""Dreiecksnetze: Ikosaeder-Kugel mit n-facher Kantenunterteilung (20 n^2 Dreiecke), aussen orientiert."""
import numpy as np
def icosahedron():
    t = (1 + 5**0.5)/2
    V = np.array([[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
                  [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]], float)
    F = np.array([[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
                  [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]])
    return V/np.linalg.norm(V, axis=1)[:, None], F
def icosphere(n, R=1.0):
    V0, F0 = icosahedron(); verts = {}; P = []; T = []
    def vid(p):
        key = tuple(np.round(p, 12))
        if key not in verts: verts[key] = len(P); P.append(p)
        return verts[key]
    for f in F0:
        A, B, C = V0[f]
        idx = {}
        for i in range(n+1):
            for j in range(n+1-i):
                p = A + (B - A)*i/n + (C - A)*j/n; p = p/np.linalg.norm(p)*R
                idx[(i, j)] = vid(p)
        for i in range(n):
            for j in range(n-i):
                T.append([idx[(i, j)], idx[(i+1, j)], idx[(i, j+1)]])
                if j < n-i-1: T.append([idx[(i+1, j)], idx[(i+1, j+1)], idx[(i, j+1)]])
    P = np.array(P); T = np.array(T)
    # Orientierung nach aussen
    c = P[T].mean(1); nrm = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
    flip = np.sum(nrm*c, 1) < 0; T[flip] = T[flip][:, [0, 2, 1]]
    return P, T
class Mesh:
    def __init__(self, P, T):
        self.P, self.T = P, T
        a, b, c = P[T[:, 0]], P[T[:, 1]], P[T[:, 2]]
        cr = np.cross(b - a, c - a); self.area = 0.5*np.linalg.norm(cr, axis=1)
        self.n = cr/(2*self.area[:, None]); self.cent = (a + b + c)/3
        self.h = np.max([np.linalg.norm(b - a, axis=1), np.linalg.norm(c - b, axis=1), np.linalg.norm(a - c, axis=1)], 0)
        self.V = np.stack([a, b, c], 1)
# Dunavant-Regel Grad 5 (7 Punkte) in baryzentrischen Koordinaten
_a1, _b1 = 0.059715871789770, 0.470142064105115
_a2, _b2 = 0.797426985353087, 0.101286507323456
BARY = np.array([[1/3, 1/3, 1/3], [_a1, _b1, _b1], [_b1, _a1, _b1], [_b1, _b1, _a1], [_a2, _b2, _b2], [_b2, _a2, _b2], [_b2, _b2, _a2]])
WQ = np.array([0.225, 0.132394152788506, 0.132394152788506, 0.132394152788506, 0.125939180544827, 0.125939180544827, 0.125939180544827])
def quad_points(mesh, sub=1):
    """Quadraturpunkte (N, q, 3) und Gewichte (N, q) inkl. Flaeche; sub>1: Unterteilung in sub^2 Teildreiecke."""
    if sub == 1:
        X = np.einsum('qk,nkd->nqd', BARY, mesh.V); W = WQ[None, :]*mesh.area[:, None]
        return X, W
    B = []; Ws = []
    for i in range(sub):
        for j in range(sub - i):
            for up in [0, 1]:
                if up == 0: tri = np.array([[i, j], [i+1, j], [i, j+1]])/sub
                elif j < sub - i - 1: tri = np.array([[i+1, j], [i+1, j+1], [i, j+1]])/sub
                else: continue
                for bq, wq in zip(BARY, WQ):
                    uv = bq@tri; B.append([1 - uv.sum(), uv[0], uv[1]]); Ws.append(wq/sub**2)
    B = np.array(B); Ws = np.array(Ws)
    X = np.einsum('qk,nkd->nqd', B, mesh.V); W = Ws[None, :]*mesh.area[:, None]
    return X, W
if __name__ == '__main__':
    for n in [2, 4, 6]:
        m = Mesh(*icosphere(n)); print(n, len(m.T), 'Dreiecke; Flaeche', m.area.sum(), '(Kugel 4pi =', 4*np.pi, '); h_max', m.h.max(),
                                       '; Normalen aussen:', np.all(np.sum(m.n*m.cent, 1) > 0))
    X, W = quad_points(m, 3); print('Quadratur (sub=3): Summe der Gewichte / Flaeche =', W.sum()/m.area.sum())
