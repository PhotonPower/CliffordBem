"""Wuerfel [-1,1]^3 mit geometrisch zu den Kanten gradiertem, spiegelsymmetrischem Dreiecksnetz.
Knoten je Halbseite: 0, 0.5, 1-2^-2, ..., 1-2^-L, 1  ->  48 (L+1)^2 Dreiecke."""
import numpy as np
from mesh import Mesh
def nodes_1d(L):
    half = [0.0, 0.5] + [1 - 0.5**k for k in range(2, L+1)] + [1.0]
    half = sorted(set(half))
    return np.array([-x for x in half[::-1]] + half[1:])
def cube_mesh(L):
    s = nodes_1d(L); P = []; T = []; key = {}
    def vid(p):
        kk = tuple(np.round(p, 13))
        if kk not in key: key[kk] = len(P); P.append(np.array(p, float))
        return key[kk]
    for ax in range(3):
        for sg in [1, -1]:
            u_ax, v_ax = [a for a in range(3) if a != ax]
            for i in range(len(s)-1):
                for j in range(len(s)-1):
                    def pt(a, b):
                        p = np.zeros(3); p[ax] = sg; p[u_ax] = a; p[v_ax] = b; return vid(p)
                    u0, u1, v0, v1 = s[i], s[i+1], s[j], s[j+1]
                    uc, vc = 0.5*(u0 + u1), 0.5*(v0 + v1)
                    A, B, C, D = pt(u0, v0), pt(u1, v0), pt(u1, v1), pt(u0, v1)
                    if np.sign(uc)*np.sign(vc) > 0: T += [[A, B, C], [A, C, D]]       # Diagonale A-C
                    else:                           T += [[A, B, D], [B, C, D]]       # Diagonale B-D (gespiegelt)
    P = np.array(P); T = np.array(T)
    c = P[T].mean(1); nrm = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
    flip = np.sum(nrm*c, 1) < 0; T[flip] = T[flip][:, [0, 2, 1]]
    return Mesh(P, T)
if __name__ == '__main__':
    for L in range(1, 7):
        m = cube_mesh(L); print(L, len(m.T), 'Dreiecke, Flaeche', round(m.area.sum(), 12), '(24), min. Kantenabstand', 0.5**L if L > 1 else 0.5)
