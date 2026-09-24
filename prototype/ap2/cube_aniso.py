"""Nichtkonformes, anisotrop gradiertes Wuerfelnetz: je Seitenflaeche ein grobes Innenquadrat,
vier Kantenstreifen (nur quer zur Kante geometrisch gradiert, entlang der Kante n_a Elemente)
und vier Eckquadrate (in beiden Richtungen gradiert). Haengende Knoten sind fuer stueckweise
konstante L^2-Ansaetze zulaessig."""
import numpy as np
from mesh import Mesh
def graded(L, w=0.5):
    """Knoten in [0, w] zur 0 hin geometrisch: 0, w 2^-(L-1), ..., w/2, w."""
    return np.array([0.0] + [w*0.5**k for k in range(L-1, -1, -1)])
def face_rects(L, na, nc, w=0.5):
    """Rechtecke (u0,u1,v0,v1) auf [-1,1]^2."""
    g = graded(L, w)                                  # Abstand zur Kante
    lo = -1 + g; hi = 1 - g[::-1]                     # zu -1 bzw. +1 hin gradiert
    mid = np.linspace(-1 + w, 1 - w, nc + 1)
    along = np.linspace(-1 + w, 1 - w, na + 1)
    R = []
    for a in range(nc):
        for b in range(nc): R.append((mid[a], mid[a+1], mid[b], mid[b+1]))
    for side in [lo, hi]:
        for a in range(len(side) - 1):
            for b in range(na):
                R.append((side[a], side[a+1], along[b], along[b+1]))   # Streifen an u = +-1
                R.append((along[b], along[b+1], side[a], side[a+1]))   # Streifen an v = +-1
    for su in [lo, hi]:
        for sv in [lo, hi]:
            for a in range(len(su) - 1):
                for b in range(len(sv) - 1): R.append((su[a], su[a+1], sv[b], sv[b+1]))
    return R
def cube_aniso(L, na=4, nc=4):
    P = []; T = []
    for ax in range(3):
        for sg in [1, -1]:
            ua, va = [x for x in range(3) if x != ax]
            for u0, u1, v0, v1 in face_rects(L, na, nc):
                def pt(u, v):
                    p = np.zeros(3); p[ax] = sg; p[ua] = u; p[va] = v; P.append(p); return len(P) - 1
                A, B, C, D = pt(u0, v0), pt(u1, v0), pt(u1, v1), pt(u0, v1)
                T += [[A, B, C], [A, C, D]]
    P = np.array(P); T = np.array(T)
    c = P[T].mean(1); nrm = np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]])
    flip = np.sum(nrm*c, 1) < 0; T[flip] = T[flip][:, [0, 2, 1]]
    return Mesh(P, T)
if __name__ == '__main__':
    for L in [2, 3, 4, 5, 6]:
        m = cube_aniso(L); print(L, len(m.T), 'Dreiecke; Flaeche', round(m.area.sum(), 10), '; kleinstes h', m.h.min().round(5))
