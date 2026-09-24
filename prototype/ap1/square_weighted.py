"""Gewichteter Ansatzraum: h = r^{-beta} g, g stueckweise konstant. Diskretisiert wird
A_beta = r^beta T r^{-beta} auf L^2 (Kollokation, exakte Mittelung von r^{-beta} auf Teilsegmenten,
geometrische Unterteilung in Eckelementen). Stabilitaet sollte vom Mellin-Symbol auf
Re a = 1/2 - beta bestimmt sein."""
import numpy as np
from square_bem import square_mesh, E1m, E2m, vec_mat, T1_matrix
from corner import Jblock

def seg_integral(X, p0, p1):
    """int_seg Phi2(x - y) dl(y) fuer alle x in X (p.v., falls x im Segmentinneren)."""
    Lv = np.linalg.norm(p1 - p0); u = (p1 - p0)/Lv; up = np.array([-u[1], u[0]])
    w = X - p0; a = w@u; b = w@up
    on = (np.abs(b) < 1e-13) & (a > 0) & (a < Lv)
    with np.errstate(divide='ignore', invalid='ignore'):
        Iu = -0.5*(np.log((a - Lv)**2 + b**2) - np.log(a**2 + b**2))
        # Oeffnungswinkel robust (ohne Verzweigungsschnitt): Winkel zwischen x-p0 und x-p1
        v0 = X - p0; v1 = X - p1
        cr = v0[:, 0]*v1[:, 1] - v0[:, 1]*v1[:, 0]; dt = np.sum(v0*v1, 1)
        Ip = np.arctan2(cr*np.sign(u[0]*up[1] - u[1]*up[0]), dt)
    Ip = np.where(on, 0.0, Ip)
    return (Iu[:, None]*u[None, :] + Ip[:, None]*up[None, :])/(2*np.pi)

def subsegments(t0, t1, h, corner_left, corner_right, m=5, K=36):
    """Unterteilung von [t0,t1] (Bogenlaenge auf der Seite, Seitenlaenge 2) und Mittelwerte von r^{-beta}."""
    mid = 0.5*(t0 + t1)
    if corner_left or corner_right:
        # lokale Koordinate s = Abstand zur Ecke
        d = 0.1*(t1 - t0)
        g = [(0.5*(t1-t0) - d)*0.5**k for k in range(K)] + [0.0]
        s_nodes = sorted(set([0.0] + g[:-1] + [0.5*(t1-t0) - d, 0.5*(t1-t0) + d, 0.75*(t1-t0), (t1-t0)]))
        s_nodes = np.array(s_nodes)
        if corner_left: nodes = t0 + s_nodes
        else: nodes = t1 - s_nodes[::-1]
    else:
        nodes = np.linspace(t0, t1, m+1)
    return np.array(nodes)

def build(N, beta, kappa, a=None, b=None):
    P0, P1, Nr, V = square_mesh(N)
    Ns = len(P0); X = 0.5*(P0 + P1); L = np.linalg.norm(P1 - P0, axis=1)
    h = 2.0/N
    rX = np.min([np.linalg.norm(X - v, axis=1) for v in V], axis=0)
    Ew = np.zeros((8*Ns, 8*Ns), complex)
    for j in range(Ns):
        side = j // N; i_loc = j % N
        A = np.array([(1, -1), (1, 1), (-1, 1), (-1, -1)][side], float)
        t0, t1 = i_loc*h, (i_loc+1)*h
        nodes = subsegments(t0, t1, h, i_loc == 0, i_loc == N-1)
        tdir = (P1[j] - P0[j])/L[j]
        Vsum = np.zeros((Ns, 2))
        for s0, s1 in zip(nodes[:-1], nodes[1:]):
            # Mittelwert von r^{-beta}, r = min(t, 2 - t)
            def F(t):  # Stammfunktion von r^{-beta} entlang der Seite
                return np.where(t <= 1, t**(1-beta)/(1-beta), 2/(1-beta) - (2-t)**(1-beta)/(1-beta))
            avg = (F(s1) - F(s0))/(s1 - s0)
            Vsum += avg*seg_integral(X, A + s0*tdir, A + s1*tdir)
        M = -2*(Vsum[:, 0, None, None]*E1m + Vsum[:, 1, None, None]*E2m)@vec_mat(Nr[j])
        Ew[:, 8*j:8*j+8] = M.reshape(8*Ns, 8)
    # Diagonale: r^{-beta}(x_i) (Punktwert der Ansatzfunktion am Kollokationspunkt)
    Dr = np.repeat(rX**(-beta), 8)
    Jb = np.zeros_like(Ew)
    aa = 1.0 if a is None else a; bb = np.sqrt(kappa) if b is None else b
    for j in range(Ns):
        Jb[8*j:8*j+8, 8*j:8*j+8] = Jblock(np.array([Nr[j][0], Nr[j][1], 0.]), kappa, 1.0, aa, bb)
    Bm = 0.5*(np.diag(Dr) + Ew) + 0.5*(np.diag(Dr) - Ew)@Jb
    # A_beta: Zeilen mit r^beta sqrt(L), Spalten mit 1/sqrt(L)  (g in L^2)
    rs = np.repeat(rX**beta*np.sqrt(L), 8); cs = np.repeat(1/np.sqrt(L), 8)
    return rs[:, None]*Bm*cs[None, :]

if __name__ == '__main__':
    import sys, time
    # Kontrolle: beta = 0 muss die ungewichtete Matrix reproduzieren (bis auf Teilsegment-Quadratur)
    from square_bem import E_matrix
    P0, P1, Nr, V = square_mesh(16); E, X, L = E_matrix(P0, P1, Nr)
    T = T1_matrix(E, Nr, -3+0.3j); d = np.repeat(np.sqrt(L), 8)
    A0 = build(16, 0.0, -3+0.3j)
    print("beta=0 Abgleich:", np.abs(A0 - d[:, None]*T/d[None, :]).max())
