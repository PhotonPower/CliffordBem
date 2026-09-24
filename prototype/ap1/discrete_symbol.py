"""Diskretes Symbol der Kollokation auf dem unendlichen geometrischen Eckgitter (Verhaeltnis q):
T ist block-Toeplitz im Schichtindex, a(z) = sum_k T_{0,k} z^k (16x16). Nullstellen im Ring
1 < |z| < 1/q entsprechen diskreten Exponenten a_h = 1 + ln z / ln q  (kontinuierlich: z = q^{a-1})."""
import numpy as np
from square_weighted import seg_integral
from square_bem import E1m, E2m, vec_mat
from corner import Jblock
def blocks(kappa, q=0.5, K1=70, K2=40, alpha=np.pi/2):
    d = [np.array([1., 0]), np.array([np.cos(alpha), np.sin(alpha)])]
    n = [np.array([0., -1]), np.array([-np.sin(alpha), np.cos(alpha)])]
    Xc = [0.5*(1 + q)*d[i] for i in range(2)]                      # Kollokation: Mitte von [q,1]
    Jb = [Jblock(np.array([n[i][0], n[i][1], 0.]), kappa, 1.0, 1.0, np.sqrt(kappa)) for i in range(2)]
    B = {}
    for k in range(-K1, K2+1):
        Tk = np.zeros((16, 16), complex)
        for j in range(2):                                          # Ansatzelement [q^{k+1}, q^k] auf Strahl j
            y0, y1 = q**(k+1)*d[j], q**k*d[j]
            V = seg_integral(np.array(Xc), y0, y1)                  # (2,2)
            for i in range(2):
                Eij = -2*(V[i, 0]*E1m + V[i, 1]*E2m)@vec_mat(n[j])
                Id = np.eye(8) if (i == j and k == 0) else np.zeros((8, 8))
                Tk[8*i:8*i+8, 8*j:8*j+8] = 0.5*(Id + Eij) + 0.5*(Id - Eij)@Jb[j]
        B[k] = Tk
    return B
def symbol(B, z):
    return sum(Tk*z**k for k, Tk in B.items())
def scan(kappa, q=0.5, nr=48, nphi=96):
    B = blocks(kappa, q); rs = np.linspace(1.03, 1/q - 0.02, nr); ph = np.linspace(-np.pi, np.pi, nphi, endpoint=False)
    D = np.array([[abs(np.linalg.det(symbol(B, r*np.exp(1j*p)))) for p in ph] for r in rs])
    return rs, ph, D, B
if __name__ == '__main__':
    from scipy.optimize import minimize
    for kappa in [-2.8+0.3j, -1.5+0.2j, -11+1.2j]:
        rs, ph, D, B = scan(kappa)
        # Normierung: det auf jedem Kreis relativ zum Maximum
        cand = []
        for i in range(1, len(rs)-1):
            for j in range(len(ph)):
                v = D[i, j]
                if v < D[i-1, j] and v < D[i+1, j] and v < D[i, j-1] and v < D[i, (j+1) % len(ph)]:
                    cand.append((rs[i], ph[j], v/np.max(D[i])))
        res = []
        for r0, p0, _ in cand:
            f = lambda v: abs(np.linalg.det(symbol(B, v[0]*np.exp(1j*v[1]))))
            o = minimize(f, [r0, p0], method='Nelder-Mead', options={'xatol': 1e-8, 'fatol': 1e-14})
            z = o.x[0]*np.exp(1j*o.x[1]); rel = o.fun/np.max(np.abs([np.linalg.det(symbol(B, o.x[0]*np.exp(1j*t))) for t in ph]))
            if rel < 1e-4 and 1.0 < o.x[0] < 2.0: res.append((z, 1 + np.log(z)/np.log(0.5), rel))
        print(f"kappa={kappa}: diskrete Nullstellen (z, a_h, rel|det|):")
        for z, ah, rel in res: print(f"     z={z:.4f}  a_h={ah:.4f}  {rel:.1e}")

def blocks_w(kappa, beta, q=0.5, K1=90, K2=50, m=25, alpha=np.pi/2):
    """wie blocks, aber Ansatz (r/q^k)^{-beta} auf Element k (skaleninvariant normiert)."""
    d = [np.array([1., 0]), np.array([np.cos(alpha), np.sin(alpha)])]
    n = [np.array([0., -1]), np.array([-np.sin(alpha), np.cos(alpha)])]
    rc = 0.5*(1 + q); Xc = np.array([rc*d[i] for i in range(2)])
    Jb = [Jblock(np.array([n[i][0], n[i][1], 0.]), kappa, 1.0, 1.0, np.sqrt(kappa)) for i in range(2)]
    B = {}
    for k in range(-K1, K2+1):
        Tk = np.zeros((16, 16), complex); s0 = q**(k+1); s1 = q**k
        nodes = np.linspace(s0, s1, m+1)
        for j in range(2):
            V = np.zeros((2, 2))
            for a0, a1 in zip(nodes[:-1], nodes[1:]):
                avg = ((a1/s1)**(1-beta) - (a0/s1)**(1-beta))*s1/((1-beta)*(a1 - a0))
                V = V + avg*seg_integral(Xc, a0*d[j], a1*d[j])
            for i in range(2):
                Eij = -2*(V[i, 0]*E1m + V[i, 1]*E2m)@vec_mat(n[j])
                Id = (rc/1.0)**(-beta)*np.eye(8) if (i == j and k == 0) else np.zeros((8, 8))
                Tk[8*i:8*i+8, 8*j:8*j+8] = 0.5*(Id + Eij) + 0.5*(Id - Eij)@Jb[j]
        B[k] = Tk
    return B

def zeros_annulus(B, q, nr=36, nph=120):
    from scipy.optimize import minimize
    zmax = 1/q; rs = np.linspace(1 + 0.15*(zmax - 1), zmax*0.995, nr); ph = np.linspace(-np.pi, np.pi, nph, endpoint=False)
    D = np.array([[abs(np.linalg.det(symbol(B, r*np.exp(1j*p)))) for p in ph] for r in rs]); res = []
    for i in range(1, nr-1):
        for j in range(nph):
            v = D[i, j]
            if v < D[i-1, j] and v < D[i+1, j] and v < D[i, j-1] and v < D[i, (j+1) % nph]:
                o = minimize(lambda w: abs(np.linalg.det(symbol(B, w[0]*np.exp(1j*w[1])))), [rs[i], ph[j]], method='Nelder-Mead', options={'xatol': 1e-9, 'fatol': 1e-15})
                z = o.x[0]*np.exp(1j*o.x[1]); ah = 1 + np.log(z)/np.log(q)
                if o.fun/np.max(D[i]) < 1e-5 and 0 < ah.real < 1 and not any(abs(ah - x) < 1e-4 for x in res): res.append(ah)
    return sorted(res, key=lambda x: x.real)

def symbol_tail(B, z, q=0.5):
    """Symbol mit analytischen Schwaenzen: T_{0,k} ~ T_{0,K2} q^{k-K2} (k -> +oo), ~ T_{0,-K1} (k -> -oo)."""
    ks = sorted(B); K1, K2 = -ks[0], ks[-1]
    S = sum(Tk*z**k for k, Tk in B.items())
    S = S + B[K2]*z**K2*(q*z)/(1 - q*z) + B[-K1]*z**(-K1)*(1/z)/(1 - 1/z)
    return S

def zeros_tail(B, q, nr=40, nph=120, rmax_frac=0.999):
    from scipy.optimize import minimize
    zmax = 1/q; rs = np.linspace(1 + 0.05*(zmax - 1), zmax*rmax_frac, nr); ph = np.linspace(-np.pi, np.pi, nph, endpoint=False)
    f = lambda z: abs(np.linalg.det(symbol_tail(B, z, q)))
    D = np.array([[f(r*np.exp(1j*p)) for p in ph] for r in rs]); res = []
    for i in range(1, nr-1):
        for j in range(nph):
            v = D[i, j]
            if v < D[i-1, j] and v < D[i+1, j] and v < D[i, j-1] and v < D[i, (j+1) % nph]:
                o = minimize(lambda w: f(w[0]*np.exp(1j*w[1])), [rs[i], ph[j]], method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-16})
                z = o.x[0]*np.exp(1j*o.x[1]); ah = 1 + np.log(z)/np.log(q)
                if o.fun/np.max(D[i]) < 1e-6 and 0 < ah.real < 1 and not any(abs(ah - x) < 1e-4 for x in res): res.append(ah)
    return sorted(res, key=lambda x: x.real)
