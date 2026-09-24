"""Schritt 1.5: Mellin-Symbol des Cauchy-Randoperators an einer 2D-Ecke (Konormalsymbol
einer 3D-Kante). Rand = zwei Strahlen vom Ursprung mit Richtungen d1, d2 (Oeffnungswinkel
alpha des Gebiets Omega = Medium 1). Dichte auf Strahl j:  h_j(t) = t^lam c_j,
lam = -1/2 + i tau  (L^2-kritische Gerade). E h = s^lam M(lam) c."""
import numpy as np
from cl3 import MV

def left(A): return A.left_matrix()

def geometry(alpha):
    d = [np.array([1., 0, 0]), np.array([np.cos(alpha), np.sin(alpha), 0])]
    n = [np.array([0., -1, 0]), np.array([-np.sin(alpha), np.cos(alpha), 0])]
    return d, n

def Phi2(z):   # 2D-Cauchy-Kern z/(2 pi |z|^2) als Vektor in der Ebene
    return z/(2*np.pi*np.sum(z**2, -1, keepdims=True))

_x = np.linspace(-70, 70, 28001); _dx = _x[1]-_x[0]
def mellin_cross(di, dj, nj, tau, c=0.5):
    """ -2 * int_0^inf Phi2(di - u dj) nj u^lam du  als 8x8 (Linksmultiplikation);
    Mellin-Gerade lam + 1 = c + i tau (c = 1/2: L^2; c = 1/2 - beta: Gewicht r^{2 beta})."""
    u = np.exp(_x)
    z = di[None, :] - u[:, None]*dj[None, :]
    P = Phi2(z)                                              # (N,3)
    w = np.exp(_x*c + 1j*tau*_x)*_dx
    v = np.sum(P*w[:, None], 0)                              # komplexer Vektor
    return -2*left(MV.vec(v))@left(MV.vec(nj))

def M_E(alpha, tau, c=0.5):
    d, n = geometry(alpha)
    M = np.zeros((16, 16), complex)
    for i in range(2):
        for j in range(2):
            if i == j:
                # -(1/pi) pv int_0^inf u^lam/(1-u) du = -cot(pi a),  a = c + i tau
                B = -1/np.tan(np.pi*(c + 1j*tau))*left(MV.vec(d[i]))@left(MV.vec(n[i]))
            else:
                B = mellin_cross(d[i], d[j], n[j], tau, c)
            M[8*i:8*i+8, 8*j:8*j+8] = B
    return M

def Jblock(n, re_, rm, a, b, inv=False):
    """J als 8x8-Matrix zur Normalen n (in der Ebene)."""
    se, sm = np.sqrt(re_), np.sqrt(rm)
    if inv: se, sm, a, b = 1/se, 1/sm, 1/a, 1/b
    Jm = np.zeros((8, 8), complex)
    for bb in range(8):
        H = np.zeros(8, complex); H[bb] = 1
        v = H[[1, 2, 4]]; Ib = H[[6, 5, 3]]*np.array([1, -1, 1])
        vn, bn = n@v, n@Ib
        v2 = se*(v - vn*n) + vn*n/se; b2 = sm*(Ib - bn*n) + bn*n/sm
        out = np.zeros(8, complex); out[0] = a*H[0]; out[7] = b*H[7]
        out[[1, 2, 4]] = v2; out[[6, 5, 3]] = b2*np.array([1, -1, 1])
        Jm[:, bb] = out
    return Jm

def Jfull(alpha, re_, rm=1.0, a=None, b=None, inv=False):
    a = np.sqrt(rm) if a is None else a; b = np.sqrt(re_) if b is None else b
    _, n = geometry(alpha)
    J = np.zeros((16, 16), complex)
    for i in range(2): J[8*i:8*i+8, 8*i:8*i+8] = Jblock(n[i], re_, rm, a, b, inv)
    return J

def T1_symbol(alpha, tau, re_, rm=1.0, a=None, b=None, ME=None):
    M = M_E(alpha, tau) if ME is None else ME
    I = np.eye(16)
    return 0.5*(I + M) + 0.5*(I - M)@Jfull(alpha, re_, rm, a, b)

if __name__ == '__main__':
    for alpha in [np.pi/2, np.pi, 3*np.pi/2, np.pi/3]:
        errs = [np.abs(M_E(alpha, t)@M_E(alpha, t) - np.eye(16)).max() for t in [0, 0.3, 1.0, 3.0]]
        print(f"alpha={alpha/np.pi:.3f} pi:  max |M(tau)^2 - 1| = {max(errs):.1e}")
