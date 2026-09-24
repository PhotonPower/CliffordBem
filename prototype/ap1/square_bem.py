"""2D-Test der Kanten-Vorhersage: Kollokation (stueckweise konstant, Mittelpunkte) von
T_1 = 1/2(1+E) + 1/2(1-E)J  im quasistatischen Grenzfall (k=0) auf dem Rand eines Quadrats.
E h(x) = -2 p.v. int Phi2(x-y) n(y) h(y) dl(y),  Phi2(z) = z/(2 pi |z|^2).
Vergleich: L^2 (beta = 0) gegen gewichtetes L^2(r^{2 beta}) via W T W^{-1}, W = r^beta."""
import numpy as np
from cl3 import MV, SIGN
from corner import Jblock

def square_mesh(N, grade=1.0):
    """N Segmente je Seite; grade > 1: symmetrische Graduierung zu den Ecken hin."""
    s = np.linspace(-1, 1, N+1)
    s = np.sign(s)*np.abs(s)**(1.0)  # gleichmaessig in [-1,1]
    if grade != 1.0:                   # Graduierung zu beiden Enden
        u = np.linspace(0, 1, N+1); s = -1 + 2*np.where(u < .5, 0.5*(2*u)**grade, 1 - 0.5*(2*(1-u))**grade)
    V = [np.array(v, float) for v in [(1, -1), (1, 1), (-1, 1), (-1, -1)]]
    P0, P1, Nr = [], [], []
    for k in range(4):
        A, B = V[k], V[(k+1) % 4]; t = (B - A)/2; n = np.array([t[1], -t[0]])   # aussen (gegen den Uhrzeigersinn)
        for i in range(N):
            P0.append(A + (s[i]+1)*t); P1.append(A + (s[i+1]+1)*t); Nr.append(n)
    return np.array(P0), np.array(P1), np.array(Nr), V

def vec_mat(v):   # Linksmultiplikation mit Vektor (v1, v2, 0)
    return MV.vec([v[0], v[1], 0]).left_matrix()
E1m, E2m = vec_mat([1, 0]), vec_mat([0, 1])

def E_matrix(P0, P1, Nr):
    Ns = len(P0); X = 0.5*(P0 + P1); L = np.linalg.norm(P1 - P0, axis=1)
    U = (P1 - P0)/L[:, None]; Up = np.stack([-U[:, 1], U[:, 0]], 1)
    E = np.zeros((8*Ns, 8*Ns), complex)
    Nm = [vec_mat(n) for n in Nr]
    for j in range(Ns):
        w = X - P0[j]; a = w@U[j]; b = w@Up[j]
        with np.errstate(divide='ignore', invalid='ignore'):
            Iu = -0.5*(np.log((a - L[j])**2 + b**2) - np.log(a**2 + b**2))       # int (a-t)/(...) dt
            Ip = -(np.arctan2(a - L[j], b) - np.arctan2(a, b))                    # int b/(...) dt
        Iu[j] = 0.0; Ip[j] = 0.0                                                   # p.v. am eigenen Segment
        V = (Iu[:, None]*U[j][None, :] + Ip[:, None]*Up[j][None, :])/(2*np.pi)     # int Phi2 dl
        M = -2*(V[:, 0, None, None]*E1m + V[:, 1, None, None]*E2m)@Nm[j]           # (Ns,8,8)
        E[:, 8*j:8*j+8] = M.reshape(8*Ns, 8)
    return E, X, L

def T1_matrix(E, Nr, kappa, a=None, b=None):
    Ns = len(Nr); J = np.zeros_like(E)
    a = 1.0 if a is None else a; b = np.sqrt(kappa) if b is None else b
    for j in range(Ns):
        J[8*j:8*j+8, 8*j:8*j+8] = Jblock(np.array([Nr[j][0], Nr[j][1], 0.]), kappa, 1.0, a, b)
    I = np.eye(len(E)); return 0.5*(I + E) + 0.5*(I - E)@J

def weighted_svals(T, X, L, V, beta, nsmall=8):
    r = np.min([np.linalg.norm(X - v, axis=1) for v in V], axis=0)
    w = np.repeat(r**beta, 8); d = np.repeat(np.sqrt(L), 8)
    A = (d*w)[:, None]*T/(d*w)[None, :]
    s = np.linalg.svd(A, compute_uv=False)
    return s[-nsmall:][::-1], s[0]

if __name__ == '__main__':
    import sys, time
    P0, P1, Nr, V = square_mesh(16); E, X, L = E_matrix(P0, P1, Nr)
    print("Kontrolle ||E^2 - 1|| (N=16, relativ, Spektralnorm):", np.linalg.norm(E@E - np.eye(len(E)), 2))

def trace(X, xs, c):
    z = np.concatenate([X - xs, np.zeros((len(X), 1))], 1)
    K = z/(2*np.pi*np.sum(z**2, 1, keepdims=True))
    return np.concatenate([(MV.vec(k)*MV(c)).c for k in K])

def check_plemelj(N):
    P0, P1, Nr, V = square_mesh(N); E, X, L = E_matrix(P0, P1, Nr)
    rng = np.random.default_rng(0); c = rng.normal(size=8) + 1j*rng.normal(size=8)
    hi = trace(X, np.array([2.1, 0.7]), c); he = trace(X, np.array([0.2, -0.3]), c)
    d = np.repeat(np.sqrt(L), 8); nrm = lambda v: np.linalg.norm(d*v)
    return nrm(E@hi - hi)/nrm(hi), nrm(E@he + he)/nrm(he)
