"""Hauptwinkel auf Symbolebene: Grenzwert der Winkel fuer hohe Moden (|xi| -> oo)."""
import numpy as np
from cl3 import MV
def sigma():
    xi, e3 = MV.vec([1,0,0]), MV.vec([0,0,1])
    return 1j*(xi*e3).left_matrix()
def Jmat(re_, rm, a, b):
    d = {0:a, 1:np.sqrt(re_), 2:np.sqrt(re_), 4:1/np.sqrt(re_), 3:1/np.sqrt(rm), 5:np.sqrt(rm), 6:np.sqrt(rm), 7:b}
    return np.diag([d[i] for i in range(8)]).astype(complex)
def rng(P):
    U, s, _ = np.linalg.svd(P); return U[:, s > 1e-10]
def sym_angle(re_, rm=1.0, a=None, b=None, dual=False, W='J', maxwell=False):
    a = np.sqrt(rm) if a is None else a; b = np.sqrt(re_) if b is None else b
    S = sigma(); pp, pm = (np.eye(8)+S)/2, (np.eye(8)-S)/2; J = Jmat(re_, rm, a, b); Ji = np.linalg.inv(J)
    if not dual: A, B = Ji@rng(pp), rng(pm)
    else: A, B = rng(pp), (Ji if W == 'J' else np.eye(8))@rng(pm)
    if maxwell:   # nur Grade 1,2: Schnitt der Raeume mit {Grad 0,3 = 0}
        Pg = np.diag([0,1,1,1,1,1,1,0]).astype(complex)
        def restrict(Q):
            # Unterraum von span(Q) mit verschwindenden Graden 0,3
            N = Q[[0,7], :]; _, s, Vh = np.linalg.svd(N); null = Vh[len(s[s>1e-10]):].conj().T
            return Q@null
        A, B = restrict(A), restrict(B)
    QA, _ = np.linalg.qr(A); QB, _ = np.linalg.qr(B)
    c = np.linalg.svd(QA.conj().T@QB, compute_uv=False)[0]
    return np.sqrt(max(0, 1-min(c,1)**2))
if __name__ == '__main__':
    for r in [1.0, 2.25, 12.0, -11+1.2j, -2+0.3j, -1+0.1j]:
        print(f"r_eps={r}: voll {sym_angle(r):.4f}  Maxwell {sym_angle(r, maxwell=True):.4f}  dual(W=J^-1) {sym_angle(r, dual=True):.4f}  dual(W=1) {sym_angle(r, dual=True, W='one'):.4f}")
