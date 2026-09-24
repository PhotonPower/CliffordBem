"""Rotationskoerper-Zerlegung (Body of Revolution) fuer die Subraumwinkel-Tests.
Alle Operatoren kommutieren mit Drehungen um e3 (U_phi h)(x) = R h(R^-1 x) R~.
Fuer Azimutalmodus m ist ein Feld durch seine Werte auf einem Meridian bestimmt;
die Moden sind in L^2(Gamma) orthogonal, also  sin theta_min = min_m sin theta_min(m)."""
import numpy as np
from fields import dirac_kernel, maxwell_dipole, gp, PSEUDO, blade
from cl3 import SIGN

def rotz(phi):
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.]])

VEC, BIV = [1, 2, 4], [6, 5, 3]            # e1,e2,e3 ; I e1 = e23, I e2 = -e13, I e3 = e12
BSIGN = np.array([1, -1, 1])

def rotate_mv(F, Rm):
    """Drehung eines Multivektorfeldes (...,8) mit Rotationsmatrix Rm (Rotor-Wirkung)."""
    out = F.copy()
    out[..., VEC] = F[..., VEC] @ Rm.T
    out[..., BIV] = ((F[..., BIV] * BSIGN) @ Rm.T) * BSIGN
    return out

class BoR:
    """Meridiankurve: Punkte (rho,z), Normalen (n_rho,n_z), Linienelemente; Quellen je Meridian."""
    def __init__(self, P, Nrm, dl, src_in, src_out, nrot=96):
        self.Xm = np.stack([P[:, 0], 0*P[:, 0], P[:, 1]], -1)
        self.N = np.stack([Nrm[:, 0], 0*Nrm[:, 0], Nrm[:, 1]], -1)
        self.W = 2*np.pi*P[:, 0]*dl                       # Gewicht inkl. Azimut-Integral
        self.src_in = np.stack([src_in[:, 0], 0*src_in[:, 0], src_in[:, 1]], -1)
        self.src_out = np.stack([src_out[:, 0], 0*src_out[:, 0], src_out[:, 1]], -1)
        self.phis = 2*np.pi*np.arange(nrot)/nrot
        self.Rs = [rotz(p) for p in self.phis]
        self.X = self.Xm
    def _raw(self, x, kind, k, eps, S):
        """Basisfelder (Q, ncol, 8) an Punkten x fuer Quellen S (auf Meridian u=0)."""
        if kind == 'dirac':
            K = dirac_kernel(x[:, None, :], S[None, :, :], k)
            out = np.zeros(K.shape[:2] + (8, 8), dtype=complex)
            for b in range(8):
                for a_ in range(8):
                    out[:, :, b, a_ ^ b] += SIGN[a_, b]*K[:, :, a_]
            return out.reshape(len(x), -1, 8)
        out = []
        for p in np.eye(3):
            F = maxwell_dipole(x[:, None, :], S[None, :, :], p, k, eps)
            out += [F, gp(PSEUDO, F)]
        return np.stack(out, 2).reshape(len(x), -1, 8)
    def traces_m(self, m, kind, k, eps, where):
        """Spuren der Modus-m-Basis  B_m = sum_j e^{-i m phi_j} U_{phi_j} F  auf dem Meridian: (Q,8,ncol)."""
        S = self.src_out if where == 'int' else self.src_in
        acc = 0
        for phi, Rm in zip(self.phis, self.Rs):
            xr = self.Xm @ Rm           # R^-1 x  (Zeilenvektoren: x @ R = R^T x)
            F = self._raw(xr, kind, k, eps, S)
            acc = acc + np.exp(-1j*m*phi)*rotate_mv(F, Rm)
        return acc.transpose(0, 2, 1)
    def traces_all(self, ms, kind, k, eps, where):
        """wie traces_m, aber fuer mehrere Moden gleichzeitig (Feldauswertung nur einmal je Drehung)."""
        S = self.src_out if where == 'int' else self.src_in
        acc = {m: 0 for m in ms}
        for phi, Rm in zip(self.phis, self.Rs):
            F = rotate_mv(self._raw(self.Xm @ Rm, kind, k, eps, S), Rm)
            for m in ms: acc[m] = acc[m] + np.exp(-1j*m*phi)*F
        return {m: acc[m].transpose(0, 2, 1) for m in ms}
    J = __import__('mfs').Sphere.J
    orth = __import__('mfs').Sphere.orth

def angle_m(g, m, kind, om, e1, m1, e2, m2, a=None, b=None, dual=False, W='J'):
    k1, k2 = om*np.sqrt(e1)*np.sqrt(m1), om*np.sqrt(e2)*np.sqrt(m2)
    re_, rm = e1/e2, m1/m2
    a = np.sqrt(rm) if a is None else a; b = np.sqrt(re_) if b is None else b
    if not dual:
        A = g.J(g.traces_m(m, kind, k1, e1, 'int'), re_, rm, a, b, inverse=True)
        B = g.traces_m(m, kind, k2, e2, 'ext')
    else:
        A = g.traces_m(m, kind, k2, e2, 'int')
        B = g.traces_m(m, kind, k1, e1, 'ext')
        if W == 'J': B = g.J(B, re_, rm, a, b, inverse=True)
    QA, QB = g.orth(A), g.orth(B)
    c = np.linalg.svd(QA.conj().T @ QB, compute_uv=False)[0]
    return np.sqrt(max(0.0, 1 - min(c, 1.0)**2))

def angles(g, ms, kind, om, e1, m1, e2, m2, a=None, b=None, dual=False, W='J'):
    """sin theta_min(m) fuer alle m in ms."""
    k1, k2 = om*np.sqrt(e1)*np.sqrt(m1), om*np.sqrt(e2)*np.sqrt(m2)
    re_, rm = e1/e2, m1/m2
    a = np.sqrt(rm) if a is None else a; b = np.sqrt(re_) if b is None else b
    if not dual:
        TA = g.traces_all(ms, kind, k1, e1, 'int'); TB = g.traces_all(ms, kind, k2, e2, 'ext')
        TA = {m: g.J(TA[m], re_, rm, a, b, inverse=True) for m in ms}
    else:
        TA = g.traces_all(ms, kind, k2, e2, 'int'); TB = g.traces_all(ms, kind, k1, e1, 'ext')
        if W == 'J': TB = {m: g.J(TB[m], re_, rm, a, b, inverse=True) for m in ms}
    out = {}
    for m in ms:
        QA, QB = g.orth(TA[m]), g.orth(TB[m])
        c = np.linalg.svd(QA.conj().T @ QB, compute_uv=False)[0]
        out[m] = np.sqrt(max(0.0, 1 - min(c, 1.0)**2))
    return out

def angle(g, M, *args, **kw):
    vals = {m: angle_m(g, m, *args, **kw) for m in range(-M, M+1)}
    mstar = min(vals, key=vals.get)
    return vals[mstar], mstar, vals

# --- Geometrien ---
def sphere_bor(nq=48, ns=24, rin=0.45, rout=2.2, nrot=96):
    t, w = np.polynomial.legendre.leggauss(nq); th = np.arccos(-t)     # theta in (0,pi)
    P = np.stack([np.sin(th), np.cos(th)], -1); dl = w/np.sin(th)       # d theta = dt/sin(theta)
    ts = np.pi*(np.arange(ns)+0.5)/ns
    Sin = rin*np.stack([np.sin(ts), np.cos(ts)], -1); Sout = rout*np.stack([np.sin(ts), np.cos(ts)], -1)
    return BoR(P, P.copy(), dl, Sin, Sout, nrot)

def torus_bor(R=1.0, r=0.4, nq=96, ns=24, fin=0.5, fout=1.75, nrot=256):
    v = 2*np.pi*np.arange(nq)/nq; vs = 2*np.pi*(np.arange(ns)+0.5)/ns
    P = np.stack([R + r*np.cos(v), r*np.sin(v)], -1); Nr = np.stack([np.cos(v), np.sin(v)], -1)
    dl = r*2*np.pi/nq*np.ones(nq)
    Sin = np.stack([R + fin*r*np.cos(vs), fin*r*np.sin(vs)], -1)
    Sout = np.stack([R + fout*r*np.cos(vs), fout*r*np.sin(vs)], -1)
    return BoR(P, Nr, dl, Sin, Sout, nrot)

if __name__ == '__main__':
    import time
    g = sphere_bor(); t0 = time.time()
    for kind, om, e1, ref in [('maxwell', 1.0, 2.25, 0.63584), ('dirac', 1.0, 2.25, 0.43293), ('dirac', 2.0, -11+1.2j, 0.19281)]:
        v, ms, _ = angle(g, 8, kind, om, e1, 1, 1.0, 1)
        print(f"Kugel-BoR {kind:8s} om={om} eps1={e1}: {v:.5f} (m*={ms})   3D-Referenz {ref}", flush=True)
    print('%.0fs' % (time.time()-t0))
