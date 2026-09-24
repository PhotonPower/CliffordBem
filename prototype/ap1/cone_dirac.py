"""Konormalsymbol von T_1 an der Kegelspitze (statisch, Mode m, Mellin-Grad lambda):
Hardy-Raeume R_+ / R_- = Spuren homogener monogener Funktionen F = grad(H) c im Kegelinneren
bzw. -aeusseren am Meridianpunkt x0 (r=1, phi=0). H = r^nu P_nu^{|mu|}(+-cos th) e^{i mu phi},
nu = lambda + 1, c aus dem Eigenraum der Drehung zum Index q = m - mu in {-1,0,1}.
Kritisch: R_+ cap J R_- != {0}, d.h. det[B_+, J B_-] = 0; Laurent-Polynom in s = sqrt(kappa)."""
import numpy as np, mpmath as mp
from cl3 import MV
from cone import P, dP, kappa as kappa_phys
from corner import Jblock
mp.mp.dps = 18
# Eigenbasis der Drehung um e3 (Rotor-Konjugation): Index 0, +1, -1
def blade(b): v = np.zeros(8, complex); v[b] = 1; return v
Q = {0: [blade(0), blade(4), blade(3), blade(7)],                 # 1, e3, e12, e123
     1: [blade(1) - 1j*blade(2), blade(5) - 1j*blade(6)],          # Index +1 unter Drehung um +psi
     -1: [blade(1) + 1j*blade(2), blade(5) + 1j*blade(6)]}
def check_indices():
    """Drehung um psi: R c R~ = e^{i q psi} c ?"""
    psi = 0.37; R = MV([np.cos(psi/2), 0, 0, -np.sin(psi/2), 0, 0, 0, 0])   # R = cos - e12 sin (Drehung um +psi)
    Rt = R.reverse()
    for q, vs in Q.items():
        for v in vs:
            w = (R*MV(v)*Rt).c
            lam = w@np.conj(v)/(v@np.conj(v))
            assert np.allclose(w, lam*v), (q, v)
            print(f"Index {q:+d}: Eigenwert {np.round(lam, 6)}  erwartet e^(i q psi) = {np.round(np.exp(1j*q*psi), 6)}")
def traces(nu, m, th0, side):
    th0 = mp.mpf(th0); s0, c0 = mp.sin(th0), mp.cos(th0)
    rh = np.array([float(s0), 0, float(c0)]); thh = np.array([float(c0), 0, -float(s0)]); phh = np.array([0, 1., 0])
    cols = []
    for q, vs in Q.items():
        mu = m + q; am = abs(mu)          # Kovarianz: R c R~ = e^{i(mu - m) psi} c
        if side == +1: Pv, dth = P(nu, am, c0), -s0*dP(nu, am, c0)
        else:          Pv, dth = P(nu, am, -c0), s0*dP(nu, am, -c0)
        Pv, dth = complex(Pv), complex(dth)
        gradH = complex(nu)*Pv*rh + dth*thh + 1j*mu*Pv/float(s0)*phh
        for v in vs:
            cols.append((MV.vec(gradH)*MV(v)).c)
    A = np.array(cols).T
    U, sv, _ = np.linalg.svd(A); r = int(np.sum(sv > 1e-9*sv[0]))
    return U[:, :r], sv
def J_s(n, s, a=1.0, b=None):
    """J mit explizitem s = sqrt(kappa) (mu1 = mu2), Standardwahl a = 1, b = s."""
    b = s if b is None else b
    Jm = np.zeros((8, 8), complex)
    for bb in range(8):
        H = np.zeros(8, complex); H[bb] = 1
        v = H[[1, 2, 4]]; Ib = H[[6, 5, 3]]*np.array([1, -1, 1])
        vn, bn = n@v, n@Ib
        v2 = s*(v - vn*n) + vn*n/s; b2 = Ib
        out = np.zeros(8, complex); out[0] = a*H[0]; out[7] = b*H[7]
        out[[1, 2, 4]] = v2; out[[6, 5, 3]] = b2*np.array([1, -1, 1]); Jm[:, bb] = out
    return Jm
def crit_kappas(nu, m, th0, choice='standard'):
    Bp, _ = traces(nu, m, th0, +1); Bm, _ = traces(nu, m, th0, -1)
    n0 = np.array([np.cos(th0), 0, -np.sin(th0)])
    N = 16; ss = np.exp(2j*np.pi*np.arange(N)/N)
    dets = []
    for s in ss:
        J = J_s(n0, s)
        dets.append(np.linalg.det(np.hstack([Bp, J@Bm])))
    c = np.fft.fft(np.array(dets))/N           # Koeffizienten von s^k, k = 0..N-1 (negative Potenzen gefaltet)
    coef = np.array([c[k % N] for k in range(-4, 5)])   # s^-4 .. s^4
    roots = np.roots(coef[::-1])
    roots = roots[(roots.real > 1e-12) | ((abs(roots.real) <= 1e-12) & (roots.imag > 0))]   # Hauptzweig s = sqrt(kappa)
    return roots**2, np.array(dets), Bp.shape[1], Bm.shape[1]
if __name__ == '__main__':
    check_indices()
    th0 = 41.41*np.pi/180
    for m in [0, 1, 2]:
        for tau in [0.35, 1.2]:
            nu = mp.mpc(0, tau)                   # L^2-Gerade: nu = i tau (Feldgrad lambda = nu - 1)
            ks, dets, rp, rm = crit_kappas(nu, m, th0)
            kp = complex(kappa_phys(nu, m, mp.mpf(th0)))
            ks = ks[np.isfinite(ks)]
            print(f"m={m}, tau={tau}: dim R+={rp}, R-={rm}; kritische kappa (Dirac): {np.round(np.sort_complex(ks[np.abs(ks)<1e6]),5)}")
            print(f"            physikalisch kappa(nu,m) = {np.round(kp,5)}", flush=True)
