"""Anreicherung der Eckelemente mit dem Singulaerprofil r^{a-1} (a = Kantenexponent) und
Richardson-Extrapolation mit bekanntem komplexem Exponenten."""
import numpy as np
from cl3 import MV
from square_geo import build_geo
from exponents import exponents
Vv = np.array([(1, -1), (1, 1), (-1, 1), (-1, -1)], float)
def far_field(ng, kappa, beta_corner=None, x=np.array([5.0, 3.0]), q=0.5, nmid=4, q0=0.5, kappa_fn=None):
    A, X, rr = build_geo(ng, 0.0, kappa, q=q, beta_corner=beta_corner, nmid=nmid, q0=q0, kappa_fn=kappa_fn)
    d = build_geo.last; L = d['L']; bet = d['bet']; T0, T1 = d['T0'], d['T1']
    Ns = len(X); hin = np.zeros((Ns, 8), complex); hin[:, 1] = 1.0
    rs = np.repeat(d['rX']**bet*np.sqrt(L), 8); cs = np.repeat(1/np.sqrt(L), 8)
    g = np.linalg.solve(A, rs*hin.ravel()); xc = (cs*g).reshape(Ns, 8)      # Koeffizienten der Ansatzfunktionen
    # Integral der Ansatzfunktion ueber das Element (Eckelemente: int_0^h r^{-bet} dr)
    I = np.where(bet == 0, L, (np.minimum(T1, 2 - T0) - 0)**(1 - bet)/(1 - bet))
    I = np.where((T0 == 0.0) | (T1 == 2.0), L**(1 - bet)/(1 - bet), L)
    M = Ns//4; n = []
    for k in range(4):
        u = (Vv[(k+1) % 4] - Vv[k])/2; n += [np.array([u[1], -u[0]])]*M
    Fs = np.zeros(8, complex)
    for j in range(Ns):
        z = np.array([x[0]-X[j][0], x[1]-X[j][1], 0.]); K = z/(2*np.pi*(z@z))
        dens = xc[j]*I[j] - hin[j]*L[j]
        Fs += (MV.vec(K)*MV.vec([n[j][0], n[j][1], 0])*MV(dens)).c
    return Fs
if __name__ == '__main__':
    import sys
    for kappa in [-2.8+0.3j, -3+0.3j, -1.5+0.2j, -11+1.2j]:
        a = exponents(kappa)[0]
        base, enr = [], []
        for ng in [8, 16, 24, 32]:
            base.append(far_field(ng, kappa)[1]); enr.append(far_field(ng, kappa, beta_corner=1 - a)[1])
        db = [abs(base[i+1]-base[i]) for i in range(3)]; de = [abs(enr[i+1]-enr[i]) for i in range(3)]
        # Richardson mit bekanntem Exponenten: F(h) = F_inf + C h^a,  h-Verhaeltnis zwischen ng-Stufen 2^-8
        rho = (2.0**-8)**a
        rich = [(base[i+1] - rho*base[i])/(1 - rho) for i in range(3)]
        dr = [abs(rich[i+1]-rich[i]) for i in range(2)]
        print(f"kappa={kappa} a={a:.4f}: ohne {np.round(db,5)} | angereichert {np.round(de,5)} | Richardson {np.round(dr,5)}  (F32 ohne {base[-1]:.5f}, angereichert {enr[-1]:.5f}, Richardson {rich[-1]:.5f})", flush=True)
