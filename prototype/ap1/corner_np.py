"""Gegenprobe: skalarer quasistatischer Transmissionsoperator (Neumann-Poincare K*) auf L^2
an derselben Ecke. Kritisch: -(kappa+1)/(2(kappa-1)) in spec(Mellin-Symbol von K*),
d.h. kappa = (2 lam - 1)/(2 lam + 1)."""
import numpy as np
from corner import geometry, _x, _dx
def np_symbol(alpha, tau):
    d, n = geometry(alpha); u = np.exp(_x); w = np.exp(0.5*_x + 1j*tau*_x)*_dx
    m = np.zeros((2, 2), complex)
    for i, j in [(0, 1), (1, 0)]:
        z = d[i][None, :] - u[:, None]*d[j][None, :]
        m[i, j] = np.sum((z@n[i])/(2*np.pi*np.sum(z**2, 1))*w)
    return m
def crit_np(alpha, tau):
    lam = np.linalg.eigvals(np_symbol(alpha, tau))
    return (2*lam - 1)/(2*lam + 1)
if __name__ == '__main__':
    from corner_crit import crit_kappas
    for alpha in [np.pi/2, np.pi/3, 2*np.pi/3]:
        for tau in [0.0, 0.25, 1.0]:
            kd = crit_kappas(alpha, tau); kd = np.unique(np.round(kd[np.abs(kd-1) > 1e-6], 6))
            print(f"alpha={alpha/np.pi:.3f}pi tau={tau}: skalar {np.round(np.sort_complex(crit_np(alpha,tau)),4)}   Dirac {np.round(np.sort_complex(kd),4)}")
