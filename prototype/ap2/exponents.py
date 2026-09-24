"""Singulaerexponenten an der Ecke: Nullstellen a von S(a)^2 - q^2 im Streifen 0 < Re a < 1,
q = (1+kappa)/(1-kappa). Die Randspur verhaelt sich wie r^{a-1}; sie liegt in L^2(r^{2 beta})
genau fuer Re a > 1/2 - beta. Behauptung: c* (Tabelle Teil IV) = min Re a_j."""
import numpy as np
from scipy.optimize import fsolve
from corner_weights import S, c_star
def exponents(kappa, alpha=np.pi/2):
    q = (1+kappa)/(1-kappa); f = lambda a: S(alpha, a)**2 - q**2
    roots = []
    for x0 in np.linspace(0.02, 0.98, 25):
        for y0 in np.linspace(-3, 3, 25):
            g = lambda v: [f(v[0]+1j*v[1]).real, f(v[0]+1j*v[1]).imag]
            v, info, ier, _ = fsolve(g, [x0, y0], full_output=True)
            a = v[0] + 1j*v[1]
            if ier == 1 and 1e-6 < a.real < 1-1e-6 and abs(f(a)) < 1e-10 and not any(abs(a-r) < 1e-7 for r in roots):
                roots.append(a)
    return sorted(roots, key=lambda a: a.real)
if __name__ == '__main__':
    for k in [-5+0.5j, -3+0.3j, -2+0.3j, -1.5+0.2j, -0.5+0.1j, -6.215+0.678j, -11+1.2j, -2.0+0j, -4.0+0j]:
        ex = exponents(k); m = min(a.real for a in ex) if ex else np.nan
        print(f"kappa={k}: Exponenten {np.round(ex,4)}  min Re a = {m:.3f}   c* = {c_star(k):.3f}")
