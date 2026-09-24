import mpmath as mp, numpy as np
from cone import kappa
mp.mp.dps = 18
def exponent(k0, m, thd):
    th = mp.mpf(thd)*mp.pi/180
    best = None
    for re0 in [-0.4, -0.2]:
        for im0 in [-0.8, -0.3, 0.3]:
            try:
                nu = mp.findroot(lambda n: kappa(n, m, th) - k0, mp.mpc(re0, im0), tol=1e-14, maxsteps=60)
                if -0.5 < mp.re(nu) < 0.5 and abs(kappa(nu, m, th) - k0) < 1e-10:
                    if best is None or mp.re(nu) < mp.re(best): best = nu
            except Exception: pass
    return best
if __name__ == '__main__':
    for thd, k0, name in [(20.0, -11+1.2j, 'Gold/Vakuum'), (20.0, -15+0.5j, 'Silber/Vakuum'), (41.41, (-11+1.2j)/1.77, 'Gold/Wasser'),
                          (41.41, -3+0.3j, '-3+0,3i'), (60.0, -2+0.3j, '-2+0,3i')]:
        nu = exponent(mp.mpc(k0), 0, thd)
        if nu is None: print(thd, name, 'kein Exponent im Streifen gefunden'); continue
        print(f"theta0={thd}: {name}: nu = {mp.nstr(nu, 6)}  -> Feld ~ r^(nu-1); noetig c3 < Re nu, d.h. Gewicht r^(2 beta) mit beta > {mp.nstr(-mp.re(nu), 4)}", flush=True)
