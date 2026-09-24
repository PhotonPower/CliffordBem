"""Kegelspitze (Schritt 1.5, Ecken): homogene Loesungen u = r^nu P_nu^m(+-cos theta) e^{i m phi}
des quasistatischen Transmissionsproblems an einem Kreiskegel mit halbem Oeffnungswinkel theta0
(Medium 1 innen). Kritischer Kontrast:
   kappa(nu, m) = - P'(-c) P(c) / (P'(c) P(-c)),   c = cos(theta0),  P = P_nu^m.
Gerade Re nu = 0: L^2-Randspuren (E ~ r^{nu-1}, dS ~ r dr); Re nu = -1/2: Energieraum."""
import mpmath as mp
import numpy as np
mp.mp.dps = 18
def _F(nu, m, z, k=0):
    """d^k/dz^k der Funktion (a)_m (b)_m/m! 2F1(a+m, b+m; 1+m; z), a=-nu, b=nu+1."""
    a, b = -nu, nu + 1
    pref = mp.rf(a, m + k)*mp.rf(b, m + k)/(mp.factorial(m)*mp.rf(1 + m, k))
    return pref*mp.hyp2f1(a + m + k, b + m + k, 1 + m + k, z)
def P(nu, m, x):
    """Ferrers-Funktion P_nu^m(x) = (-1)^m (1-x^2)^{m/2} d^m/dx^m P_nu(x),  P_nu = 2F1(-nu, nu+1; 1; (1-x)/2)."""
    z = (1 - x)/2
    return (-1)**m*(1 - x**2)**(mp.mpf(m)/2)*(-mp.mpf(1)/2)**m*_F(nu, m, z)
def dP(nu, m, x):
    z = (1 - x)/2; w = (1 - x**2)
    F0 = _F(nu, m, z); F1 = _F(nu, m, z, 1)
    base = (-1)**m*(-mp.mpf(1)/2)**m
    return base*(-m*x*w**(mp.mpf(m)/2 - 1)*F0 + w**(mp.mpf(m)/2)*F1*(-mp.mpf(1)/2))
def kappa(nu, m, th0):
    c = mp.cos(th0)
    return -dP(nu, m, -c)*P(nu, m, c)/(dP(nu, m, c)*P(nu, m, -c))
if __name__ == '__main__':
    # Kontrolle: ebener Rand theta0 = pi/2 -> kappa = -1 fuer alle nu, m
    print("Kontrolle eben:", [mp.nstr(kappa(mp.mpc(-0.5, t), m, mp.pi/2), 8) for t, m in [(0.3, 0), (1.2, 1), (0.1, 2)]])
    # Symmetrie nu <-> -1-nu
    print("Symmetrie:", mp.nstr(kappa(mp.mpc(0.2, 0.7), 1, 0.7), 10), mp.nstr(kappa(mp.mpc(-1.2, -0.7), 1, 0.7), 10))
    import json
    out = {}
    for thd in [20.0, 41.41, 60.0, 120.0]:
        th = mp.mpf(thd)*mp.pi/180
        for m in [0, 1, 2, 3]:
            for line in [-0.5, 0.0]:
                taus = np.concatenate([np.linspace(1e-4, 2, 60), np.linspace(2.1, 8, 30)])
                ks = []
                for t in taus:
                    try: k = kappa(mp.mpc(line, t), m, th); ks.append([float(mp.re(k)), float(mp.im(k))])
                    except ZeroDivisionError: pass
                out[f"{thd}|{m}|{line}"] = ks
        e0 = [out[f"{thd}|{m}|-0.5"][0][0] for m in range(4)]
        l0 = [out[f"{thd}|{m}|0.0"][0][0] for m in range(4)]
        print(f"theta0={thd:6.2f}: Energie-Gerade tau->0: kappa(m=0..3) =", np.round(e0, 4), " | L2-Gerade tau->0:", np.round(l0, 4), flush=True)
    json.dump(out, open('cone_curves.json', 'w'))
