"""Zerfallsraten eines elektrischen Dipols vor einer (geschichteten) Kugel, radial und tangential (Reihenloesung nach Mie;
vgl. Ruppin 1982, Chew 1987; Koeffizienten a_n, b_n nach Bohren-Huffman aus tools/mie_coated.py).

rho = k r0 (Abstand r0 vom Kugelmittelpunkt, Wellenzahl k aussen), h_n = spharische Hankelfunktion erster Art,
psi_n = rho j_n, xi_n = rho h_n, Strich = d/d rho:
  radial:      gamma/gamma0 = 1 + s (3/2) Re sum n(n+1)(2n+1) a_n (h_n/rho)^2
               rad          = (3/2) sum n(n+1)(2n+1) |(j_n + s a_n h_n)/rho|^2
  tangential:  gamma/gamma0 = 1 + s (3/4) Re sum (2n+1) [b_n h_n^2 + a_n (xi_n'/rho)^2]
               rad          = (3/4) sum (2n+1) [|j_n + s b_n h_n|^2 + |(psi_n' + s a_n xi_n')/rho|^2]
Das Vorzeichen s wird nicht vorausgesetzt, sondern ueber die verlustfreie Kugel bestimmt (gesamt = strahlend); siehe __main__.
"""
import os, sys
import numpy as np
from scipy.special import spherical_jn, spherical_yn
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mie_coated as mc

S = -1   # Vorzeichen der Koeffizienten: bestimmt ueber die verlustfreie Kugel (gesamt = strahlend), s = +1 ergaebe negative Raten


def rates(om, radii, eps, eps2, r0, orient='radial', nmax=None, s=None):
    s = S if s is None else s
    k = om * np.sqrt(eps2); rho = k * r0
    R = radii[-1]
    # Reihe des freien Dipols bis etwa k r0, die der Kugel nur so weit, wie sie beitraegt (sonst Ueber-/Unterlauf der a_n, b_n)
    ns = min(int(20 + 3 * abs(k) * R + 15 * R / (r0 - R)), 120)
    nmax = nmax or max(ns, int(30 + 2 * abs(k) * r0))
    nn_s, a_s, b_s = mc.coeffs_layered(om, radii, eps, eps2, ns)
    a = np.zeros(nmax, complex); b = np.zeros(nmax, complex)
    a[:len(a_s)] = np.nan_to_num(a_s); b[:len(b_s)] = np.nan_to_num(b_s)
    nn = np.arange(1, nmax + 1)
    tot = 1.0 + 0j; rad = 0.0
    for i, n in enumerate(nn):
        j = spherical_jn(n, rho); dj = spherical_jn(n, rho, True); y = spherical_yn(n, rho); dy = spherical_yn(n, rho, True)
        h = j + 1j * y; dh = dj + 1j * dy
        dpsi = j + rho * dj; dxi = h + rho * dh
        if orient == 'radial':
            tot += s * 1.5 * n * (n + 1) * (2 * n + 1) * a[i] * (h / rho) ** 2
            rad += 1.5 * n * (n + 1) * (2 * n + 1) * abs((j + s * a[i] * h) / rho) ** 2
        else:
            tot += s * 0.75 * (2 * n + 1) * (b[i] * h ** 2 + a[i] * (dxi / rho) ** 2)
            rad += 0.75 * (2 * n + 1) * (abs(j + s * b[i] * h) ** 2 + abs((dpsi + s * a[i] * dxi) / rho) ** 2)
    return np.real(tot), rad


def rates_mp(om, R, eps1, eps2, r0, orient='radial', nmax=None, dps=60):
    """Homogene Kugel in hoher Genauigkeit (mpmath): dicht an der Kugel ueberlaufen a_n h_n(k r0)^2 in doppelter Genauigkeit
    (a_n winzig, h_n riesig), obwohl gerade diese Ordnungen n ~ R/(r0 - R) das Quenching bestimmen. Koeffizienten nach
    Bohren-Huffman (4.53) mit mu = 1, Vorzeichen s = S."""
    import mpmath as mp
    mp.mp.dps = dps
    k = mp.mpf(om) * mp.sqrt(mp.mpc(eps2)); m = mp.sqrt(mp.mpc(eps1) / mp.mpc(eps2)); x = k * R; rho = k * mp.mpf(r0)
    nmax = nmax or int(40 + 2 * abs(complex(k)) * r0 + 25 * R / (r0 - R))
    jn = lambda n, z: mp.sqrt(mp.pi / (2 * z)) * mp.besselj(n + mp.mpf(1) / 2, z)
    yn = lambda n, z: mp.sqrt(mp.pi / (2 * z)) * mp.bessely(n + mp.mpf(1) / 2, z)
    def fun(n, z, kind):                                  # f_n(z) und d/dz [z f_n(z)]
        f = jn(n, z) if kind == 'j' else jn(n, z) + 1j * yn(n, z)
        fm = jn(n - 1, z) if kind == 'j' else jn(n - 1, z) + 1j * yn(n - 1, z)
        return f, z * fm - n * f                          # (z f_n)' = z f_{n-1} - n f_n
    tot = mp.mpc(1); rad = mp.mpf(0); s = S
    for n in range(1, nmax + 1):
        jx, dpx = fun(n, x, 'j'); hx, dxx = fun(n, x, 'h'); jm, dpm = fun(n, m * x, 'j')
        a = (m ** 2 * jm * dpx - jx * dpm) / (m ** 2 * jm * dxx - hx * dpm)
        b = (jm * dpx - jx * dpm) / (jm * dxx - hx * dpm)
        jr, dpr = fun(n, rho, 'j'); hr, dxr = fun(n, rho, 'h')
        if orient == 'radial':
            tot += s * mp.mpf(3) / 2 * n * (n + 1) * (2 * n + 1) * a * (hr / rho) ** 2
            rad += mp.mpf(3) / 2 * n * (n + 1) * (2 * n + 1) * abs((jr + s * a * hr) / rho) ** 2
        else:
            tot += s * mp.mpf(3) / 4 * (2 * n + 1) * (b * hr ** 2 + a * (dxr / rho) ** 2)
            rad += mp.mpf(3) / 4 * (2 * n + 1) * (abs(jr + s * b * hr) ** 2 + abs((dpr + s * a * dxr) / rho) ** 2)
    return float(mp.re(tot)), float(rad)


if __name__ == "__main__":
    print("Vorzeichen s ueber die verlustfreie Kugel (eps = 4, eps2 = 1, omega a = 0,8, r0 = 1,3): gesamt / strahlend")
    for s in [+1, -1]:
        for o in ['radial', 'tangential']:
            t, r = rates(0.8, [1.0], [4.0], 1.0, 1.3, o, s=s)
            print(f"  s = {s:+d}, {o:10s}: {t:.10f} / {r:.10f}  (Differenz {t - r:+.1e})")
    print("hohe Genauigkeit gegen doppelte (Goldkugel, r0 = 1,5, beide stabil):",
          ["%.8f" % v for v in rates(0.5, [1.0], [-11 + 1.2j], 1.7689, 1.5, 'radial')], ["%.8f" % v for v in rates_mp(0.5, 1.0, -11 + 1.2j, 1.7689, 1.5, 'radial')])
    print("Grenzfall grosser Abstand (r0 = 60): gesamt, strahlend -> 1")
    for o in ['radial', 'tangential']:
        print(f"  {o:10s}:", *[f"{v:.6f}" for v in rates(0.5, [1.0], [-11 + 1.2j], 1.0, 60.0, o)])
