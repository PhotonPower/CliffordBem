"""Mie-Referenz fuer die beschichtete Kugel (Aden/Kerker 1951; Bohren/Huffman, Abschn. 8.1), nichtmagnetisch.

Kern (Radius a, eps_core), Schale bis Radius b (eps_shell), Aussenmedium eps2. Laengen in Einheiten, omega = k0.
Beliebig viele Schalen ueber die Transfermatrix der Riccati-Bessel-Funktionen (radii aufsteigend, eps je Schicht).

Aufruf:
  python3 tools/mie_coated.py --omega 0.5 --radii 1,1.05 --eps " -11,1.2;2.25" [--eps2 1.0]
  -> Q_ext (bezogen auf pi b^2), Vorwaertsamplitude S(0) = sum (2n+1)/2 (a_n + b_n), |S|, arg S.
"""
import argparse
import numpy as np
from scipy.special import spherical_jn, spherical_yn


def _psi(n, z):  return z * spherical_jn(n, z)
def _dpsi(n, z): return spherical_jn(n, z) + z * spherical_jn(n, z, derivative=True)
def _chi(n, z):  return z * spherical_yn(n, z)
def _dchi(n, z): return spherical_yn(n, z) + z * spherical_yn(n, z, derivative=True)


def coeffs_layered(om, radii, eps, eps2=1.0, nmax=None):
    """a_n, b_n fuer eine geschichtete Kugel. radii: aufsteigend (Kern zuerst), eps: Permittivitaet je Schicht.
    Innen: Feld = c_n psi(k_l r) + d_n chi(k_l r); Stetigkeit von E_t, H_t an jeder Grenzflaeche."""
    radii = np.asarray(radii, float); eps = np.asarray(eps, complex)
    k2 = om * np.sqrt(complex(eps2)); x_out = k2 * radii[-1]
    if nmax is None: nmax = int(abs(x_out) + 4 * abs(x_out) ** (1 / 3) + 12)
    n = np.arange(1, nmax + 1)
    ks = om * np.sqrt(eps.astype(complex))
    m = np.sqrt(eps.astype(complex) / eps2)
    res = {}
    for mode in ("a", "b"):
        # Koeffizienten im Kern: nur psi
        C = np.ones_like(n, dtype=complex); D = np.zeros_like(n, dtype=complex)
        for l in range(len(radii) - 1):
            r = radii[l]; z1 = ks[l] * r; z2 = ks[l + 1] * r
            f1 = C * _psi(n, z1) + D * _chi(n, z1); g1 = C * _dpsi(n, z1) + D * _dchi(n, z1)
            # Stetigkeit bis auf einen gemeinsamen Faktor: TE (b_n): (f, m f'), TM (a_n): (m f, f')
            if mode == "b":  u, v = f1, m[l] * g1; wu, wv = 1.0, m[l + 1]
            else:            u, v = m[l] * f1, g1; wu, wv = m[l + 1], 1.0
            # loese  wu (C' psi + D' chi)(z2) = u,  wv (C' psi' + D' chi')(z2) = v
            P, Q, dP, dQ = _psi(n, z2), _chi(n, z2), _dpsi(n, z2), _dchi(n, z2)
            det = P * dQ - Q * dP
            Cn = (u / wu * dQ - Q * v / wv) / det
            Dn = (P * v / wv - dP * u / wu) / det
            C, D = Cn, Dn
        r = radii[-1]; z = ks[-1] * r; x = k2 * r
        f = C * _psi(n, z) + D * _chi(n, z); g = C * _dpsi(n, z) + D * _dchi(n, z)
        xi = x * (spherical_jn(n, x) + 1j * spherical_yn(n, x))
        dxi = (spherical_jn(n, x) + 1j * spherical_yn(n, x)) + x * (spherical_jn(n, x, True) + 1j * spherical_yn(n, x, True))
        ps, dps = _psi(n, x), _dpsi(n, x); ml = m[-1]
        if mode == "a":   # TM
            num = f * ml * dps - g * ps
            den = f * ml * dxi - g * xi
        else:             # TE
            num = f * dps - ml * g * ps
            den = f * dxi - ml * g * xi
        res[mode] = num / den
    return n, res["a"], res["b"]


def qext(om, radii, eps, eps2=1.0):
    n, a, b = coeffs_layered(om, radii, eps, eps2); x = om * np.sqrt(eps2) * radii[-1]
    return float(2 / x ** 2 * np.sum((2 * n + 1) * np.real(a + b)))


def forward_amplitude(om, radii, eps, eps2=1.0):
    n, a, b = coeffs_layered(om, radii, eps, eps2)
    return complex(np.sum((2 * n + 1) / 2 * (a + b)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--omega", type=float, required=True)
    ap.add_argument("--radii", required=True, help="aufsteigend, z. B. 1,1.05")
    ap.add_argument("--eps", required=True, help="je Schicht 're,im', getrennt durch ';'")
    ap.add_argument("--eps2", type=float, default=1.0)
    A = ap.parse_args()
    radii = [float(v) for v in A.radii.split(",")]
    eps = []
    for s in A.eps.split(";"):
        p = [float(v) for v in s.split(",")]; eps.append(complex(p[0], p[1] if len(p) > 1 else 0.0))
    S = forward_amplitude(A.omega, radii, eps, A.eps2)
    print(f"Q_ext = {qext(A.omega, radii, eps, A.eps2):.8f}  S(0) = {S.real:.8f}{S.imag:+.8f}i  |S| = {abs(S):.8f}  arg S = {np.angle(S):.8f}")
