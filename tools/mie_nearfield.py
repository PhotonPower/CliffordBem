"""Nahfeld einer (geschichteten) Kugel nach Mie ausserhalb der aeusseren Flaeche (Bohren-Huffman, Kap. 4; e^{-i omega t}).

Einfall entlang z, E_0 = x (Amplitude 1) im Medium eps2 (mu = 1), omega = k0. Aussen gilt fuer jede geschichtete Kugel dieselbe
Entwicklung mit ihren Koeffizienten a_n, b_n (tools/mie_coated.py):
  E_inc = sum E_n (M_o1n^(1) - i N_e1n^(1)),          H_inc = -(k/omega) sum E_n (M_e1n^(1) + i N_o1n^(1)),
  E_s   = sum E_n (i a_n N_e1n^(3) - b_n M_o1n^(3)),   H_s   =  (k/omega) sum E_n (i b_n N_o1n^(3) + a_n M_e1n^(3)),
  E_n = i^n (2n+1) / (n (n+1)); (3) = spharische Hankelfunktion h_n^(1). H in der Normierung des Kerns (F = sqrt(eps) E + I H,
  ebene Welle H = sqrt(eps) d x E). y-Polarisation durch Drehung um 90 Grad um z; zirkular p = x + i s y (wie
  circular_polarization im Kern).
Aufruf: python3 tools/mie_nearfield.py   (Selbsttest: das einfallende Feld aus der Entwicklung ist die ebene Welle)
"""
import os, sys
import numpy as np
from scipy.special import spherical_jn, spherical_yn
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mie_coated as mc


def _pi_tau(nmax, mu):
    pi_ = np.zeros(nmax + 1); tau = np.zeros(nmax + 1); pi_[1] = 1.0
    for n in range(2, nmax + 1): pi_[n] = (2 * n - 1) / (n - 1) * mu * pi_[n - 1] - n / (n - 1) * pi_[n - 2]
    for n in range(1, nmax + 1): tau[n] = n * mu * pi_[n] - (n + 1) * pi_[n - 1]
    return pi_, tau


def _vsh(n, z, dz, rho, th, ph, pin, taun):
    """M_o1n, M_e1n, N_o1n, N_e1n in (r, theta, phi) fuer Radialfunktion z(rho), dz = d(rho z)/d rho / rho."""
    st, cp, sp = np.sin(th), np.cos(ph), np.sin(ph)
    Mo = np.array([0, cp * pin * z, -sp * taun * z]); Me = np.array([0, -sp * pin * z, -cp * taun * z])
    No = np.array([sp * n * (n + 1) * st * pin * z / rho, sp * taun * dz, cp * pin * dz])
    Ne = np.array([cp * n * (n + 1) * st * pin * z / rho, cp * taun * dz, -sp * pin * dz])
    return Mo, Me, No, Ne


def fields_xpol(om, radii, eps, eps2, x, total=True, nmax=None):
    """E, H (kartesisch, komplex) am Punkt x ausserhalb der Kugel fuer x-polarisierten Einfall entlang z."""
    k = om * np.sqrt(eps2)
    nn, a, b = mc.coeffs_layered(om, radii, eps, eps2, nmax)
    r = np.linalg.norm(x); th = np.arccos(np.clip(x[2] / r, -1, 1)); ph = np.arctan2(x[1], x[0]); rho = k * r
    pi_, tau = _pi_tau(len(nn), np.cos(th))
    E = np.zeros(3, complex); H = np.zeros(3, complex)
    for i, n in enumerate(nn):
        En = 1j ** n * (2 * n + 1) / (n * (n + 1))
        j = spherical_jn(n, rho); dj = spherical_jn(n, rho, True); y = spherical_yn(n, rho); dy = spherical_yn(n, rho, True)
        h = j + 1j * y; dh = dj + 1j * dy
        Mo, Me, No, Ne = _vsh(n, h, (h + rho * dh) / rho, rho, th, ph, pi_[n], tau[n])
        E += En * (1j * a[i] * Ne - b[i] * Mo); H += (k / om) * En * (1j * b[i] * No + a[i] * Me)
        if total:
            Mo1, Me1, No1, Ne1 = _vsh(n, j, (j + rho * dj) / rho, rho, th, ph, pi_[n], tau[n])
            E += En * (Mo1 - 1j * Ne1); H += -(k / om) * En * (Me1 + 1j * No1)
    st, ct, sp, cp = np.sin(th), np.cos(th), np.sin(ph), np.cos(ph)
    R = np.array([[st * cp, ct * cp, -sp], [st * sp, ct * sp, cp], [ct, -st, 0]])   # (r, theta, phi) -> (x, y, z)
    return R @ E, R @ H


def fields(om, radii, eps, eps2, x, s=0, total=True):
    """s = 0: x-Polarisation; s = +-1: zirkular p = x + i s y."""
    Ex, Hx = fields_xpol(om, radii, eps, eps2, x, total)
    if s == 0: return Ex, Hx
    Rz = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])                   # Drehung um +90 Grad: x -> y
    Ey, Hy = fields_xpol(om, radii, eps, eps2, Rz.T @ x, total); Ey, Hy = Rz @ Ey, Rz @ Hy
    return Ex + 1j * s * Ey, Hx + 1j * s * Hy


if __name__ == "__main__":
    om, eps2 = 0.5, 1.7689; k = om * np.sqrt(eps2)
    x = np.array([0.7, -0.4, 1.3])
    # Selbsttest: Kugel aus dem Aussenmedium -> a_n = b_n = 0, Gesamtfeld = ebene Welle
    E, H = fields(om, [1.0], [eps2], eps2, x, s=+1)
    Ew = np.array([1, 1j, 0]) * np.exp(1j * k * x[2]); Hw = np.sqrt(eps2) * np.cross([0, 0, 1], [1, 1j, 0]) * np.exp(1j * k * x[2])
    print("ebene Welle aus der Entwicklung: |E - E_w| = %.1e, |H - H_w| = %.1e" % (np.linalg.norm(E - Ew), np.linalg.norm(H - Hw)))
    E, H = fields(1.0, [1.0], [-2 + 0.3j], 1.0, np.array([0.0, 0.0, 1.5]), s=+1)
    print("Goldaehnliche Kugel, Punkt (0, 0, 1.5): |E|^2 =", np.sum(abs(E) ** 2) / 2)
