"""Strahlungsdruck auf eine (geschichtete) Kugel nach Mie (Bohren-Huffman 4.61, 4.62): F_z = sigma_pr 1/2 eps2 |E0|^2,
sigma_pr = sigma_ext - <cos theta> sigma_sca mit
  sigma_ext = (2 pi/k^2) sum (2n+1) Re(a_n + b_n),
  <cos theta> sigma_sca = (4 pi/k^2) [sum n(n+2)/(n+1) Re(a_n a*_{n+1} + b_n b*_{n+1}) + sum (2n+1)/(n(n+1)) Re(a_n b*_n)].
Aufruf: python3 tools/mie_force.py   (Selbsttest: verlustfreie Kugel, sigma_pr = sigma_sca (1 - <cos theta>))
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mie_coated as mc


def sigmas(om, radii, eps, eps2=1.0, nmax=None):
    k = om * np.sqrt(eps2)
    n, a, b = mc.coeffs_layered(om, radii, eps, eps2, nmax)
    ext = 2 * np.pi / k ** 2 * np.sum((2 * n + 1) * np.real(a + b))
    sca = 2 * np.pi / k ** 2 * np.sum((2 * n + 1) * (abs(a) ** 2 + abs(b) ** 2))
    g = 4 * np.pi / k ** 2 * (np.sum(n[:-1] * (n[:-1] + 2) / (n[:-1] + 1) * np.real(a[:-1] * np.conj(a[1:]) + b[:-1] * np.conj(b[1:])))
                             + np.sum((2 * n + 1) / (n * (n + 1)) * np.real(a * np.conj(b))))
    return ext, sca, ext - g


def force_z(om, radii, eps, eps2=1.0, E0sq=1.0):
    return sigmas(om, radii, eps, eps2)[2] * 0.5 * np.real(eps2) * E0sq


if __name__ == "__main__":
    e, s, pr = sigmas(0.8, [1.0], [4.0], 1.0)
    print(f"verlustfrei: sigma_ext {e:.10f}, sigma_sca {s:.10f}, sigma_pr {pr:.10f}")
    print("Goldkugel in Wasser (omega a = 0,5): F_z =", force_z(0.5, [1.0], [-11 + 1.2j], 1.7689))
