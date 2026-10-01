"""Dipol-Polarisierbarkeiten einer kleinen (auch chiralen) Kugel fuer die Kraft in Dipolnaeherung (optical_force.hpp, v0.33):
p/sqrt(eps) = A_e e + i A_c h, m/sqrt(mu) = A_m h - i A_c e mit e = sqrt(eps) E, h = sqrt(mu) H.
Aus den Koeffizienten n = 1 der chiralen Mie-Loesung (mie_chiral_layered.coefficients; gestreut an M + bn N, einfallend
C (M + s N)). Die Abbildung auf Bohren-Huffman ist am achiralen Grenzfall bestimmt: an/C = -b_1, bn/C = -s a_1. Damit
  A_s = 6 pi i [-(an + s bn)/(2C)] / k^3   (gleiche Helizitaet s, Beschriftung wie circular_polarization),
  (A_e - A_m)/2 = 6 pi i (an - s bn)/(2C) / k^3,   A_c = (A_+ - A_-)/2,   (A_e + A_m)/2 = (A_+ + A_-)/2.
Achirale Kugel: A_e = 6 pi i a_1/k^3, A_m = 6 pi i b_1/k^3, A_c = 0.
Aufruf: python3 tools/mie_polarizability.py omega radius eps chi eps2   (gibt Ae Am Ac als Re Im aus)
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mie_chiral_layered as mcl


def polarizability(om, R, eps, chi, eps2):
    k = om * np.sqrt(eps2); c6 = 6j * np.pi / k ** 3
    As, Ax = {}, []
    for s in [+1, -1]:
        C, an, bn = mcl.coefficients(1, om, [R], [(eps, 1.0, chi)], s, eps2)
        As[s] = c6 * (-(an + s * bn) / (2 * C)); Ax.append(c6 * (an - s * bn) / (2 * C))
    Ax = np.mean(Ax); Ap = (As[1] + As[-1]) / 2; Ac = (As[1] - As[-1]) / 2
    return Ap + Ax, Ap - Ax, Ac


if __name__ == "__main__":
    a = [float(x) for x in sys.argv[1:6]] if len(sys.argv) > 5 else [0.5, 0.08, 2.25, 0.2, 1.7689]
    Ae, Am, Ac = polarizability(*a)
    print(Ae.real, Ae.imag, Am.real, Am.imag, Ac.real, Ac.imag)
