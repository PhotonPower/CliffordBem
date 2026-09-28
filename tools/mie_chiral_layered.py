"""Mie-Loesung fuer geschichtete Kugeln mit chiralen (Pasteur-)Schichten in Vakuum, zirkular polarisierte Anregung.

Verallgemeinert tools/mie_chiral.py (homogene chirale Kugel) auf beliebig viele Schichten.
Konstitutiv wie im Kern: D = eps E + i chi H, B = mu H - i chi E (e^{-i omega t}).
Basis je Gebiet und Ordnung n (m = s), Radialfunktion z = j (regulaer) bzw. y (irregulaer):
  achiral:  M(z, k), N(z, k),                    H = -(i/eta) (a N + b M),     eta = sqrt(mu/eps)
  chiral:   W_A = M + N (k_A = omega(sqrt(eps mu) + chi)), H = -(i/eta) W_A,
            W_B = M - N (k_B = omega(sqrt(eps mu) - chi)), H = +(i/eta) W_B.
Tangentiale Komponenten auf der Kugel (entlang X_nm und r x X_nm): M -> (z, 0), N -> (0, psi'(u)/u), psi = u z.
Kern: nur regulaer; Schalen: regulaer und irregulaer; aussen: einfallend (j) und gestreut (h).
Je n ein lineares System aus der Stetigkeit von E_tan und H_tan an allen Grenzflaechen.

Aufruf: python3 tools/mie_chiral_layered.py   (Selbsttests gegen mie_coated.py und mie_chiral.py)
"""
import numpy as np
from scipy.special import spherical_jn, spherical_yn


def _rad(n, u, kind):
    """z(u), psi'(u)/u fuer z = j, y oder h (u komplex)."""
    u = complex(u)
    if kind == 'j':
        z = spherical_jn(n, u); zp = spherical_jn(n, u, derivative=True)
    elif kind == 'y':
        z = spherical_yn(n, u); zp = spherical_yn(n, u, derivative=True)
    else:
        z = spherical_jn(n, u) + 1j * spherical_yn(n, u); zp = spherical_jn(n, u, True) + 1j * spherical_yn(n, u, True)
    return z, (z + u * zp) / u


def _basis(n, om, med, r, kind):
    """Tangentiale (E_X, E_rX, H_X, H_rX) der beiden Basisfelder eines Gebiets bei Radius r."""
    eps, mu, chi = med
    eta = np.sqrt(mu / eps) if chi == 0 or True else None
    if chi == 0:
        k = om * np.sqrt(eps * mu)
        z, q = _rad(n, k * r, kind)
        M = np.array([z, 0, 0, -1j / eta * q], complex)          # E = M, H = -(i/eta) N
        N = np.array([0, q, -1j / eta * z, 0], complex)          # E = N, H = -(i/eta) M
        return [M, N]
    kA, kB = om * (np.sqrt(eps * mu) + chi), om * (np.sqrt(eps * mu) - chi)
    hA, hB = -1j / eta, 1j / eta
    zA, qA = _rad(n, kA * r, kind); zB, qB = _rad(n, kB * r, kind)
    WA = np.array([zA, qA, hA * zA, hA * qA], complex)
    WB = np.array([zB, -qB, hB * zB, -hB * qB], complex)
    return [WA, WB]


def coefficients(n, om, radii, media, s):
    """radii aufsteigend, media je Schicht (eps, mu, chi) von innen nach aussen; aussen Vakuum.
    Unbekannte: Kern (2), je Schale (4), aussen (a_n, b_n) gestreut."""
    L = len(radii)
    nunk = 2 + 4 * (L - 1) + 2
    A = np.zeros((4 * L, nunk), complex); rhs = np.zeros(4 * L, complex)
    C = 1j ** n * np.sqrt(4 * np.pi * (2 * n + 1))
    def cols(layer):                          # Spaltenindex der Unbekannten von Gebiet layer (0 = Kern, L = aussen)
        if layer == 0: return 0
        if layer == L: return 2 + 4 * (L - 1)
        return 2 + 4 * (layer - 1)
    for i, r in enumerate(radii):             # Grenzflaeche i zwischen Gebiet i (innen) und i+1 (aussen)
        rows = slice(4 * i, 4 * i + 4)
        # innen
        kinds_in = ['j'] if i == 0 else ['j', 'y']
        c = cols(i)
        for kd in kinds_in:
            for f in _basis(n, om, media[i], r, kd):
                A[rows, c] -= f; c += 1
        # aussen
        if i + 1 == L:
            vac = (1.0, 1.0, 0.0)
            Mh, Nh = _basis(n, om, vac, r, 'h'); c = cols(L)
            A[rows, c] += Mh; A[rows, c + 1] += Nh
            Mj, Nj = _basis(n, om, vac, r, 'j')
            rhs[rows] -= C * (Mj + s * Nj)                        # E_inc = sum C (M + s N)
        else:
            c = cols(i + 1)
            for kd in ['j', 'y']:
                for f in _basis(n, om, media[i + 1], r, kd):
                    A[rows, c] += f; c += 1
    x = np.linalg.solve(A, rhs)
    return C, x[cols(L)], x[cols(L) + 1]


def cross_sections(om, radii, media, s=+1, nmax=None):
    """sigma_ext, sigma_sca (Vakuum aussen, |E_inc|^2 = 2 wie in mie_chiral.py)."""
    R = radii[-1]; nmax = nmax or int(om * R + 4 * (om * R) ** (1 / 3) + 12)
    ext = sca = 0.0
    for n in range(1, nmax + 1):
        C, an, bn = coefficients(n, om, radii, media, s)
        ext += -np.real(np.conj(C) * an + np.conj(s * C) * bn); sca += abs(an) ** 2 + abs(bn) ** 2
    return ext / (2 * om ** 2), sca / (2 * om ** 2)


if __name__ == '__main__':
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import mie_coated as mc, mie_chiral as mch
    om = 0.5; gold = (-11 + 1.2j, 1.0, 0.0)
    print("achiral beschichtet gegen mie_coated.py (sigma_ext):")
    for d in [0.05, 0.2]:
        a = cross_sections(om, [1.0, 1.0 + d], [gold, (2.25, 1.0, 0.0)])[0]
        b = mc.qext(om, [1.0, 1.0 + d], [-11 + 1.2j, 2.25]) * np.pi * (1 + d) ** 2
        print(f"  d={d}: {a:.10f} {b:.10f}  rel {abs(a - b) / b:.1e}")
    print("chirale Schale aus demselben Material wie der chirale Kern = homogene chirale Kugel (mie_chiral.py):")
    ch = (2.25, 1.0, 0.2)
    for s in [+1, -1]:
        a = cross_sections(1.0, [0.7, 1.0], [ch, ch], s)[0]; b = mch.cross_sections(1.0, 2.25, 1.0, 0.2, s)[0]
        print(f"  s={s:+d}: {a:.10f} {b:.10f}  rel {abs(a - b) / b:.1e}")
    print("Symmetrie sigma_s(chi) = sigma_{-s}(-chi), Goldkern mit chiraler Schale:")
    a = cross_sections(om, [1.0, 1.05], [gold, (2.25, 1.0, 0.1)], +1)[0]; b = cross_sections(om, [1.0, 1.05], [gold, (2.25, 1.0, -0.1)], -1)[0]
    print(f"  {a:.12f} {b:.12f}")
    print("verlustfrei (Glaskern, chirale Schale): Extinktion = Streuung")
    e, sc = cross_sections(1.0, [0.8, 1.0], [(2.25, 1.0, 0.0), (2.0, 1.0, 0.15)], +1)
    print(f"  ext {e:.10f} sca {sc:.10f} rel {abs(e - sc) / e:.1e}")
    for d in [0.02, 0.05, 0.1]:
        sp = cross_sections(om, [1.0, 1 + d], [gold, (2.25, 1.0, 0.1)], +1)[0]; sm = cross_sections(om, [1.0, 1 + d], [gold, (2.25, 1.0, 0.1)], -1)[0]
        print(f"Goldkern, Schale d={d}, eps 2,25, chi 0,1: sigma+ {sp:.8f} sigma- {sm:.8f} CD {sp - sm:+.8e}")
