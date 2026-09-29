"""Mie-Loesung fuer geschichtete Kugeln mit chiralen (Pasteur-)Schichten in einem Aussenmedium (Standard Vakuum, eps2
waehlbar; seit v0.21 auch chiral, chi2), zirkular polarisierte Anregung (Helizitaet s).
Chirales Aussenmedium: Die einfallende Welle der Helizitaet s ist eine reine Beltrami-Welle, s = +1: W_A(j) mit k_A,
s = -1: W_B(j) mit k_B (E_inc = sum C W_s); gestreut alpha W_A(h) + beta W_B(h). Mit der einfallenden Welle interferiert nur
der Kanal derselben Helizitaet: sigma_ext^s = -sum Re(conj(C) 2 gamma) / (2 k_s^2), gamma = alpha bzw. beta;
sigma_sca = sum [2 |alpha|^2 / (2 k_A^2) + 2 |beta|^2 / (2 k_B^2)]. Achiraler Grenzfall: alpha = (a + b)/2, beta = (a - b)/2.

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


def coefficients(n, om, radii, media, s, eps2=1.0, chi2=0.0):
    """radii aufsteigend, media je Schicht (eps, mu, chi) von innen nach aussen; aussen Vakuum.
    Unbekannte: Kern (2), je Schale (4), aussen (a_n, b_n) gestreut; Aussenmedium eps2 (achiral, mu = 1)."""
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
            host = (eps2, 1.0, chi2)
            f1, f2 = _basis(n, om, host, r, 'h'); c = cols(L)          # achiral: M, N; chiral: W_A, W_B
            A[rows, c] += f1; A[rows, c + 1] += f2
            g1, g2 = _basis(n, om, host, r, 'j')
            rhs[rows] -= C * ((g1 + s * g2) if chi2 == 0 else (g1 if s > 0 else g2))   # E_inc = sum C (M + s N) bzw. C W_s
        else:
            c = cols(i + 1)
            for kd in ['j', 'y']:
                for f in _basis(n, om, media[i + 1], r, kd):
                    A[rows, c] += f; c += 1
    x = np.linalg.solve(A, rhs)
    return C, x[cols(L)], x[cols(L) + 1]


def cross_sections(om, radii, media, s=+1, nmax=None, eps2=1.0, chi2=0.0):
    """sigma_ext, sigma_sca (|E_inc|^2 = 2 wie in mie_chiral.py; Wellenzahl aussen k = omega sqrt(eps2), chiral k_A, k_B)."""
    k = om * np.sqrt(eps2); R = radii[-1]; nmax = nmax or int(abs(k) * R + 4 * (abs(k) * R) ** (1 / 3) + 12)
    if chi2 != 0:
        kA, kB = om * (np.sqrt(eps2) + chi2), om * (np.sqrt(eps2) - chi2); ks = kA if s > 0 else kB
        ext = sca = 0.0
        for n in range(1, nmax + 1):
            C, al, be = coefficients(n, om, radii, media, s, eps2, chi2)
            ext += -np.real(np.conj(C) * 2 * (al if s > 0 else be)) / (2 * ks ** 2)
            sca += 2 * abs(al) ** 2 / (2 * kA ** 2) + 2 * abs(be) ** 2 / (2 * kB ** 2)
        return ext, sca
    ext = sca = 0.0
    for n in range(1, nmax + 1):
        C, an, bn = coefficients(n, om, radii, media, s, eps2)
        ext += -np.real(np.conj(C) * an + np.conj(s * C) * bn); sca += abs(an) ** 2 + abs(bn) ** 2
    return ext / (2 * k ** 2), sca / (2 * k ** 2)


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
    print("achiral beschichtet in Wasser (eps2 = 1,7689) gegen mie_coated.py:")
    a = cross_sections(om, [1.0, 1.05], [gold, (2.25, 1.0, 0.0)], eps2=1.7689)[0]
    b = mc.qext(om, [1.0, 1.05], [-11 + 1.2j, 2.25], 1.7689) * np.pi * 1.05 ** 2
    print(f"  {a:.10f} {b:.10f}  rel {abs(a - b) / b:.1e}")
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
    print("chirales Aussenmedium: Grenzfall chi2 -> 0, Energieerhaltung, Spiegelsymmetrie, Kugel aus dem Aussenmedium:")
    for s_ in [+1, -1]:
        a = cross_sections(om, [1.0], [gold], s_, eps2=1.7689, chi2=1e-9)[0]; b = cross_sections(om, [1.0], [gold], s_, eps2=1.7689)[0]
        print(f"  s={s_:+d}, chi2 = 1e-9: {a:.10f} gegen achiral {b:.10f}  rel {abs(a - b) / b:.1e}")
    e, sc = cross_sections(1.0, [0.8, 1.0], [(2.25, 1.0, 0.0), (2.0, 1.0, 0.05)], +1, eps2=1.5, chi2=0.1)
    print(f"  verlustfrei, chirales Aussenmedium: ext {e:.10f} sca {sc:.10f} rel {abs(e - sc) / e:.1e}")
    a = cross_sections(om, [1.0], [gold], +1, eps2=1.7689, chi2=0.05)[0]; b = cross_sections(om, [1.0], [gold], -1, eps2=1.7689, chi2=-0.05)[0]
    print(f"  sigma_+(chi2) = sigma_-(-chi2): {a:.12f} {b:.12f}")
    a = cross_sections(om, [1.0], [(1.7689, 1.0, 0.05)], +1, eps2=1.7689, chi2=0.05)
    print(f"  Kugel aus dem Aussenmedium: ext {a[0]:.1e} sca {a[1]:.1e} (0 erwartet)")
    for d in [0.02, 0.05, 0.1]:
        sp = cross_sections(om, [1.0, 1 + d], [gold, (2.25, 1.0, 0.1)], +1)[0]; sm = cross_sections(om, [1.0, 1 + d], [gold, (2.25, 1.0, 0.1)], -1)[0]
        print(f"Goldkern, Schale d={d}, eps 2,25, chi 0,1: sigma+ {sp:.8f} sigma- {sm:.8f} CD {sp - sm:+.8e}")
