"""GLMT-Referenz (v0.39): Kraft auf eine homogene Kugel in einem fokussierten Strahl (Richards-Wolf) oder Gaussstrahl,
derselbe Strahl wie BeamField in C++ (gleiche Quadraturknoten, gleiche Normierung auf die Leistung 1).

Weg (unabhaengig von der BEM):
  1. Strahlkoeffizienten: Der Strahl ist eine endliche Summe ebener Wellen; jede Teilwelle e^{ik k.r} e hat die geschlossene
     Entwicklung c_M = 4 pi i^n X*_nm(k).e,  c_N = alpha 4 pi i^n (k x X*_nm(k)).e  nach M_nm = j_n X_nm, N_nm = (1/k) rot M_nm
     (X_nm = L Y_nm / sqrt(n(n+1)), L = -i r x grad); alpha wird im Selbsttest bestimmt, nicht vorausgesetzt.
  2. Mie (Bohren-Huffman): einlaufend c_M, c_N -> auslaufend (h_n) s_M = -b_n c_M, s_N = -a_n c_N.
  3. Kraft: Maxwellscher Spannungstensor auf einer Kugel um das Teilchen; einfallendes Feld direkt aus dem Strahl, Streufeld aus
     der Multipolsumme, H_s = (k/(i omega mu)) sum (s_M N + s_N M).
Einheiten: Laengen in Teilchenradien, omega = 2 pi R / lambda0, eps0 = mu0 = 1; Effizienz Q = F / n (Leistung 1).

Aufruf: python3 tools/glmt.py --selftest
        python3 tools/glmt.py --radius 0.25 --n-particle 1.59 --nbg 1.33 --lambda 1.064 --NA 1.2 --f0 1 --pol x --z "-1:2:13"
"""
import argparse, os, sys
import numpy as np
from scipy.special import sph_harm_y, spherical_jn, spherical_yn
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mie_coated as mc

ALPHA = None   # Vorfaktor der N-Koeffizienten der ebenen Welle, im Selbsttest bestimmt


# ---------------------------------------------------------------- vektorielle Kugelflaechenfunktionen und Kugelwellen
def _frame(th, ph):
    st, ct, sp, cp = np.sin(th), np.cos(th), np.sin(ph), np.cos(ph)
    rh = np.stack([st * cp, st * sp, ct], -1); thh = np.stack([ct * cp, ct * sp, -st], -1); phh = np.stack([-sp, cp, np.zeros_like(ph)], -1)
    return rh, thh, phh


def Xnm(n, m, th, ph):
    """X_nm = L Y_nm / sqrt(n(n+1)), L Y = -i phi_hat dY/dtheta - theta_hat m Y / sin(theta) (kartesisch, Form (..., 3))."""
    h = 1e-6
    Y = sph_harm_y(n, m, th, ph)
    dY = (sph_harm_y(n, m, th + h, ph) - sph_harm_y(n, m, th - h, ph)) / (2 * h)
    _, thh, phh = _frame(th, ph)
    LY = -1j * dY[..., None] * phh - (m * Y / np.sin(th))[..., None] * thh
    return LY / np.sqrt(n * (n + 1)), Y


def zn(n, x, kind):
    if kind == 'j': return spherical_jn(n, x), spherical_jn(n, x, derivative=True)
    return spherical_jn(n, x) + 1j * spherical_yn(n, x), spherical_jn(n, x, derivative=True) + 1j * spherical_yn(n, x, derivative=True)


def MN(n, m, k, pts, kind):
    """M = z_n(kr) X_nm,  N = (1/k) rot M = i sqrt(n(n+1)) z_n/(kr) Y r_hat + (1/(kr)) d(kr z_n)/d(kr) (r_hat x X)."""
    r = np.linalg.norm(pts, axis=-1); th = np.arccos(np.clip(pts[..., 2] / r, -1, 1)); ph = np.arctan2(pts[..., 1], pts[..., 0])
    X, Y = Xnm(n, m, th, ph); rh, _, _ = _frame(th, ph)
    z, dz = zn(n, k * r, kind); x = k * r
    M = z[..., None] * X
    N = (1j * np.sqrt(n * (n + 1)) * z / x * Y)[..., None] * rh + ((z + x * dz) / x)[..., None] * np.cross(rh, X)
    return M, N


# ---------------------------------------------------------------- Strahl wie BeamField (C++)
def gauss_legendre_ab(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (a + b) + 0.5 * (b - a) * x, 0.5 * (b - a) * w


class Beam:
    """Ebene Teilwellen dir (K,3), amp (K,3) = a dOmega, dw (K,), Brennpunkt focus; auf Leistung 1 normiert."""
    def __init__(self, eps, mu, omega, focus, dirs, amps, dws):
        self.eps, self.mu, self.k = eps, mu, omega * np.sqrt(eps * mu)
        self.focus, self.dir, self.amp, self.dw = np.asarray(focus, float), np.asarray(dirs), np.asarray(amps, complex), np.asarray(dws)
        P = 0.5 * np.sqrt(eps / mu) * (2 * np.pi / self.k) ** 2 * np.sum(np.sum(abs(self.amp) ** 2, 1) / self.dw)
        self.amp /= np.sqrt(P)

    @staticmethod
    def focused(eps, mu, omega, focus, NA, f0, pol, nt=40, nphi=80):
        n = np.sqrt(eps * mu); tmax = np.arcsin(NA / n); smax = np.sin(tmax)
        th, wth = gauss_legendre_ab(nt, 0.0, tmax); dirs, amps, dws = [], [], []
        for t, wt in zip(th, wth):
            st, ct = np.sin(t), np.cos(t); f = np.exp(-st ** 2 / (f0 ** 2 * smax ** 2)) * np.sqrt(ct)
            for j in range(nphi):
                p = 2 * np.pi * (j + 0.5) / nphi; cp, sp = np.cos(p), np.sin(p)
                phi = np.array([-sp, cp, 0.0]); the = np.array([ct * cp, ct * sp, -st])
                er = pol[0] * cp + pol[1] * sp; ep = -pol[0] * sp + pol[1] * cp; w = wt * st * 2 * np.pi / nphi
                dirs.append([st * cp, st * sp, ct]); amps.append(w * f * (ep * phi + er * the)); dws.append(w)
        return Beam(eps, mu, omega, focus, dirs, amps, dws)

    @staticmethod
    def plane(eps, mu, omega, d, pol):          # eine ebene Welle mit Amplitude pol (fuer den Selbsttest; keine Normierung)
        b = Beam.__new__(Beam); b.eps, b.mu, b.k = eps, mu, omega * np.sqrt(eps * mu)
        b.focus = np.zeros(3); b.dir = np.array([d], float); b.amp = np.array([pol], complex); b.dw = np.array([1.0]); return b

    def fields(self, pts):
        ph = np.exp(1j * self.k * (pts - self.focus) @ self.dir.T)              # (P, K)
        E = ph @ self.amp; H = np.sqrt(self.eps / self.mu) * (ph @ np.cross(self.dir, self.amp))
        return E, H

    def coefficients(self, nmax):
        """c_M, c_N je (n, m) aus den Teilwellen (exakt fuer die endliche Summe)."""
        th = np.arccos(np.clip(self.dir[:, 2], -1, 1)); ph = np.arctan2(self.dir[:, 1], self.dir[:, 0])
        shift = np.exp(-1j * self.k * self.dir @ self.focus)                   # Brennpunkt relativ zum Kugelmittelpunkt
        cM, cN = {}, {}
        for n in range(1, nmax + 1):
            for m in range(-n, n + 1):
                X, _ = Xnm(n, m, th, ph); Xc = np.conj(X)
                cM[n, m] = 4 * np.pi * 1j ** n * np.sum(shift * np.sum(Xc * self.amp, 1))
                cN[n, m] = ALPHA * 4 * np.pi * 1j ** n * np.sum(shift * np.sum(np.cross(self.dir, Xc) * self.amp, 1))
        return cM, cN


# ---------------------------------------------------------------- Streufeld und Kraft
def scattered_fields(cM, cN, a, b, k, omega, mu, pts):
    E = np.zeros(pts.shape, complex); H = np.zeros(pts.shape, complex)
    for (n, m), c in cM.items():
        if n > len(a): continue
        sM, sN = -b[n - 1] * c, -a[n - 1] * cN[n, m]
        M, N = MN(n, m, k, pts, 'h')
        E += sM * M + sN * N; H += (k / (1j * omega * mu)) * (sM * N + sN * M)
    return E, H


def stress_force(E, H, nrm, w, eps, mu):
    En = np.sum(E * nrm, 1); Hn = np.sum(H * nrm, 1)
    u = 0.5 * (eps * np.sum(abs(E) ** 2, 1) + mu * np.sum(abs(H) ** 2, 1))
    T = 0.5 * (np.real(eps * E * np.conj(En)[:, None] + mu * H * np.conj(Hn)[:, None]) - u[:, None] * nrm)
    return np.sum(T * w[:, None], 0)


def sphere_quadrature(R, nt=40):
    ct, wt = np.polynomial.legendre.leggauss(nt); nphi = 2 * nt
    P, Nv, W = [], [], []
    for c, w in zip(ct, wt):
        s = np.sqrt(1 - c * c)
        for j in range(nphi):
            p = 2 * np.pi * (j + 0.5) / nphi; v = np.array([s * np.cos(p), s * np.sin(p), c]); P.append(R * v); Nv.append(v); W.append(w * 2 * np.pi / nphi * R * R)
    return np.array(P), np.array(Nv), np.array(W)


def force(beam, eps_p, omega, nmax=None, Rs=1.3, nt=40):
    k = beam.k
    nmie = nmax or int(np.ceil(k + 4 * k ** (1 / 3) + 6))
    nn, a, b = mc.coeffs_layered(omega, [1.0], [eps_p], beam.eps, nmie)
    cM, cN = beam.coefficients(nmie)
    P, Nv, W = sphere_quadrature(Rs, nt)
    Ei, Hi = beam.fields(P); Es, Hs = scattered_fields(cM, cN, a, b, k, omega, beam.mu, P)
    return stress_force(Ei + Es, Hi + Hs, Nv, W, beam.eps, beam.mu)


# ---------------------------------------------------------------- Selbsttest
def selftest():
    global ALPHA
    rng = np.random.default_rng(3); k = 1.7; ok = True
    # 1. N = (1/k) rot M (numerische Rotation)
    x0 = np.array([[0.4, -0.3, 0.9]]); h = 1e-5; err = 0
    for (n, m) in [(1, 0), (2, 1), (3, -2), (5, 4)]:
        M, N = MN(n, m, k, x0, 'h'); rot = np.zeros(3, complex)
        d = [(MN(n, m, k, x0 + h * e, 'h')[0] - MN(n, m, k, x0 - h * e, 'h')[0])[0] / (2 * h) for e in np.eye(3)]   # d[i][j] = d_i M_j
        rot = np.array([d[1][2] - d[2][1], d[2][0] - d[0][2], d[0][1] - d[1][0]])
        err = max(err, np.linalg.norm(rot / k - N[0]) / np.linalg.norm(N[0]))
    print(f"  N = rot M / k: relativer Fehler {err:.1e}"); ok &= err < 1e-5   # numerische Ableitung von Y und M
    # 2. Entwicklung der ebenen Welle: alpha bestimmen
    d = np.array([0.3, -0.4, 0.866]); d /= np.linalg.norm(d); e = np.cross(d, [0, 0, 1.0]); e /= np.linalg.norm(e); e = e + 0.5j * np.cross(d, e)
    pts = rng.normal(size=(5, 3)) * 0.8; Eex = np.exp(1j * k * pts @ d)[:, None] * e
    best = None
    for al in [1, -1, 1j, -1j]:
        ALPHA = al; B = Beam.plane(1.0, 1.0, k, d, e); cM, cN = B.coefficients(22)
        E = np.zeros(pts.shape, complex)
        for (n, m), c in cM.items():
            M, N = MN(n, m, k, pts, 'j'); E += c * M + cN[n, m] * N
        r = np.linalg.norm(E - Eex) / np.linalg.norm(Eex)
        if best is None or r < best[1]: best = (al, r)
    ALPHA = best[0]; print(f"  ebene Welle: alpha = {ALPHA}, Rekonstruktionsfehler {best[1]:.1e}"); ok &= best[1] < 1e-8
    # 3. Kraftweg gegen Mie: ebene Welle auf eine Kugel (sigma_pr 1/2 eps |E0|^2)
    import mie_force
    om, eps2, epsp = 0.8, 1.7689, 2.5281
    # schraeg einfallend (entlang z waere sin(theta) = 0 in X_nm; die Gauss-Knoten der Strahlen treffen theta = 0 nie)
    d = np.array([0.3, -0.4, 0.866]); d /= np.linalg.norm(d); e = np.cross(d, [0, 0, 1.0]); e /= np.linalg.norm(e)
    B = Beam.plane(eps2, 1.0, om, d, e); F = force(B, epsp, om)
    ref = mie_force.force_z(om, [1.0], [epsp], eps2); Fd = F @ d
    print(f"  ebene Welle auf Kugel: F.d {Fd:.10f}, Mie {ref:.10f}, quer {np.linalg.norm(F - Fd * d):.1e}"); ok &= abs(Fd / ref - 1) < 1e-8
    # 4. Strahl: Feld aus den Koeffizienten im Teilchengebiet = Strahl direkt
    Bf = Beam.focused(1.7689, 1.0, 0.39, np.array([0.3, -0.2, 0.5]), 1.2, 1.0, np.array([1.0, 0.0]))
    cM, cN = Bf.coefficients(24); pts = rng.normal(size=(6, 3)) * 0.5; Ed, _ = Bf.fields(pts)
    E = np.zeros(pts.shape, complex)
    for (n, m), c in cM.items():
        M, N = MN(n, m, Bf.k, pts, 'j'); E += c * M + cN[n, m] * N
    r = np.linalg.norm(E - Ed) / np.linalg.norm(Ed); print(f"  Strahlkoeffizienten: Feld im Teilchengebiet {r:.1e}"); ok &= r < 1e-8
    print("alle Pruefungen bestanden" if ok else "PRUEFUNG FEHLGESCHLAGEN"); return ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--selftest', action='store_true'); ap.add_argument('--radius', type=float, default=0.25)
    ap.add_argument('--n-particle', type=float, default=1.59); ap.add_argument('--nbg', type=float, default=1.33)
    ap.add_argument('--lambda', dest='lam', type=float, default=1.064); ap.add_argument('--NA', type=float, default=1.2)
    ap.add_argument('--f0', type=float, default=1.0); ap.add_argument('--pol', default='x')
    ap.add_argument('--z', default=None); ap.add_argument('--x', default=None); ap.add_argument('--at-z', type=float, default=0.0)
    ap.add_argument('--csv', default=None)
    A = ap.parse_args()
    if not selftest(): sys.exit(1)
    if A.selftest: sys.exit(0)
    R = A.radius; om = 2 * np.pi * R / A.lam; eps2 = A.nbg ** 2; epsp = A.n_particle ** 2
    s2 = 1 / np.sqrt(2); pol = {'x': [1, 0], 'y': [0, 1], 'circ+': [s2, 1j * s2], 'circ-': [s2, -1j * s2]}[A.pol]
    rng = lambda s: np.linspace(float(s.split(':')[0]), float(s.split(':')[1]), int(s.split(':')[2]))
    pos = [np.array([0, 0, z]) for z in rng(A.z)] if A.z else [np.array([x, 0, A.at_z]) for x in rng(A.x)]
    out = open(A.csv, 'w') if A.csv else None
    if out: out.write("x_um,y_um,z_um,Qx,Qy,Qz\n")
    for p in pos:
        B = Beam.focused(eps2, 1.0, om, -p / R, A.NA, A.f0, np.array(pol, complex))
        Q = force(B, epsp, om) / A.nbg
        print(f"  x {p[0]:6.3f}  z {p[2]:6.3f}   Q = ({Q[0]:+.5e}, {Q[1]:+.5e}, {Q[2]:+.5e})", flush=True)
        if out: out.write(f"{p[0]},{p[1]},{p[2]},{Q[0]},{Q[1]},{Q[2]}\n"); out.flush()
