"""Mie-Loesung fuer eine chirale (Pasteur-)Kugel in Vakuum, zirkular polarisierte Anregung.
Konstitutiv (wie im Kern, AP 1): D = eps E + i chi H,  B = mu H - i chi E  (e^{-i omega t}, eps0 = mu0 = 1).
Innen Beltrami-Felder: W_A = M(k_A) + N(k_A), rot W_A = +k_A W_A, k_A = omega (sqrt(eps mu) + chi);
                       W_B = M(k_B) - N(k_B), rot W_B = -k_B W_B, k_B = omega (sqrt(eps mu) - chi).
Aussen: E_inc = (e_x + i s e_y) e^{ikz} = sum_n C_n (M_n + s N_n) (m = s), C_n = i^n sqrt(4 pi (2n+1)).
Je n ein 4x4-System aus der Stetigkeit von E_tan und H_tan (Komponenten entlang X_nm und r x X_nm).
Normierung der Querschnitte ueber den achiralen Grenzfall (Vergleich mit mie.py) und die Energiebilanz."""
import numpy as np
from scipy.special import spherical_jn, spherical_yn
def zfun(n, u, kind):
    j = spherical_jn(n, u); jp = spherical_jn(n, u, derivative=True)
    if kind == 'j': return j, j + u*jp                        # z_n, psi'(u) = d/du [u z_n] = z_n + u z_n'
    y = spherical_yn(n, u); yp = spherical_yn(n, u, derivative=True)
    h = j + 1j*y; hp = jp + 1j*yp; return h, h + u*hp
def jcx(n, u):
    """sphaerische Bessel j_n fuer komplexes Argument ueber Rekursion aus sin/cos (fuer kleine n stabil genug)"""
    u = complex(u)
    if abs(u.imag) < 1e-14 and False: pass
    # Aufwaertsrekursion ist fuer |u| > n stabil; sonst Rueckwaertsrekursion (Miller)
    N = n + 30 + int(abs(u))
    jn = np.zeros(N + 2, complex); jn[N + 1] = 0; jn[N] = 1e-30
    for l in range(N, 0, -1): jn[l - 1] = (2*l + 1)/u*jn[l] - jn[l + 1]
    j0 = np.sin(u)/u; jn *= j0/jn[0]
    jv = jn[n]; jd = jn[n - 1] - (n + 1)/u*jn[n] if n > 0 else -jn[1]
    return jv, jv + u*jd
def coefficients(n, om, eps, mu, chi, s, a=1.0):
    k0 = om; u0 = k0*a
    kA, kB = om*(np.sqrt(eps*mu) + chi), om*(np.sqrt(eps*mu) - chi)
    eta = np.sqrt(mu/eps); hA, hB = -1j/eta, 1j/eta
    jA, pA = jcx(n, kA*a); jB, pB = jcx(n, kB*a)
    j0, p0 = zfun(n, u0, 'j'); h0, ph0 = zfun(n, u0, 'h')
    uA, uB = kA*a, kB*a
    C = 1j**n*np.sqrt(4*np.pi*(2*n + 1))
    # Unbekannte x = (a_n, b_n, alpha, beta)
    A = np.array([[h0, 0, -jA, -jB],
                  [0, ph0/u0, -pA/uA, pB/uB],
                  [0, -1j*h0, -hA*jA, -hB*jB],
                  [-1j*ph0/u0, 0, -hA*pA/uA, hB*pB/uB]], complex)
    rhs = -np.array([C*j0, C*s*p0/u0, -1j*C*s*j0, -1j*C*p0/u0], complex)
    x = np.linalg.solve(A, rhs)
    return C, x
def cross_sections(om, eps, mu=1.0, chi=0.0, s=+1, a=1.0, nmax=None):
    nmax = nmax or int(om*a + 4*(om*a)**(1/3) + 12)
    ext = 0; sca = 0
    for n in range(1, nmax + 1):
        C, (an, bn, al, be) = coefficients(n, om, eps, mu, chi, s, a)
        ext += -np.real(np.conj(C)*an + np.conj(s*C)*bn)
        sca += abs(an)**2 + abs(bn)**2
    # Normierung: |E_inc|^2 = 2; Fernfeld der normierten VSH: |M|^2 r^2 -> |c|^2/k^2 (k = om)
    k = om
    return ext/(2*k**2), sca/(2*k**2)
if __name__ == '__main__':
    import mie
    print("achiraler Grenzfall gegen mie.py (Q_ext = sigma_ext / (pi a^2)):")
    for om, e in [(1.0, 2.25), (0.5, -11+1.2j), (0.3, -2.2+0.3j)]:
        se, ss = cross_sections(om, e, chi=0.0, s=+1)
        print(f"  omega={om}, eps={e}: chiral-Code {se/np.pi:.6f}, mie.py {mie.qext(om, e):.6f}  (Q_sca {ss/np.pi:.6f})")
    print("verlustfrei chiral: Extinktion = Streuung?")
    for chi in [0.1, 0.3]:
        for s in [+1, -1]:
            se, ss = cross_sections(1.0, 2.25, chi=chi, s=s)
            print(f"  chi={chi}, s={s:+d}: Q_ext {se/np.pi:.8f}, Q_sca {ss/np.pi:.8f}, Diff {abs(se-ss)/se:.1e}")
    print("Symmetrie: sigma_s(chi) = sigma_{-s}(-chi):")
    a1 = cross_sections(1.0, 2.25, chi=0.2, s=+1)[0]; a2 = cross_sections(1.0, 2.25, chi=-0.2, s=-1)[0]
    print(f"  {a1:.10f} {a2:.10f}")
