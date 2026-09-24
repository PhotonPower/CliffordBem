"""Transmissionsabbildung J fuer chirale (Pasteur-)Medien und Pruefung an der exakten
Fresnel-Loesung einer ebenen Grenzflaeche (Medium 1 chiral in z<0, Medium 2 in z>0).
Konstitutiv: D = eps E + i chi H,  B = mu H - i chi E.
Normal: C_j (E_n, H_n) stetig mit C_j = [[eps_j, i chi_j], [-i chi_j, mu_j]]."""
import numpy as np
from cl3 import MV
from fields import gp, vec, PSEUDO
def Jchiral(n, m1, m2, a=None, b=None):
    """J als 8x8 (Blade-Basis) zur Normalen n; m = (eps, mu, chi). Standardwahl der Hilfsparameter."""
    e1, u1, c1 = m1; e2, u2, c2 = m2
    se, sm = np.sqrt(e1)/np.sqrt(e2), np.sqrt(u1)/np.sqrt(u2)
    a = sm if a is None else a; b = se if b is None else b
    C1 = np.array([[e1, 1j*c1], [-1j*c1, u1]]); C2 = np.array([[e2, 1j*c2], [-1j*c2, u2]])
    D1 = np.diag([np.sqrt(e1), np.sqrt(u1)]); D2 = np.diag([np.sqrt(e2), np.sqrt(u2)])
    Nn = D1@np.linalg.solve(C1, C2)@np.linalg.inv(D2)     # (e_n, h_n)_1 = Nn (e_n, h_n)_2
    Jm = np.zeros((8, 8), complex)
    for bb in range(8):
        H = np.zeros(8, complex); H[bb] = 1
        v = H[[1, 2, 4]]; Ib = H[[6, 5, 3]]*np.array([1, -1, 1])
        vn, bn = n@v, n@Ib
        vn2, bn2 = Nn@np.array([vn, bn])
        v2 = se*(v - vn*n) + vn2*n; b2 = sm*(Ib - bn*n) + bn2*n
        out = np.zeros(8, complex); out[0] = a*H[0]; out[7] = b*H[7]
        out[[1, 2, 4]] = v2; out[[6, 5, 3]] = b2*np.array([1, -1, 1]); Jm[:, bb] = out
    return Jm
def P(sign):  # Helizitaetsprojektor (1 + sign i I)/2 als Multivektor
    p = np.zeros(8, complex); p[0] = 0.5; p[7] = 0.5j*sign; return p
def fresnel_check(m1, m2, theta, om=1.0, x0=np.array([0.3, 0.2, 0.0])):
    e1, u1, c1 = m1; e2, u2, c2 = m2
    k2 = om*np.sqrt(e2*u2); kx = k2*np.sin(theta)
    kp, km = om*(np.sqrt(e1)*np.sqrt(u1) - c1), om*(np.sqrt(e1)*np.sqrt(u1) + c1)
    def pw_medium2(dvec, pvec):   # achiral: F = sqrt(eps)(p + I d x p) e^{ik d.x}
        ph = np.exp(1j*k2*(dvec@x0)); return (vec(np.sqrt(e2)*pvec*ph) + gp(PSEUDO, vec(np.sqrt(e2)*np.cross(dvec, pvec)*ph)))
    def pw_chiral(sign, kk, cvec):
        kz = np.sqrt(kk**2 - kx**2 + 0j); kz = kz if kz.imag >= 0 else -kz
        d = np.array([kx, 0, -kz])/kk                         # nach -z laufend (z<0 Medium)
        cvec = cvec - (d@cvec)*d                               # transversal (d.d = 1 komplex)
        g = gp(P(sign), (np.eye(8)[0] + vec(d)).astype(complex)); g = gp(g, vec(cvec))
        return g*np.exp(1j*kk*(d@x0))
    din = np.array([kx, 0, -np.sqrt(k2**2 - kx**2)])/k2; dre = din*np.array([1, 1, -1])
    s_in = np.cross(din, [0, 0, 1.]); s_in /= np.linalg.norm(s_in); p_in = np.cross(s_in, din)
    s_re = np.cross(dre, [0, 0, 1.]); s_re /= np.linalg.norm(s_re); p_re = np.cross(s_re, dre)
    Finc = pw_medium2(din, 0.6*s_in + 0.8*p_in)
    basis2 = [pw_medium2(dre, s_re), pw_medium2(dre, p_re)]
    basis1 = [pw_chiral(+1, kp, np.array([1., 0.3, 0])), pw_chiral(-1, km, np.array([1., -0.4, 0.2]))]
    def fields(F, eps, mu):
        E = F[[1, 2, 4]]/np.sqrt(eps); Hm = -gp(PSEUDO, F); H = Hm[[1, 2, 4]]/np.sqrt(mu); return E, H
    def resid(c):
        F2 = Finc + c[0]*basis2[0] + c[1]*basis2[1]; F1 = c[2]*basis1[0] + c[3]*basis1[1]
        E1, H1 = fields(F1, e1, u1); E2, H2 = fields(F2, e2, u2)
        return np.concatenate([(E1 - E2)[:2], (H1 - H2)[:2]]), F1, F2
    r0 = resid(np.zeros(4))[0]; A = np.array([resid(np.eye(4)[j])[0] - r0 for j in range(4)]).T
    c = np.linalg.solve(A, -r0); _, F1, F2 = resid(c)
    J = Jchiral(np.array([0, 0, 1.]), m1, m2)
    err = np.linalg.norm(F1 - J@F2)/np.linalg.norm(F1)
    # Kontrolle der konstitutiven Normalbedingung direkt
    E1, H1 = fields(F1, e1, u1); E2, H2 = fields(F2, e2, u2)
    Dn = abs((e1*E1[2] + 1j*c1*H1[2]) - (e2*E2[2] + 1j*c2*H2[2])); Bn = abs((u1*H1[2] - 1j*c1*E1[2]) - (u2*H2[2] - 1j*c2*E2[2]))
    return err, Dn, Bn
if __name__ == '__main__':
    for m1 in [(2.25, 1.0, 0.1), (2.25, 1.0, 0.4), (-11+1.2j, 1.0, 0.3), (3.0, 1.5, 0.2+0.05j)]:
        for th in [0.05, 0.5, 1.1]:
            err, Dn, Bn = fresnel_check(m1, (1.0, 1.0, 0.0), th)
            print(f"eps1,mu1,chi1={m1}, theta={th}: |h1 - J h2|/|h1| = {err:.1e}   (Kontrolle D_n: {Dn:.1e}, B_n: {Bn:.1e})")
