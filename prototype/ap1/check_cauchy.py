"""Numerische Pruefung der Cauchy-Darstellungsformel auf der Einheitskugel.
C_k h(x) = - int_Gamma E_k(x-y) n(y) h(y) dS(y)
Erwartung:  innere Loesung F:   C_k(F|_G) = F  innen,  0 aussen
            aeussere (strahlende) Loesung F: C_k(F|_G) = 0 innen, -F aussen
"""
import numpy as np
from cl3 import MV, SIGN

def Ek(z, k):
    r = np.linalg.norm(z); zh = z/r
    Phi = np.exp(1j*k*r)/(4*np.pi*r)
    return MV.vec(Phi*(1/r - 1j*k)*zh) + MV.blade(0, -1j*k*Phi)

def sphere_quad(nt, nphi, R=1.0):
    t, wt = np.polynomial.legendre.leggauss(nt)          # t = cos theta
    ph = 2*np.pi*np.arange(nphi)/nphi
    pts, wts = [], []
    for ti, wi in zip(t, wt):
        st = np.sqrt(1-ti**2)
        for p in ph:
            n = np.array([st*np.cos(p), st*np.sin(p), ti])
            pts.append(R*n); wts.append(R**2*wi*2*np.pi/nphi)
    return np.array(pts), np.array(wts)

def cauchy(h, x, k, pts, wts):
    acc = MV()
    for y, wq in zip(pts, wts):
        acc = acc + (Ek(x-y, k) * MV.vec(y/np.linalg.norm(y)) * h(y)) * (-wq)
    return acc

rng = np.random.default_rng(1)
c = MV(rng.normal(size=8) + 1j*rng.normal(size=8))       # beliebiger konstanter Multivektor
pts, wts = sphere_quad(40, 80)
xin, xout = np.array([0.2,-0.1,0.3]), np.array([1.8,0.5,-0.4])

for k in [2.0, 2.0+0.3j, 1e-3]:
    x1 = np.array([2.5, -1.0, 0.7])                        # Quelle aussen -> innere Loesung
    x0 = np.array([0.1, 0.25, -0.2])                       # Quelle innen  -> aeussere Loesung
    Fint = lambda y: Ek(y-x1, k)*c
    Fext = lambda y: Ek(y-x0, k)*c
    e1 = (cauchy(Fint, xin, k, pts, wts) - Fint(xin)).norm2()**.5 / Fint(xin).norm2()**.5
    e2 = cauchy(Fint, xout, k, pts, wts).norm2()**.5 / Fint(xout).norm2()**.5
    e3 = cauchy(Fext, xin, k, pts, wts).norm2()**.5 / Fext(xin).norm2()**.5
    e4 = (cauchy(Fext, xout, k, pts, wts) + Fext(xout)).norm2()**.5 / Fext(xout).norm2()**.5
    print(f"k={k}:  innen/int {abs(e1):.1e}  aussen/int {abs(e2):.1e}  innen/ext {abs(e3):.1e}  aussen/ext {abs(e4):.1e}")

# Maxwell-Feld: ebene Welle in Medium (eps=2.25, mu=1), omega=1.3  -> innere Loesung
eps, mu, om = 2.25, 1.0, 1.3; k = om*np.sqrt(eps*mu)
def Fpw(y):
    ph = np.exp(1j*k*y[2])
    return MV.vec(np.sqrt(eps)*ph*np.array([1,0,0])) + MV.blade(7)*MV.vec(np.sqrt(mu)*np.sqrt(eps/mu)*ph*np.array([0,1,0]))
e = (cauchy(Fpw, xin, k, pts, wts) - Fpw(xin)).norm2()**.5/Fpw(xin).norm2()**.5
print(f"ebene Maxwell-Welle, innen: rel. Fehler {abs(e):.1e}")
