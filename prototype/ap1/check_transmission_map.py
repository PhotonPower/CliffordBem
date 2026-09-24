"""Pruefung der Transmissionsabbildung J an der exakten Fresnel-Loesung einer ebenen
Grenzflaeche z=0 (Medium 1 in z<0 = 'innen', n = e3; Medium 2 in z>0 = 'aussen')."""
import numpy as np
from cl3 import MV
om = 1.0
def pw(eps, mu, kvec, Evec, x):
    ph = np.exp(1j*np.dot(kvec, x))
    E = Evec*ph; H = np.cross(kvec, Evec)/(om*mu)*ph
    return E, H
def polar(kvec):
    kh = kvec/np.sqrt(np.dot(kvec,kvec))   # komplexe 'Einheits'-Richtung
    s = np.cross(kh, [0,0,1.]); s = s/np.sqrt(np.dot(s,s))
    p = np.cross(s, kh)
    return s, p
def run(eps1, mu1, eps2, mu2, theta):
    k1, k2 = om*np.sqrt(eps1*mu1), om*np.sqrt(eps2*mu2)
    kx = k2*np.sin(theta)
    kz2 = np.sqrt(k2**2-kx**2); kz1 = np.sqrt(k1**2-kx**2+0j)
    if kz1.imag < 0: kz1 = -kz1
    kin, kre, ktr = np.array([kx,0,-kz2]), np.array([kx,0,kz2]), np.array([kx,0,-kz1])
    sin_, pin = polar(kin); Ein = 0.6*sin_ + 0.8*pin
    Bre, Btr = polar(kre), polar(ktr)
    x0 = np.array([0.3,0.2,0.0])
    def fields(c):
        E2,H2 = pw(eps2,mu2,kin,Ein,x0)
        for j in range(2):
            e,h = pw(eps2,mu2,kre,Bre[j]*c[j],x0); E2=E2+e; H2=H2+h
        E1,H1 = np.zeros(3,complex), np.zeros(3,complex)
        for j in range(2):
            e,h = pw(eps1,mu1,ktr,Btr[j]*c[2+j],x0); E1=E1+e; H1=H1+h
        return E1,H1,E2,H2
    # lineares System fuer Tangentialstetigkeit (x,y-Komponenten von E und H)
    def resid(c):
        E1,H1,E2,H2 = fields(c)
        return np.concatenate([(E1-E2)[:2],(H1-H2)[:2]])
    r0 = resid(np.zeros(4)); A = np.array([resid(np.eye(4)[j])-r0 for j in range(4)]).T
    c = np.linalg.solve(A, -r0)
    E1,H1,E2,H2 = fields(c)
    I = MV.blade(7)
    h1 = MV.vec(np.sqrt(eps1)*E1) + I*MV.vec(np.sqrt(mu1)*H1)
    h2 = MV.vec(np.sqrt(eps2)*E2) + I*MV.vec(np.sqrt(mu2)*H2)
    re_, rm = eps1/eps2, mu1/mu2
    Jd = {0:1,1:np.sqrt(re_),2:np.sqrt(re_),4:1/np.sqrt(re_),3:1/np.sqrt(rm),5:np.sqrt(rm),6:np.sqrt(rm),7:1}
    Jh2 = MV([Jd[b]*h2.c[b] for b in range(8)])
    return np.sqrt(abs((h1-Jh2).norm2()))/np.sqrt(abs(h1.norm2()))
for (e1,m1,e2,m2) in [(2.25,1,1,1),(-11+1.2j,1,1.77,1),(4,2.5,1.3,0.7)]:
    for th in [0.05, 0.4, 1.1]:
        print(f"eps1={e1}, mu1={m1}, eps2={e2}, mu2={m2}, theta={th}:  |h1 - J h2|/|h1| = {run(e1,m1,e2,m2,th):.1e}")
