"""Pasteur-Medium: D = eps E + i chi H,  B = mu H - i chi E  (Zeit e^{-i w t}).
Behauptung: F = sqrt(eps)E + I sqrt(mu)H erfuellt  nabla F = (ik + w chi I) F,
und P_pm F erfuellen (nabla - i k_pm) F_pm = 0 mit k_pm = w (sqrt(eps mu) -+ chi)."""
import sympy as sp
from cl3 import MV, nabla_left
x,y,z = X = sp.symbols('x y z', real=True)
w, eps, mu, chi = sp.symbols('omega epsilon mu chi', positive=True)
Ih = sp.I; S = lambda b,c: MV.blade(b,c,True); Ip = S(7,1)
k = w*sp.sqrt(eps*mu)
Pp = (S(0,1) + Ih*Ip)*sp.Rational(1,2); Pm = (S(0,1) - Ih*Ip)*sp.Rational(1,2)
def curl(v): return [sp.diff(v[2],y)-sp.diff(v[1],z), sp.diff(v[0],z)-sp.diff(v[2],x), sp.diff(v[1],x)-sp.diff(v[0],y)]
def split(F):   # F -> (E, H)
    Ev = [F.c[1]/sp.sqrt(eps), F.c[2]/sp.sqrt(eps), F.c[4]/sp.sqrt(eps)]
    Bi = F.grade(2); Hm = (Ip*Bi)*(-1)   # grade2 = I sqrt(mu) H  ->  sqrt(mu) H = -I * grade2
    Hv = [Hm.c[1]/sp.sqrt(mu), Hm.c[2]/sp.sqrt(mu), Hm.c[4]/sp.sqrt(mu)]
    return Ev, Hv
for sgn, P, name in [(1,Pp,'+'),(-1,Pm,'-')]:
    kpm = w*(sp.sqrt(eps*mu) - sgn*chi)
    # Ansatz: beliebiger transversaler Vektor, Projektion auf Helizitaet
    F0 = ((S(0,1) + MV.vec([0,0,1],True))*MV.vec([1,0,0],True))*sp.exp(Ih*kpm*z)   # e3 F0 = F0
    F = P*F0
    assert (nabla_left(F,X) - Ih*kpm*F).simplify().is_zero()
    Ev, Hv = split(F)
    D = [eps*e + Ih*chi*h for e,h in zip(Ev,Hv)]; B = [mu*h - Ih*chi*e for e,h in zip(Ev,Hv)]
    r1 = [sp.simplify(c - Ih*w*b) for c,b in zip(curl(Ev),B)]
    r2 = [sp.simplify(c + Ih*w*d) for c,d in zip(curl(Hv),D)]
    print(f"Helizitaet {name}: (nabla - i k_{name})F = 0 und chirale Maxwell-Gl. erfuellt:", all(v==0 for v in r1+r2))
