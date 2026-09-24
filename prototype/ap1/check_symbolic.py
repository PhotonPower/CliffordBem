import sympy as sp
from cl3 import *
x,y,z = X = sp.symbols('x y z', real=True)
k, w, eps, mu = sp.symbols('k omega epsilon mu', positive=True)
Ih = sp.I
S = lambda b,c: MV.blade(b,c,True)
Ipse = S(7,1)

# --- (A) Grad-Zerlegung von nabla F fuer allgemeines F = a0 + a + I b + I b0 ---
fs = [sp.Function(n)(*X) for n in ['al','a1','a2','a3','b1','b2','b3','be']]
a = MV.vec(fs[1:4],True); b = MV.vec(fs[4:7],True)
F = S(0,fs[0]) + a + Ipse*b + Ipse*S(0,fs[7])
DF = nabla_left(F,X)
div = lambda v: sum(sp.diff(v[i],X[i]) for i in range(3))
def curl(v): return [sp.diff(v[2],y)-sp.diff(v[1],z), sp.diff(v[0],z)-sp.diff(v[2],x), sp.diff(v[1],x)-sp.diff(v[0],y)]
grad = lambda f: [sp.diff(f,t) for t in X]
av, bv = fs[1:4], fs[4:7]
pred = S(0,div(av)) + MV.vec([g - c for g,c in zip(grad(fs[0]),curl(bv))],True) \
     + Ipse*MV.vec([c+g for c,g in zip(curl(av),grad(fs[7]))],True) + Ipse*S(0,div(bv))
assert (DF - pred).simplify().is_zero()
print("(A) nabla F = div a + (grad al - curl b) + I(curl a + grad be) + I div b : OK")

# --- (B) Maxwell-Aequivalenz fuer ebene Welle und Dipolfeld ---
kk = w*sp.sqrt(eps*mu)
E = [sp.exp(Ih*kk*z), 0, 0]; H = [0, sp.sqrt(eps/mu)*sp.exp(Ih*kk*z), 0]
Fpw = MV.vec([sp.sqrt(eps)*c for c in E],True) + Ipse*MV.vec([sp.sqrt(mu)*c for c in H],True)
assert (nabla_left(Fpw,X) - Ih*kk*Fpw).simplify().is_zero()
print("(B1) ebene Welle: nabla F = ik F : OK")
# elektrischer Dipol p=e3: H ~ curl(Phi e3), E = i/(omega eps) curl H
r = sp.sqrt(x**2+y**2+z**2)
Phi = sp.exp(Ih*kk*r)/(4*sp.pi*r)
Hd = curl([0,0,Phi]); Ed = [Ih/(w*eps)*c for c in curl(Hd)]
Fd = MV.vec([sp.sqrt(eps)*c for c in Ed],True) + Ipse*MV.vec([sp.sqrt(mu)*c for c in Hd],True)
res = nabla_left(Fd,X) - Ih*kk*Fd
pt = {x:sp.Rational(3,10), y:sp.Rational(-1,5), z:sp.Rational(7,10), w:sp.Rational(13,10), eps:sp.Rational(9,4), mu:1}
print("(B2) Dipolfeld: max|nabla F - ikF| am Testpunkt =", max(abs(complex(sp.N(c.subs(pt)))) for c in res.c))

# --- (C) Fundamentalloesung E_k = -(nabla + ik) Phi_k ---
Phk = sp.exp(Ih*k*r)/(4*sp.pi*r)
Ek = -(nabla_left(S(0,Phk),X) + Ih*k*S(0,Phk))
xh = MV.vec([x/r,y/r,z/r],True)
Ek_expl = S(0,Phk)*((S(0,1/r - Ih*k))*xh - S(0,Ih*k))
assert (Ek - Ek_expl).simplify().is_zero()
L = nabla_left(Ek,X) - Ih*k*Ek; R = nabla_right(Ek,X) - Ih*k*Ek
assert L.simplify().is_zero() and R.simplify().is_zero()
print("(C) E_k explizit, (nabla-ik)E_k = E_k(nabla-ik) = 0 fuer x!=0 : OK")
E0 = Ek.subs(k,0) if False else MV([sp.simplify(c.subs(k,0)) for c in Ek.c],True)
assert (E0 - MV.vec([x/(4*sp.pi*r**3),y/(4*sp.pi*r**3),z/(4*sp.pi*r**3)],True)).simplify().is_zero()
print("    statischer Grenzfall E_0 = x/(4 pi |x|^3) : OK")

# --- (D) Strahlungsbedingung (1 - xh)F = o(1/r) fuer Dipolfeld (numerisch in R) ---
import numpy as np
par = {w:sp.Rational(13,10), eps:sp.Rational(9,4), mu:1}
fF = [sp.lambdify(X, sp.sympify(c).subs(par), 'numpy') for c in Fd.c]
def Fnum(p): return MV([complex(f(*p)) if callable(f) else 0 for f in fF])
d = np.array([2,1,2])/3.0
for Rv in [10,100,1000,10000]:
    p = Rv*d; Fv = Fnum(p); G = Fv - MV.vec(d)*Fv
    print(f"(D) R={Rv:6d}:  R*|F| = {Rv*np.sqrt(Fv.norm2().real):.4e}   R*|(1-xh)F| = {Rv*np.sqrt(G.norm2().real):.4e}")
