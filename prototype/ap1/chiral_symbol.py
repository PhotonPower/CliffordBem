"""Hauptsymbol (glatter Rand) des Transmissionsproblems mit chiralem Innenmedium:
Ausnahmemenge in Abhaengigkeit von (eps1, mu1, chi1), Aussenmedium eps2 = mu2 = 1, chi2 = 0."""
import sympy as sp
from cl3 import MV
s1, s2, c = sp.symbols('s1 s2 chi', nonzero=True)
e, u = s1**2, s2**2
def Lm(A):
    M = sp.zeros(8, 8)
    for b in range(8):
        col = (A*MV.blade(b, 1, True)).c
        for r in range(8): M[r, b] = col[r]
    return M
sig = sp.I*Lm(MV.vec([1, 0, 0], True)*MV.vec([0, 0, 1], True))
pp, pm = (sp.eye(8) + sig)/2, (sp.eye(8) - sig)/2
se, sm = s1, s2
C1 = sp.Matrix([[e, sp.I*c], [-sp.I*c, u]]); C2 = sp.eye(2)
Nn = sp.diag(se, sm)*C1.inv()*C2
# J in Blade-Basis fuer n = e3: tangentiale E (e1,e2): se; tangentiale H (e23 = I e1, e13 = -I e2): sm;
# normale (e3, e12 = I e3) gemischt ueber Nn; Hilfsgrade: a = sm, b = se (Standardwahl)
J = sp.zeros(8, 8)
J[0, 0] = sm; J[7, 7] = se; J[1, 1] = se; J[2, 2] = se; J[6, 6] = sm; J[5, 5] = sm
J[4, 4] = Nn[0, 0]; J[4, 3] = Nn[0, 1]; J[3, 4] = Nn[1, 0]; J[3, 3] = Nn[1, 1]
def cols(P): return sp.Matrix.hstack(*P.columnspace())
Bp, Bm = cols(pp), cols(pm)
D = sp.factor(sp.cancel(sp.Matrix.hstack(J.inv()*Bp, Bm).det(method='berkowitz')))
print("Zerlegungsbedingung: det ~", D)
T1 = pp + pm*J
DT = sp.factor(sp.cancel(T1.det(method='berkowitz')))
print("det sigma(T_1) ~", DT)
# Kontrolle chi = 0
print("chi=0:", sp.factor(D.subs(c, 0)))
