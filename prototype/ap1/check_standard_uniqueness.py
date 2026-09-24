"""Symbolische Pruefung der Identitaeten im Eindeutigkeitsbeweis fuer die Standardwahl:
Hilfsteil G_aux = (1/(ik))(nabla + ik)(alpha + I beta) liefert
  E_aux = grad(alpha)/(i k sqrt(eps)),  H_aux = grad(beta)/(i k sqrt(mu)),
und mit phi = alpha/sqrt(mu), psi = beta/sqrt(eps):
  E_aux = grad(phi)/(i omega eps),  H_aux = grad(psi)/(i omega mu)."""
import sympy as sp
from cl3 import MV, nabla_left
x, y, z = X = sp.symbols('x y z', real=True)
w, eps, mu = sp.symbols('omega epsilon mu', positive=True)
k = w*sp.sqrt(eps)*sp.sqrt(mu)
al = sp.Function('alpha')(*X); be = sp.Function('beta')(*X)
Ip = MV.blade(7, 1, True)
u = MV.blade(0, al, True) + Ip*MV.blade(0, be, True)
G = nabla_left(u, X)*(1/(sp.I*k)) + u
E = [G.c[1]/sp.sqrt(eps), G.c[2]/sp.sqrt(eps), G.c[4]/sp.sqrt(eps)]
Hm = (Ip*G.grade(2))*(-1)
H = [Hm.c[1]/sp.sqrt(mu), Hm.c[2]/sp.sqrt(mu), Hm.c[4]/sp.sqrt(mu)]
gp = lambda f: [sp.diff(f, t) for t in X]
phi = al/sp.sqrt(mu); psi = be/sp.sqrt(eps)
okE = all(sp.simplify(e - g/(sp.I*w*eps)) == 0 for e, g in zip(E, gp(phi)))
okH = all(sp.simplify(h - g/(sp.I*w*mu)) == 0 for h, g in zip(H, gp(psi)))
print("E_aux = grad(phi)/(i omega eps):", okE)
print("H_aux = grad(psi)/(i omega mu):", okH)
print("Grade 0/3 von G_aux gleich alpha, I beta:", sp.simplify(G.c[0] - al) == 0, sp.simplify(G.c[7] - be) == 0)
# Oberflaechenidentitaet: n.curl(grad f) = 0 (fuer beliebiges f) -- trivial, hier als Vollstaendigkeitscheck
f = sp.Function('f')(*X); g_ = gp(f)
curl = [sp.diff(g_[2], y) - sp.diff(g_[1], z), sp.diff(g_[0], z) - sp.diff(g_[2], x), sp.diff(g_[1], x) - sp.diff(g_[0], y)]
print("curl grad f = 0:", all(sp.simplify(c) == 0 for c in curl))
