"""Identitaeten fuer den chiralen Eindeutigkeitsbeweis (SymPy):
Hilfsanteil G_aux = sum_pm (1/(i k_pm)) (nabla + i k_pm)(gamma_pm P_pm),  k_pm = w(sqrt(eps mu) -+ chi).
Behauptung: D_aux = eps E_aux + i chi H_aux = grad(phi)/(i w),  B_aux = mu H_aux - i chi E_aux = grad(psi)/(i w),
mit alpha = (g+ + g-)/2, beta = i(g+ - g-)/2, phi = alpha/sqrt(mu), psi = beta/sqrt(eps).
Ausserdem: |mu| d(phi) conj(phi) + |eps| d(psi) conj(psi) = (d(g+) conj(g+) + d(g-) conj(g-))/2."""
import sympy as sp
from cl3 import MV, nabla_left
x, y, z = X = sp.symbols('x y z', real=True)
w, se, sm, chi = sp.symbols('omega s_e s_m chi')       # se = sqrt(eps), sm = sqrt(mu) (komplex zugelassen)
eps, mu = se**2, sm**2
gp_, gm_ = sp.Function('gp')(*X), sp.Function('gm')(*X)
Ip = MV.blade(7, 1, True); one = MV.blade(0, 1, True)
Pp = (one + Ip*sp.I)*sp.Rational(1, 2); Pm = (one - Ip*sp.I)*sp.Rational(1, 2)
G = MV(symbolic=True)
for g, P, s in [(gp_, Pp, 1), (gm_, Pm, -1)]:
    ks = w*(se*sm - s*chi)
    U = P*MV.blade(0, g, True)
    G = G + nabla_left(U, X)*(1/(sp.I*ks)) + U
E = sp.Matrix([G.c[1], G.c[2], G.c[4]])/se
Hm = (Ip*G.grade(2))*(-1); H = sp.Matrix([Hm.c[1], Hm.c[2], Hm.c[4]])/sm
alpha = sp.simplify(G.c[0]); beta = sp.simplify(G.c[7])
print("alpha = (g+ + g-)/2:", sp.simplify(alpha - (gp_ + gm_)/2) == 0, "  beta = i(g+ - g-)/2:", sp.simplify(beta - sp.I*(gp_ - gm_)/2) == 0)
phi, psi = alpha/sm, beta/se
grad = lambda f: sp.Matrix([sp.diff(f, t) for t in X])
D = eps*E + sp.I*chi*H; B = mu*H - sp.I*chi*E
print("D_aux = grad(phi)/(i w):", sp.simplify(D - grad(phi)/(sp.I*w)) == sp.zeros(3, 1))
print("B_aux = grad(psi)/(i w):", sp.simplify(B - grad(psi)/(sp.I*w)) == sp.zeros(3, 1))
# Gewichtsidentitaet (punktweise, mit unabhaengigen Symbolen fuer Werte und Ableitungen)
a1, a2, d1, d2 = sp.symbols('a1 a2 d1 d2')              # Werte g+, g- und Normalableitungen
sE, sM = sp.symbols('sE sM')                             # sqrt(eps), sqrt(mu)
ph, dph = (a1 + a2)/(2*sM), (d1 + d2)/(2*sM); ps, dps = sp.I*(a1 - a2)/(2*sE), sp.I*(d1 - d2)/(2*sE)
W = sM*sp.conjugate(sM)*dph*sp.conjugate(ph) + sE*sp.conjugate(sE)*dps*sp.conjugate(ps)
print("Gewichtsidentitaet:", sp.simplify(sp.expand(W - (d1*sp.conjugate(a1) + d2*sp.conjugate(a2))/2)) == 0)
