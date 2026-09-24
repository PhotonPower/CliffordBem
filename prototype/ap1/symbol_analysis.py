"""Symbolanalyse (Tangentialebene, n = e3) fuer das Transmissionsproblem."""
import sympy as sp, numpy as np
from cl3 import MV, SIGN, GRADE

def left_mat_sym(A):
    M = sp.zeros(8,8)
    for b in range(8):
        col = (A*MV.blade(b,1,True)).c
        for r in range(8): M[r,b] = col[r]
    return M

re_, rm, a, b = sp.symbols('r_epsilon r_mu a b', nonzero=True)
xi = MV.vec([1,0,0],True); e3 = MV.vec([0,0,1],True)
sigma = sp.I*left_mat_sym(xi*e3)
assert sp.simplify(sigma*sigma - sp.eye(8)) == sp.zeros(8,8)
pp, pm = (sp.eye(8)+sigma)/2, (sp.eye(8)-sigma)/2

# Transmissionsabbildung J (n = e3): Blades 1=e1,2=e2,4=e3,3=e12,5=e13,6=e23
# Bivektor I*v: I e3 = e12 (normal), I e1 = e23, I e2 = -e13 (tangential)
Jd = {0:a, 1:sp.sqrt(re_), 2:sp.sqrt(re_), 4:1/sp.sqrt(re_),
      3:1/sp.sqrt(rm), 5:sp.sqrt(rm), 6:sp.sqrt(rm), 7:b}
J = sp.diag(*[Jd[i] for i in range(8)])

def colbasis(P):
    return sp.Matrix.hstack(*P.columnspace())
Bp, Bm = colbasis(pp), colbasis(pm)
print("rank p+ =", Bp.shape[1], " rank p- =", Bm.shape[1])

D = sp.factor(sp.simplify((sp.Matrix.hstack(J.inv()*Bp, Bm)).det()))
print("\n[1] Intrinsische Bedingung (Zerlegung C^8 = J^-1 R(p+) (+) R(p-)):")
print("    det ~", D)

T = pp + J.inv()*pm*J
DT = sp.factor(sp.simplify(T.det()))
print("\n[2] Naive Gleichung T = E2^+ + J^-1 E1^- J,  det sigma(T) ~", DT)

# Welche Bloecke koppeln? Grade-Struktur von sigma
print("\n[3] sigma koppelt Blades:")
for c in range(8):
    rows = [r for r in range(8) if sigma[r,c]!=0]
    print("   ", ["1","e1","e2","e12","e3","e13","e23","e123"][c], "->", [["1","e1","e2","e12","e3","e13","e23","e123"][r] for r in rows])

print("\n[4] Vereinfachte Determinanten:")
s_e, s_m = sp.symbols('s_epsilon s_mu', positive=True)   # s = sqrt(r)
sub = {re_: s_e**2, rm: s_m**2}
print("    det[1] =", sp.factor(sp.simplify(D.subs(sub))))
print("    det[2] =", sp.factor(sp.simplify(DT.subs(sub))))

print("\n[5] Eigenwerte von sigma(T) mit a = sqrt(r_mu), b = sqrt(r_epsilon):")
Tn = T.subs({a: sp.sqrt(rm), b: sp.sqrt(re_)})
ev = Tn.eigenvals()
for e_, m_ in ev.items(): print("    ", sp.simplify(e_), " (Vielfachheit", m_, ")")
print("    sigma(T) - Id =", "Null" if sp.simplify(Tn - sp.eye(8)) == sp.zeros(8,8) else "ungleich Null -> T nicht Id+kompakt")

D2 = sp.factor(sp.simplify(sp.Matrix.hstack(Bp, J.inv()*Bm).det().subs(sub)))
print("\n[6] Duale Bedingung R(p+) cap J^-1 R(p-) = {0}:  det ~", D2)

# Numerik: Eigenwerte fuer typische Kontraste
print("\n[7] Symbol-Eigenwerte (nichttriviale) fuer typische Kontraste r_eps (r_mu = 1):")
for r in [2.25, 12.0, -1.0+0.1j, -11.0+1.2j, -0.5+0.05j]:
    s = np.sqrt(complex(r)); d = 0.5j*(s - 1/s)
    print(f"    r_eps = {r!s:>12}:  lambda = {1+d:.3f}, {1-d:.3f}")

print("\n[8] Modifizierte Gleichung T_W = E2^+ + W E1^- J  mit W = 1:")
TW = pp + pm*J
DW = sp.factor(sp.simplify(TW.det().subs(sub)))
print("    det sigma(T_1) ~", DW)
DWd = sp.factor(sp.simplify(sp.Matrix.hstack(Bp, Bm).det()))
print("    duale Bedingung ran p+ cap ran p- (W=1): det ~", DWd)

print("\n[9] Eigenwerte von sigma(T_1) (Standardwahl a, b):")
T1 = (pp + pm*J).subs({a: sp.sqrt(rm), b: sp.sqrt(re_)})
for e_, m_ in T1.eigenvals().items(): print("    ", sp.simplify(e_), "(x%d)" % m_)
for r in [2.25, 12.0, -11+1.2j, -1+0.1j]:
    M = np.array(T1.subs({re_: r, rm: 1}).evalf(), dtype=complex)
    print(f"    r_eps={r}: Eigenwerte", np.round(np.sort_complex(np.linalg.eigvals(M)), 3))

print("\n[10] Eigenwerte von sigma(T_1) fuer allgemeine a, b:")
T1g = pp + pm*J
for e_, m_ in T1g.eigenvals().items(): print("    ", sp.simplify(e_), "(x%d)" % m_)
print("    Energiewahl a = conj(r_eps) sqrt(r_mu), b = conj(r_mu) sqrt(r_eps):")
for r in [2.25, 12.0, -11+1.2j, -2+0.3j]:
    aE, bE = np.conj(r)*1.0, 1.0*np.sqrt(r)
    M = np.array(T1g.subs({re_: r, rm: 1, a: aE, b: bE}).evalf(), dtype=complex)
    print(f"    r_eps={r}: ", np.round(np.sort_complex(np.linalg.eigvals(M)), 3))

print("\n[11] Eigenwerte von 2(1+J)^{-1} sigma(T_1):")
P1 = 2*(sp.eye(8)+J).inv()*T1g
for e_, m_ in P1.eigenvals().items(): print("    ", sp.simplify(e_), "(x%d)" % m_)
for r in [2.25, 12.0, -11+1.2j, -2+0.3j]:
    for lab, aa, bb_ in [('Standard', 1.0, np.sqrt(r)), ('Energie', np.conj(r), np.sqrt(r))]:
        M = np.array(P1.subs({re_: r, rm: 1, a: aa, b: bb_}).evalf(), dtype=complex)
        print(f"    r_eps={r}, {lab}: ", np.round(np.sort_complex(np.linalg.eigvals(M)), 3))
