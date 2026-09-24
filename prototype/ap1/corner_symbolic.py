"""Symbolischer Nachweis der Blockfaktorisierung des Eck-Mellin-Symbols (SymPy).
Kreuzintegrale in geschlossener Form: int_0^oo u^{mu-1}/(u^2+2u cos th+1) du = pi sin((1-mu)th)/(sin th sin(pi mu)),
th = pi - alpha.  Bloecke je Strahl: {e1,e2} (TM), {e13,e23} (TE), {1,e12}, {e3,e123} (Hilfsgrade)."""
import sympy as sp
from cl3 import MV
a, al = sp.symbols('a alpha')
s, t = sp.symbols('s t', nonzero=True)      # s = sqrt(Kontrast) fuer J, t = Hilfsparameter
th = sp.pi - al
Icf = lambda mu: sp.pi*sp.sin((1 - mu)*th)/(sp.sin(th)*sp.sin(sp.pi*mu))
d = [sp.Matrix([1, 0, 0]), sp.Matrix([sp.cos(al), sp.sin(al), 0])]
n = [sp.Matrix([0, -1, 0]), sp.Matrix([-sp.sin(al), sp.cos(al), 0])]
def Lmat(A):
    M = sp.zeros(8, 8)
    for b in range(8):
        col = (A*MV.blade(b, 1, True)).c
        for r in range(8): M[r, b] = col[r]
    return M
vec = lambda v: MV.vec(list(v), True)
def block_symbol(ids, Jdiag):
    """4x4-Mellin-Symbol von T_1 auf dem Block ids (2 Blades) fuer beide Strahlen."""
    Mf = sp.zeros(16, 16)
    for i in range(2):
        for j in range(2):
            if i == j:
                B = -sp.cot(sp.pi*a)*Lmat(vec(d[i]))*Lmat(vec(n[j]))
            else:
                V = (d[i]*Icf(a) - d[j]*Icf(a + 1))/(2*sp.pi)
                B = -2*Lmat(vec(V))*Lmat(vec(n[j]))
            Mf[8*i:8*i+8, 8*j:8*j+8] = B
    ii = [ids[0], ids[1], 8 + ids[0], 8 + ids[1]]
    M4 = Mf.extract(ii, ii)
    J4 = sp.diag(*Jdiag)
    return sp.Rational(1, 2)*(sp.eye(4) + M4) + sp.Rational(1, 2)*(sp.eye(4) - M4)*J4
def S(x): return sp.sin((sp.pi - al)*x)/sp.sin(sp.pi*x)
def ghat(c, x): return (1 - c)**2*S(x)**2 - (1 + c)**2
if __name__ == '__main__':
    import random
    # J-Diagonalen je Block (Strahl 1 Normale e2-Richtung, Strahl 2 allgemein):
    # TM {e1,e2}: tangential sqrt(r) , normal 1/sqrt(r). Fuer Strahl 1 (d=e1, n=-e2): e1 tangential, e2 normal.
    # Fuer Strahl 2 ist die Blade-Basis nicht an (d,n) ausgerichtet -> J nicht diagonal; daher alpha = pi/2
    # symbolisch in a, und J mit allgemeinem s.
    res = {}
    for name, ids, J in [('TM {e1,e2}', [1, 2], lambda: [s, 1/s, 1/s, s]),
                         ('Hilfs {1,e12}', [0, 3], lambda: [t, 1, t, 1]),
                         ('Hilfs {e3,e123}', [4, 7], lambda: [1, t, 1, t]),
                         ('TE {e13,e23}', [5, 6], lambda: [1/s, s, s, 1/s])]:
        T4 = block_symbol(ids, J()).subs(al, sp.pi/2)
        D = sp.simplify(T4.det())
        c = {'TM {e1,e2}': s**2, 'TE {e13,e23}': s**2, 'Hilfs {1,e12}': t, 'Hilfs {e3,e123}': t}[name]
        shift = 0 if name.startswith('T') else -1
        G = ghat(c, a + shift).subs(al, sp.pi/2)
        ratio = sp.simplify(sp.expand_trig(sp.simplify(D/G)))
        res[name] = ratio
        print(f"{name:16s}: det / ghat = {ratio}", flush=True)
