"""Symbolischer Beweis der Blockfaktorisierung fuer allgemeinen Winkel alpha (TM-Block;
TE analog durch Dualitaet, Hilfsbloecke mit allgemeinen Parametern).
Strategie: exakte Mellin-Bloecke mit geschlossenen Kreuzintegralen, J je Strahl in globalen
Koordinaten, Determinante in Exponentialform x = e^{i pi a}, y = e^{i alpha} als rationale Funktion."""
import sympy as sp
from corner_symbolic import block_symbol, Lmat, vec, Icf, a, al, S, ghat
s, t, r = sp.symbols('s t r', nonzero=True)
x, y = sp.symbols('x y', nonzero=True)
def to_xy(expr):
    e = sp.expand_trig(expr).rewrite(sp.exp)
    e = e.subs({sp.exp(sp.I*sp.pi*a): x}).subs({sp.exp(-sp.I*sp.pi*a): 1/x})
    e = sp.powsimp(sp.expand(e))
    e = e.replace(lambda z: z.func == sp.exp, lambda z: sp.exp(sp.expand(z.args[0])))
    return e
def J_TM(s):
    # Strahl 1: d1 = e1, n1 = -e2 -> Komponenten (e1, e2): diag(s, 1/s)
    # Strahl 2: d2 = (cos, sin), n2 = (-sin, cos):  Q diag(s, 1/s) Q^T
    Q = sp.Matrix([[sp.cos(al), -sp.sin(al)], [sp.sin(al), sp.cos(al)]])
    J2 = Q*sp.diag(s, 1/s)*Q.T
    J = sp.zeros(4, 4); J[0, 0] = s; J[1, 1] = 1/s; J[2:, 2:] = J2
    return J
def check(name, ids, J4, c, shift, pref):
    from corner_symbolic import d, n
    import corner_symbolic as cs
    # Mellin-Block ohne J aufbauen (block_symbol erwartet Diagonale -> hier allgemeine J4)
    Mf = sp.zeros(16, 16)
    for i in range(2):
        for j in range(2):
            if i == j: B = -sp.cot(sp.pi*a)*Lmat(vec(cs.d[i]))*Lmat(vec(cs.n[j]))
            else:
                V = (cs.d[i]*Icf(a) - cs.d[j]*Icf(a + 1))/(2*sp.pi); B = -2*Lmat(vec(V))*Lmat(vec(cs.n[j]))
            Mf[8*i:8*i+8, 8*j:8*j+8] = B
    ii = [ids[0], ids[1], 8 + ids[0], 8 + ids[1]]; M4 = Mf.extract(ii, ii)
    T4 = sp.Rational(1, 2)*(sp.eye(4) + M4) + sp.Rational(1, 2)*(sp.eye(4) - M4)*J4
    D = T4.det(method='berkowitz'); G = pref*ghat(c, a + shift)
    diff = D - G
    # numerische Pruefung mit hoher Genauigkeit an zufaelligen Punkten
    import mpmath as mp; mp.mp.dps = 40
    f = sp.lambdify((a, al, s, t), diff, 'mpmath'); g = sp.lambdify((a, al, s, t), G, 'mpmath')
    worst = 0
    for av, alv, sv, tv in [(0.31+0.7j, 1.1, 1.7+0.4j, 0.6-0.3j), (0.62-1.2j, 2.3, -0.4+1.3j, 1.9+0.2j), (0.45+0.1j, 0.7, 2.2-0.5j, -0.7+0.8j)]:
        worst = max(worst, abs(f(mp.mpc(av), mp.mpf(alv), mp.mpc(sv), mp.mpc(tv)))/abs(g(mp.mpc(av), mp.mpf(alv), mp.mpc(sv), mp.mpc(tv))))
    print(f"{name:14s}: max |det - Vorhersage|/|Vorhersage| (40 Stellen) = {mp.nstr(worst, 5)}", flush=True)

if __name__ == '__main__':
    check('TM {e1,e2}', [1, 2], J_TM(s), s**2, 0, -1/(4*s**2))
    # Hilfsblock {1,e12}: J = diag(t, s) je Strahl (t = a_H, s = sqrt(r_mu)), Kontrast t/s, Vorfaktor -s^2/4
    check('Hilfs {1,e12}', [0, 3], sp.diag(t, s, t, s), t/s, -1, -s**2/4)
    check('Hilfs {e3,e123}', [4, 7], sp.diag(s, t, s, t), t/s, -1, -s**2/4)
    # TE-Block {e13,e23}: Komponenten (e13, e23) = (-H_2, H_1) bis auf Vorzeichen; J wie TM mit H.
    # In (e23, e13)-Koordinaten entspricht e23 <-> H_1, e13 <-> -H_2. Konjugation mit diag(1,-1) und Vertauschung:
    P = sp.Matrix([[0, -1], [1, 0]])        # (e13, e23) -> (H1, H2)
    Q = sp.Matrix([[sp.cos(al), -sp.sin(al)], [sp.sin(al), sp.cos(al)]])
    JH1 = sp.diag(s, 1/s); JH2 = Q*sp.diag(s, 1/s)*Q.T
    J4 = sp.zeros(4, 4); J4[:2, :2] = P.inv()*JH1*P; J4[2:, 2:] = P.inv()*JH2*P
    check('TE {e13,e23}', [5, 6], J4, s**2, 0, -1/(4*s**2))
