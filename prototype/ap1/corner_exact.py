"""Exakter Nachweis der TM-Blockformel fuer allgemeines alpha: Umschreiben in rationale
Funktion von x = e^{i pi a}, y = e^{i alpha}, z = e^{i a alpha} und Kuerzen."""
import sympy as sp, time, sys
from corner_symbolic import Lmat, vec, Icf, a, al, ghat
import corner_symbolic as cs
from corner_symbolic_general import J_TM
s = sp.Symbol('s', nonzero=True); x, y, z = sp.symbols('x y z', nonzero=True)
def rat(expr):
    e = expr.rewrite(sp.exp)
    e = sp.expand(e)
    def repl(ex):
        arg = sp.expand(ex.args[0]/sp.I)          # reelles Argument
        p = sp.Poly(arg, a, al)
        res = 1
        for (da, dal), c in p.terms():
            if (da, dal) == (1, 0): res *= x**sp.nsimplify(c/sp.pi)
            elif (da, dal) == (0, 1): res *= y**sp.nsimplify(c)
            elif (da, dal) == (1, 1): res *= z**sp.nsimplify(c)
            elif (da, dal) == (0, 0): res *= sp.exp(sp.I*c)
            else: raise ValueError((da, dal))
        return res
    return e.replace(lambda q: q.func == sp.exp, repl)
t0 = time.time()
Mf = sp.zeros(16, 16)
for i in range(2):
    for j in range(2):
        if i == j: B = -sp.cot(sp.pi*a)*Lmat(vec(cs.d[i]))*Lmat(vec(cs.n[j]))
        else:
            V = (cs.d[i]*Icf(a) - cs.d[j]*Icf(a + 1))/(2*sp.pi); B = -2*Lmat(vec(V))*Lmat(vec(cs.n[j]))
        Mf[8*i:8*i+8, 8*j:8*j+8] = B
ii = [1, 2, 9, 10]; M4 = Mf.extract(ii, ii).applyfunc(rat)
J4 = J_TM(s).applyfunc(rat)
T4 = sp.Rational(1, 2)*(sp.eye(4) + M4) + sp.Rational(1, 2)*(sp.eye(4) - M4)*J4
D = sp.together(T4.det(method='berkowitz'))
G = rat(-ghat(s**2, a)/(4*s**2))
diff = sp.cancel(sp.together(D - G))
print("TM-Block, allgemeines alpha: det - Vorhersage =", diff, " (%.0fs)" % (time.time() - t0))

def exact_block(ids, J4, pref_c_shift):
    pref, c, shift = pref_c_shift
    M4 = Mf.extract([ids[0], ids[1], 8+ids[0], 8+ids[1]], [ids[0], ids[1], 8+ids[0], 8+ids[1]]).applyfunc(rat)
    T4 = sp.Rational(1, 2)*(sp.eye(4) + M4) + sp.Rational(1, 2)*(sp.eye(4) - M4)*J4.applyfunc(rat)
    D = sp.together(T4.det(method='berkowitz')); G = rat(pref*ghat(c, a + shift))
    return sp.cancel(sp.together(D - G))
if len(sys.argv) > 1:
    t = sp.Symbol('t', nonzero=True)
    which = sys.argv[1]
    if which == 'aux1': print('Hilfs {1,e12}:', exact_block([0, 3], sp.diag(t, s, t, s), (-s**2/4, t/s, -1)))
    if which == 'aux2': print('Hilfs {e3,e123}:', exact_block([4, 7], sp.diag(s, t, s, t), (-s**2/4, t/s, -1)))
    if which == 'te':
        P = sp.Matrix([[0, -1], [1, 0]]); Q = sp.Matrix([[sp.cos(al), -sp.sin(al)], [sp.sin(al), sp.cos(al)]])
        J4 = sp.zeros(4, 4); J4[:2, :2] = P.inv()*sp.diag(s, 1/s)*P; J4[2:, 2:] = P.inv()*Q*sp.diag(s, 1/s)*Q.T*P
        print('TE {e13,e23}:', exact_block([5, 6], J4, (-1/(4*s**2), s**2, 0)))
