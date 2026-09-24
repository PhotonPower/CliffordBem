"""Chirale Eckformel: gekoppelter TM+TE-Block (8x8) mit exakten Mellin-Eintraegen.
Bestimmt den Vorfaktor K in det = K (A - B S^2 + C S^4) und prueft die Formel mit 40 Stellen."""
import sympy as sp, mpmath as mp
from corner_symbolic import Lmat, vec, Icf, a, al, S
import corner_symbolic as cs
s1, s2, ch = sp.symbols('s1 s2 chi')
def Jray(d, n):
    """J auf (e1, e2, e13, e23) fuer Strahl mit Tangente d, Normale n (in-plane)."""
    C1 = sp.Matrix([[s1**2, sp.I*ch], [-sp.I*ch, s2**2]]); Nn = sp.diag(s1, s2)*C1.inv()
    # Koordinaten: E = (c1, c2); H-Vektor b = (c_e23, -c_e13)
    M = sp.zeros(4, 4)
    for k in range(4):
        v = sp.zeros(4, 1); v[k] = 1
        E = sp.Matrix([v[0], v[1]]); b = sp.Matrix([v[3], -v[2]])
        dd = sp.Matrix([d[0], d[1]]); nn = sp.Matrix([n[0], n[1]])
        Et, En = (E.T*dd)[0], (E.T*nn)[0]; bt, bn = (b.T*dd)[0], (b.T*nn)[0]
        En2, bn2 = list(Nn*sp.Matrix([En, bn]))
        E2 = s1*Et*dd + En2*nn; b2 = s2*bt*dd + bn2*nn
        M[:, k] = sp.Matrix([E2[0], E2[1], -b2[1], b2[0]])
    return M
Mf = sp.zeros(16, 16)
for i in range(2):
    for j in range(2):
        if i == j: B = -sp.cot(sp.pi*a)*Lmat(vec(cs.d[i]))*Lmat(vec(cs.n[j]))
        else:
            V = (cs.d[i]*Icf(a) - cs.d[j]*Icf(a + 1))/(2*sp.pi); B = -2*Lmat(vec(V))*Lmat(vec(cs.n[j]))
        Mf[8*i:8*i+8, 8*j:8*j+8] = B
ids = [1, 2, 5, 6]; ii = ids + [8 + k for k in ids]
M8 = Mf.extract(ii, ii)
J8 = sp.zeros(8, 8); J8[:4, :4] = Jray(cs.d[0], cs.n[0]); J8[4:, 4:] = Jray(cs.d[1], cs.n[1])
T8 = sp.Rational(1, 2)*(sp.eye(8) + M8) + sp.Rational(1, 2)*(sp.eye(8) - M8)*J8
e, u = s1**2, s2**2
A = ((1 + e)*(1 + u) - ch**2)**2; C = ((1 - e)*(1 - u) - ch**2)**2
Bq = 2*(1 + e**2)*(1 + u**2) - 8*e*u + 4*ch**2*(3 - e*u) + 2*ch**4
Pq = A - Bq*S(a)**2 + C*S(a)**4
fD = sp.lambdify((a, al, s1, s2, ch), T8, 'mpmath'); fP = sp.lambdify((a, al, s1, s2, ch), Pq, 'mpmath')
mp.mp.dps = 40
def K(av, alv, sv1, sv2, cv):
    Tn = mp.matrix(fD(av, alv, sv1, sv2, cv)); return mp.det(Tn)/fP(av, alv, sv1, sv2, cv)
media = [(mp.sqrt(mp.mpc(-3, 0.3)), mp.mpf(1), mp.mpf('0.4')), (mp.sqrt(mp.mpc(2.25, 0.1)), mp.sqrt(mp.mpf(1.3)), mp.mpc(0.25, 0.02))]
for sv1, sv2, cv in media:
    ks = [K(mp.mpc(av), mp.mpf(alv), sv1, sv2, cv) for av, alv in [(0.31+0.7j, 1.1), (0.62-1.2j, 2.3), (0.45+0.1j, 0.7)]]
    e_, u_ = sv1**2, sv2**2
    print('K-Schwankung ueber (a, alpha):', mp.nstr(max(abs(k/ks[0]-1) for k in ks), 5),
          '  K*16*(e u - chi^2)^2/(e u) =', mp.nstr(ks[0]*16*(e_*u_ - cv**2)**2/(e_*u_), 12),
          '  K*16*(e u - chi^2) =', mp.nstr(ks[0]*16*(e_*u_ - cv**2), 12))
