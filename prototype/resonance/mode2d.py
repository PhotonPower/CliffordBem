"""Transparentes Eckelement mit der exakten Eckmode: Auf den beiden Eckelementen einer Ecke ist die Dichte
gamma * r^{a-1} * c_Flaeche mit dem Nullvektor (c_1, c_2) des Mellin-Symbols von T_1 zum physikalischen
Exponenten a (ein Unbekannter je Ecke). Getestet mit c~_Flaeche * 1 auf denselben Elementen."""
import numpy as np
from cl3 import MV
from corner import M_E, T1_symbol
from enrich2d import exponent, build as build_full
from pml2d import V
def rotor_apply(phi, c):
    R = MV(np.array([np.cos(phi/2), 0, 0, -np.sin(phi/2), 0, 0, 0, 0], complex))   # e^{-e12 phi/2}
    return (R*MV(np.asarray(c, complex))*R.reverse()).c
def mode(kappa, a):
    M = M_E(np.pi/2, a.imag, a.real); T = T1_symbol(np.pi/2, a.imag, kappa, ME=M)
    U, s, Vh = np.linalg.svd(T); c = Vh[-1].conj()
    return c[:8], c[8:], s[-1]
def far_field(ng, kappa, q=0.7, test='conj', x=np.array([5.0, 3.0])):
    a = exponent(kappa); c1, c2, smin = mode(kappa, a)
    G, Lr, g, P, trial = build_full(ng, kappa, a, q); Ns = len(Lr)
    first = {k: np.nonzero((g.SD == k) & (g.S0 == 0.0))[0][0] for k in range(4)}
    last = {k: np.nonzero((g.SD == k) & (g.S1 == 2.0))[0][0] for k in range(4)}
    corner_panels = set(first.values()) | set(last.values())
    reg = [j for j in range(Ns) if j not in corner_panels]
    nu = 8*len(reg) + 4
    Pm = np.zeros((8*Ns, nu), complex); Rm = np.zeros((nu, 8*Ns), complex)
    for r_i, j in enumerate(reg):
        for b in range(8): Pm[8*j + b, 8*r_i + b] = 1; Rm[8*r_i + b, 8*j + b] = 1
    faces = {}
    for m in range(4):                                     # Ecke V_m: Strahl 1 = Seite m, Strahl 2 = Seite m-1
        u = g.U[m]; phi = np.arctan2(u[1], u[0])
        cA, cB = rotor_apply(phi, c1), rotor_apply(phi, c2)
        col = 8*len(reg) + m
        for j, cf in [(first[m], cA), (last[(m - 1) % 4], cB)]:
            Pm[8*j:8*j + 8, col] = cf
            Rm[col, 8*j:8*j + 8] = cf.conj() if test == 'conj' else cf
            faces[j] = cf
    hin = np.zeros((Ns, 8), complex); hin[:, 1] = 1.0
    b = np.repeat(Lr, 8)*hin.ravel()
    gam = np.linalg.solve(Rm@G@Pm, Rm@b); coef = (Pm@gam).reshape(Ns, 8)
    # Fernfeld (wie enrich2d.far_field)
    F = np.zeros(8, complex)
    for j in range(Ns):
        Y, w = P[j][0], P[j][1]; Z = x[None, :] - Y; zz = np.einsum('gd,gd->g', Z, Z); n = g.Nn[g.SD[j]]
        Kt = np.einsum('g,gd->d', trial(j, P[j]), Z/zz[:, None])/(2*np.pi); K1 = np.einsum('g,gd->d', w, Z/zz[:, None])/(2*np.pi)
        nv = MV.vec([n[0], n[1], 0])
        F += (MV.vec([Kt[0], Kt[1], 0])*nv*MV(coef[j])).c - (MV.vec([K1[0], K1[1], 0])*nv*MV(hin[j])).c
    return F, smin
if __name__ == '__main__':
    import sys
    assert np.allclose(rotor_apply(np.pi/2, np.eye(8)[1]), np.eye(8)[2]), "Rotor dreht e1 nicht auf e2"
    k = complex(sys.argv[1]); q = float(sys.argv[2]); depths = [float(v) for v in sys.argv[3].split(',')]; test = sys.argv[4] if len(sys.argv) > 4 else 'conj'
    for dep in depths:
        ng = int(np.ceil(np.log(dep/0.5)/np.log(q))); F, smin = far_field(ng, k, q, test)
        print(f"  kappa={k}, q={q}, Tiefe {dep:.0e} (ng={ng}), Test {test}: F1 = {F[1].real:+.7f}{F[1].imag:+.7f}i  (Symbol-sigma_min {smin:.1e})", flush=True)
