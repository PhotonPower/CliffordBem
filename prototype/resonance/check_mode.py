import numpy as np
import pml2d_galerkin as pg
from mode2d import mode, rotor_apply
from enrich2d import exponent
k = -2.8 + 0.3j; a = exponent(k); c1, c2, _ = mode(k, a)
G, Mass, g, P = pg.build(30, k, 2.0, q=0.7)          # PML-Loesung (theta = 2), Streckung ab r0 = 0.25
Ns = len(Mass); hin = np.zeros((Ns, 8), complex); hin[:, 1] = 1
h = np.linalg.solve(G, np.repeat(Mass, 8)*hin.ravel()).reshape(Ns, 8)
m = 1; u = g.U[m]; phi = np.arctan2(u[1], u[0])       # Ecke V_1: Strahl 1 = Seite 1, Strahl 2 = Seite 0
def cosang(x, y): return abs(np.vdot(x, y))/np.linalg.norm(x)/np.linalg.norm(y)
for name, (cA, cB) in {"Strahl1=Seite m, Strahl2=Seite m-1, phi": (rotor_apply(phi, c1), rotor_apply(phi, c2)),
                       "vertauscht": (rotor_apply(phi, c2), rotor_apply(phi, c1)),
                       "phi + pi/2": (rotor_apply(phi + np.pi/2, c1), rotor_apply(phi + np.pi/2, c2)),
                       "phi - pi/2": (rotor_apply(phi - np.pi/2, c1), rotor_apply(phi - np.pi/2, c2))}.items():
    # Elemente auf Seite 1 nahe V_1 (s klein) und Seite 0 nahe V_1 (s nahe 2), Eckabstand 0.02 .. 0.2
    selA = [j for j in range(Ns) if g.SD[j] == 1 and 0.02 < g.Sm[j] < 0.2]
    selB = [j for j in range(Ns) if g.SD[j] == 0 and 0.02 < 2 - g.Sm[j] < 0.2]
    cA_fit = np.mean([cosang(h[j] - hin[j], cA) for j in selA]); cB_fit = np.mean([cosang(h[j] - hin[j], cB) for j in selB])
    print(f"{name:42s}: |cos| Seite A {cA_fit:.3f}, Seite B {cB_fit:.3f}")
