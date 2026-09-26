"""Absorber im Material statt in der Geometrie: Kontrast je Element kappa_j = kappa + i delta(r_j), mit einem
glatten Anstieg in log r unterhalb von r_abs:  delta = delta_max * S(ln(r_abs/r)/L),  S(x) = 3x^2 - 2x^3 (0..1).
Quasistatisch geht kappa nur ueber J ein, daher kein zusaetzlicher Rand. Referenz: Streckung (Galerkin)."""
import numpy as np, sys
from cl3 import MV
import pml2d_galerkin as pg
def panel_r(g):
    r = np.minimum(g.Sm, 2 - g.Sm); corner = (g.S0 == 0.0) | (g.S1 == 2.0)
    r[corner] = 0.5*(g.S1[corner] - g.S0[corner]); return r
def profile(r, kappa, r_abs, L, dmax):
    x = np.clip(np.log(r_abs/r)/L, 0, 1); return kappa + 1j*dmax*(3*x**2 - 2*x**3)
def solve(parts, kappas, x=np.array([5.0, 3.0])):
    Eg, Mass, g, P = parts; Ns = len(Mass)
    G = pg.assemble(Eg, Mass, g, kappas)
    hin = np.zeros((Ns, 8), complex); hin[:, 1] = 1.0
    hs = np.linalg.solve(G, np.repeat(Mass, 8)*hin.ravel()).reshape(Ns, 8) - hin
    F = np.zeros(8, complex)
    for j in range(Ns):
        Y, wd, _ = P[j]; Z = x[None, :] - Y; zz = np.einsum('gd,gd->g', Z, Z)
        K = np.einsum('g,gd->d', wd, Z/zz[:, None])/(2*np.pi); n = g.Nn[g.SD[j]]
        F += (MV.vec([K[0], K[1], 0])*MV.vec([n[0], n[1], 0])*MV(hs[j])).c
    return F[1]
if __name__ == '__main__':
    kappa = complex(sys.argv[1]); q = float(sys.argv[2]); depth = float(sys.argv[3])
    ng = int(np.ceil(np.log(depth/0.5)/np.log(q)))
    parts = pg.build(ng, kappa, 0.0, q=q, parts_only=True); g = parts[2]; r = panel_r(g)
    print(f"kappa = {kappa}, q = {q}, Tiefe {depth:.0e} (ng = {ng}, {len(r)} Elemente)")
    print(f"  ohne Absorber: F1 = {solve(parts, [kappa]*len(r)):.7f}")
    for r_abs, dec, dmax in [(0.1, 2, 1.5), (0.1, 3, 1.5), (0.01, 2, 1.5), (0.01, 3, 1.5), (0.01, 3, 3.0), (1e-3, 3, 1.5)]:
        L = dec*np.log(10)
        F = solve(parts, profile(r, kappa, r_abs, L, dmax))
        print(f"  r_abs {r_abs:.0e}, Anstieg {dec} Dekaden, delta_max {dmax}: F1 = {F:.7f}", flush=True)
