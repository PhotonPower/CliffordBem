"""Zerlegung G = G_M + (1/(ik)) (nabla + ik)(alpha + I beta) einer allgemeinen Dirac-Loesung:
G_M hat keine Grade 0,3 und ist Dirac-Loesung (also Maxwell-Feld)."""
import numpy as np
from fields import gp, dirac_kernel, blade, dirac_residual
rng = np.random.default_rng(11); c = rng.normal(size=8) + 1j*rng.normal(size=8)
xs = np.array([0.1, -0.2, 0.05]); x = np.array([[0.9, 0.4, -0.6], [-0.7, 1.1, 0.3]])
for k in [1.3, 0.9+0.4j]:
    G = lambda y: gp(dirac_kernel(y, xs, k), c)
    def U(y):                                   # u = alpha + I beta
        g = G(y); u = np.zeros_like(g); u[..., 0] = g[..., 0]; u[..., 7] = g[..., 7]; return u
    def Gaux(y, h=1e-5):
        r = 1j*k*U(y)
        for j, b in enumerate([1, 2, 4]):
            d = np.zeros(3); d[j] = h; r = r + gp(blade(b), (U(y+d) - U(y-d))/(2*h))
        return r/(1j*k)
    GM = lambda y: G(y) - Gaux(y)
    gm = GM(x)
    print(f"k={k}: |<G_M>_0,3|/|G| = {np.abs(gm[:, [0, 7]]).max()/np.abs(G(x)).max():.1e},  "
          f"|(nabla-ik)G_M|/|G_M| = {np.abs(dirac_residual(GM, x, k, h=1e-3)).max()/np.abs(gm).max():.1e}")
