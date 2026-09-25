"""3D-Uebertragung der Streckung: Nullstellen von z.z auf dem komplex gestreckten Rand.
Kante entlang e3, Flaechen A (x2 = 0, x1 = r >= 0) und B (Winkel alpha). Gestreckt wird der Kantenabstand r,
die Koordinate t entlang der Kante bleibt reell:  z.z = (t1 - t2)^2 + w_T,  w_T = transversaler Anteil.
z.z = 0 fuer ein reelles t1 - t2 genau dann, wenn w_T auf der negativen reellen Achse liegt."""
import numpy as np
def psi_spiral(r, th, r0): return np.where(r < r0, th*np.log(r/r0), 0.0)
def psi_capped(r, th, r0, cap): return np.maximum(psi_spiral(r, th, r0), -cap)
def rt(r, psi): return r*np.exp(1j*psi)
def worst(psi_fn, alpha, n=1500, rmin=1e-6, r0=0.25):
    r = np.geomspace(rmin, 1.0, n); R1, R2 = np.meshgrid(r, r, indexing='ij')
    a, b = rt(R1, psi_fn(R1)), rt(R2, psi_fn(R2))
    wT = a*a + b*b - 2*a*b*np.cos(alpha)                  # alpha = 0: dieselbe Flaeche
    dist2 = R1**2 + R2**2 - 2*R1*R2*np.cos(alpha)         # reeller Abstand zum Vergleich
    off = dist2 > 1e-12*np.maximum(R1, R2)**2
    ang = np.abs(np.angle(wT[off]))                        # pi = negative reelle Achse
    return ang.max(), np.count_nonzero(ang > np.pi - 1e-3)
for name, fn in [("Spirale theta = 2", lambda r: psi_spiral(r, 2.0, 0.25)),
                 ("Spirale theta = 0,5", lambda r: psi_spiral(r, 0.5, 0.25)),
                 ("begrenzt |psi| <= 1,5 (< pi/2)", lambda r: psi_capped(r, 2.0, 0.25, 1.5)),
                 ("begrenzt |psi| <= 1,6 (> pi/2)", lambda r: psi_capped(r, 2.0, 0.25, 1.6))]:
    for alpha, an in [(0.0, "gleiche Flaeche"), (np.pi/2, "Nachbarflaeche 90 Grad")]:
        mx, cnt = worst(fn, alpha)
        print(f"{name:32s} {an:24s}: max |arg w_T| = {mx:.4f} (pi = {np.pi:.4f}), Paare mit z.z = 0 fuer ein reelles dt: {cnt}")
# Ecke (Kegelspitze): isotrope Streckung aller Koordinaten um die Spitze, rho -> rho e^{i psi(rho)}
rng = np.random.default_rng(0); P = rng.normal(size=(4000, 3)); P /= np.linalg.norm(P, axis=1)[:, None]
rho = np.geomspace(1e-6, 1, 4000); rng.shuffle(rho)
X = (rho*np.exp(1j*psi_spiral(rho, 2.0, 0.25)))[:, None]*P
i, j = rng.integers(0, 4000, 200000), rng.integers(0, 4000, 200000); m = i != j
Z = X[i[m]] - X[j[m]]; zz = np.einsum('nd,nd->n', Z, Z)
real = np.linalg.norm((rho[:, None]*P)[i[m]] - (rho[:, None]*P)[j[m]], axis=1)**2
print(f"Spitze, isotrope Spirale theta = 2: min |z.z| / |z|^2_reell = {np.min(np.abs(zz)/real):.3e} (200 000 Zufallspaare)")
