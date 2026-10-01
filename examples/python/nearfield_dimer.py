"""Nahfeldkarte |E|^2/|E0|^2 eines Gold-Dimers (zwei Kugeln, Durchmesser 40 nm, Spalt 2 nm) in Wasser, Polarisation
entlang der Dimerachse; dazu die optische Kraft auf jede Kugel (Spannungstensor auf Parallelflaechen) und die Summe gegen
die Impulsbilanz im Fernfeld.

    PYTHONPATH=build/python python3 examples/python/nearfield_dimer.py --n 8 --lambda 580 --plot nf.png
"""
import argparse
import time

import numpy as np

import cliffordbem as cb

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--n", type=int, default=8)
ap.add_argument("--lambda", dest="lam", type=float, default=580.0)
ap.add_argument("--unit", type=float, default=20.0, help="Kugelradius in nm")
ap.add_argument("--gap", type=float, default=2.0, help="Spalt in nm")
ap.add_argument("--grid", default="121,61", help="Punkte in x und z")
ap.add_argument("--plot")
a = ap.parse_args()

om = cb.omega_from_wavelength(a.lam, a.unit)
water = cb.Medium(eps=1.33 ** 2)
gold = cb.make_material("Au").medium(a.lam)
c = 1 + 0.5 * a.gap / a.unit                                              # Mittelpunkte bei +-c
s = cb.make_icosphere(a.n)
bodies = [cb.translated(s, (-c, 0, 0)), cb.translated(s, (c, 0, 0))]

t0 = time.time()
P = cb.ScatteringProblem(bodies, [gold, gold], om, outer=water)
d, p = (0, 0, 1), (1, 0, 0)
r = P.solve_plane_wave(d, p)
print(f"{len(P.mesh)} Dreiecke, sigma_ext = {r.sigma_ext * a.unit ** 2:.1f} nm^2, {r.iterations} Iterationen, {time.time() - t0:.1f} s")

# Karte in der xz-Ebene (ab 2000 Punkten automatisch als rechteckige H-Matrix)
nx, nz = (int(v) for v in a.grid.split(","))
x = np.linspace(-2.5, 2.5, nx)
z = np.linspace(-1.5, 1.5, nz)
X, Zg = np.meshgrid(x, z)
pts = np.column_stack([X.ravel(), np.zeros(X.size), Zg.ravel()])
t0 = time.time()
nf = cb.exterior_near_field(P.mesh, r.h, water, om, d, p, pts)
enh = np.where(nf["inside"] | nf["too_close"], np.nan, nf["enhancement"]).reshape(nz, nx)
gap = cb.exterior_near_field(P.mesh, r.h, water, om, d, p, np.array([0.0, 0.0, 0.0]))
print(f"Nahfeld an {len(pts)} Punkten in {time.time() - t0:.1f} s; Spaltmitte |E|^2/|E0|^2 = {gap['enhancement'][0]:.0f}")

# Kraefte: je Kugel ueber die Parallelflaeche (Abstand kleiner als der halbe Spalt), Summe gegen das Fernfeld
ev = cb.plane_wave_evaluator(P.mesh, r.h, water, om, d, p)
delta = 0.25 * a.gap / a.unit
F = [cb.force_on_offset(ev, water, b, delta) for b in bodies]
Ff = cb.force_from_far_field(P.mesh, r.h, water, om, d, p, r.sigma_ext)
for i, f in enumerate(F):
    print(f"Kraft auf Kugel {i}: ({f[0]:+.4f}, {f[1]:+.4f}, {f[2]:+.4f})")
print(f"Summe: z = {F[0][2] + F[1][2]:.4f}, Fernfeld-Impulsbilanz: z = {Ff[2]:.4f}  (Kugeln ziehen sich an: F0_x > 0, F1_x < 0)")

if a.plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    im = ax.pcolormesh(x * a.unit, z * a.unit, np.log10(enh), shading="auto", cmap="inferno")
    fig.colorbar(im, ax=ax, label=r"$\log_{10}|E|^2/|E_0|^2$")
    ax.set_aspect("equal")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("z (nm)")
    ax.set_title(f"Au-Dimer, Spalt {a.gap:g} nm, {a.lam:g} nm")
    fig.tight_layout()
    fig.savefig(a.plot, dpi=150)
    print("Abbildung:", a.plot)
