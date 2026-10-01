"""Eigene einfallende Felder als Python-Klassen (v0.43): eine optische Stehwelle (Gitterfalle) und ein Bessel-Strahl auf einer
Goldkugel (Durchmesser 40 nm) in Wasser.

1) Stehwelle aus zwei gegenlaeufigen ebenen Wellen, Bauch bei z = s. Die Kraft ist quadratisch im Feld, das Feld
   e^{-iks} F_auf + e^{iks} F_ab; daraus F_z(s) = C + A cos 2ks + B sin 2ks, und die Symmetrie erzwingt C = A = 0:
   F_z(s) = B sin(2ks) exakt (auch fuer grosse Kugeln, nicht nur in Dipolnaeherung). Das Skript prueft das.
2) Bessel-Strahl (Kegel ebener Wellen, AngularSpectrumField.bessel): Maxwell-Pruefung, Querkraft F_x bei seitlich
   versetzter Kugel (Strahl verschoben, gleiches Problem, Krylov-Recycling).

    PYTHONPATH=build/python python3 examples/python/custom_field.py --n 6 --plot falle.png
"""
import argparse
import time

import numpy as np

import cliffordbem as cb

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--n", type=int, default=6)
ap.add_argument("--lambda", dest="lam", type=float, default=600.0)
ap.add_argument("--unit", type=float, default=20.0, help="Kugelradius in nm")
ap.add_argument("--points", type=int, default=9, help="Positionen je Kurve")
ap.add_argument("--plot")
a = ap.parse_args()

om = cb.omega_from_wavelength(a.lam, a.unit)
water = cb.Medium(eps=1.33 ** 2)
gold = cb.make_material("Au").medium(a.lam)
k = water.k(om).real
mesh = cb.make_icosphere(a.n)
t0 = time.time()
P = cb.ScatteringProblem(mesh, gold, om, outer=water)
P.use_recycling()                                                     # viele rechte Seiten fuer dasselbe Problem
print(f"{len(mesh)} Dreiecke, lambda = {a.lam:g} nm, k = {k:.4f} / Einheit, Aufbau {time.time() - t0:.1f} s")


def force(field, R=1.5):
    """Kraft auf die Kugel im Feld: rechte Seite, Loesung, Spannungstensor auf der Kugel mit Radius R."""
    b = field.project(P.mesh, water)
    r = P.solve_rhs(b)
    return cb.force_on_sphere(cb.near_field_evaluator(P.mesh, r.h, b, water, om, field), water, (0, 0, 0), R), r.iterations


# --- 1) Stehwelle --------------------------------------------------------------------------------------------------------
class StandingWave(cb.CustomField):
    """E = 2 x cos(k (z - s)), H = 2i Z y sin(k (z - s)): Bauch bei z = s; Bezug: eine einzelne Welle (|E0| = 1)."""

    def __init__(self, medium, omega, s):
        self.k, self.Z, self.s = medium.k(omega), np.sqrt(medium.eps / medium.mu), s

    def fields(self, x):
        u = self.k * (x[:, 2] - self.s)
        E = np.zeros((len(x), 3), complex)
        H = np.zeros_like(E)
        E[:, 0] = 2 * np.cos(u)
        H[:, 1] = 2j * self.Z * np.sin(u)
        return E, H


sw0 = StandingWave(water, om, 0.0)
print(f"Maxwell-Residuum der Stehwelle: {cb.maxwell_residual(sw0, water, om, np.random.default_rng(0).normal(size=(10, 3))):.1e}")
# Gegenprobe: dieselbe Welle aus zwei Feldern des Kerns
ref = cb.PlaneWaveField(water, om, (0, 0, 1), (1, 0, 0)) + cb.PlaneWaveField(water, om, (0, 0, -1), (1, 0, 0))
x = np.random.default_rng(1).normal(size=(5, 3))
print(f"|StandingWave - (PlaneWaveField + PlaneWaveField)| = {np.abs(sw0.eval(x)[0] - ref.eval(x)[0]).max():.1e}")

s = np.linspace(-np.pi / (2 * k), np.pi / (2 * k), a.points)        # eine halbe Gitterperiode lambda/2
Fz = []
t0 = time.time()
for si in s:
    F, its = force(StandingWave(water, om, si))
    Fz.append(F[2])
    print(f"  Bauch bei s = {si * a.unit:+7.2f} nm: F_z = {F[2]:+.5f}   ({its} It.)", flush=True)
Fz = np.array(Fz)
B = np.linalg.lstsq(np.sin(2 * k * s)[:, None], Fz, rcond=None)[0][0]
dev = np.abs(Fz - B * np.sin(2 * k * s)).max() / abs(B)
print(f"F_z(s) = B sin(2ks) mit B = {B:+.5f}; groesste Abweichung {dev:.1e} von |B|  ({time.time() - t0:.1f} s)")
# Bauch oberhalb der Kugel (s > 0) und F_z > 0: die Kugel wird zum Bauch gezogen (sucht hohe Intensitaet)
print("Die Kugel wird zum " + ("Bauch (hohe Intensitaet)" if B > 0 else "Knoten (niedrige Intensitaet)") + " gezogen.")

# --- 2) Bessel-Strahl ------------------------------------------------------------------------------------------------------
alpha = np.deg2rad(35)
bes0 = cb.AngularSpectrumField.bessel(water, om, alpha, nphi=64)
print(f"\nBessel-Strahl, Kegelwinkel 35 Grad, Kernradius ~ 2.405/(k sin a) = {2.405 / (k * np.sin(alpha)) * a.unit:.0f} nm; "
      f"Maxwell-Residuum {cb.maxwell_residual(bes0, water, om, np.random.default_rng(2).normal(size=(10, 3))):.1e}")
xs = np.linspace(0, 2.405 / (k * np.sin(alpha)), a.points)
Fx = []
for xi in xs:
    beam = cb.AngularSpectrumField.bessel(water, om, alpha, nphi=64, focus=(-xi, 0, 0))   # Kugel bei x = xi relativ zum Strahl
    F, its = force(beam)
    Fx.append(F[0])
    print(f"  Versatz {xi * a.unit:6.1f} nm: F = ({F[0]:+.5f}, {F[1]:+.5f}, {F[2]:+.5f})   ({its} It.)", flush=True)

if a.plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ss = np.linspace(s[0], s[-1], 200)
    ax[0].plot(ss * a.unit, B * np.sin(2 * k * ss), "k-", label=r"$B\,\sin 2ks$")
    ax[0].plot(s * a.unit, Fz, "o", label="BEM")
    ax[0].set_xlabel("Lage des Bauchs s (nm)")
    ax[0].set_ylabel("$F_z$")
    ax[0].set_title("Stehwelle")
    ax[0].legend()
    ax[1].plot(xs * a.unit, Fx, "o-")
    ax[1].axhline(0, color="0.6", lw=0.8)
    ax[1].set_xlabel("seitlicher Versatz (nm)")
    ax[1].set_ylabel("$F_x$")
    ax[1].set_title("Bessel-Strahl")
    fig.tight_layout()
    fig.savefig(a.plot, dpi=150)
    print("Abbildung:", a.plot)
