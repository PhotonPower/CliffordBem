"""Extinktionsspektrum einer Goldkugel (Durchmesser 40 nm) in Wasser gegen Mie -- wie
``spectrum --sphere 12 --unit 20 --materials Au --nbg 1.33 --lambda 450:650:10``, aber als Skript.

    PYTHONPATH=build/python python3 examples/python/spectrum_gold_sphere.py --n 8 --plot spektrum.png
"""
import argparse
import os
import sys

import numpy as np

import cliffordbem as cb

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
import mie  # noqa: E402  Referenz (Bohren-Huffman), Kugelradius 1

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--n", type=int, default=8, help="Ikosaederkugel mit 20 n^2 Dreiecken")
ap.add_argument("--lambda", dest="lam", default="450:650:20", help="von:bis:Schritt in nm")
ap.add_argument("--unit", type=float, default=20.0, help="Kugelradius in nm (Laengeneinheit)")
ap.add_argument("--nbg", type=float, default=1.33, help="Brechzahl des Hintergrunds")
ap.add_argument("--csv", help="Ergebnisse als CSV")
ap.add_argument("--plot", help="Abbildung (benoetigt matplotlib)")
a = ap.parse_args()

lo, hi, step = (float(x) for x in a.lam.split(":"))
lams = np.arange(lo, hi + 0.5 * step, step)
gold = cb.make_material("Au")
sphere = cb.make_icosphere(a.n)                                          # Radius 1 = unit nm

print(f"cliffordbem {cb.__version__}, {len(sphere)} Dreiecke, {cb.omp_threads()} Thread(s)")
print(f"{'lambda nm':>10} {'BEM nm^2':>12} {'Mie nm^2':>12} {'Abw. %':>8} {'It.':>5}")


def report(i, row):
    lam = row["lambda_nm"]
    q_mie = mie.qext(cb.omega_from_wavelength(lam, a.unit), gold.eps(lam), a.nbg ** 2)
    s_mie = q_mie * np.pi * a.unit ** 2
    print(f"{lam:10.1f} {row['sigma_nm2']:12.2f} {s_mie:12.2f} {100 * (row['sigma_nm2'] / s_mie - 1):+8.2f} {row['iterations']:5d}", flush=True)


res = cb.spectrum(sphere, gold, lams, unit=a.unit, n_bg=a.nbg, callback=report)
s_mie = np.array([mie.qext(cb.omega_from_wavelength(l, a.unit), gold.eps(l), a.nbg ** 2) for l in lams]) * np.pi * a.unit ** 2

if a.csv:
    np.savetxt(a.csv, np.column_stack([lams, res["sigma_nm2"], s_mie, res["forward"].real, res["forward"].imag]), delimiter=",",
               header="lambda_nm,sigma_bem_nm2,sigma_mie_nm2,S_re,S_im", comments="", fmt="%.10g")
if a.plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(lams, s_mie, "k-", label="Mie")
    ax.plot(lams, res["sigma_nm2"], "o", label=f"BEM, {len(sphere)} Dreiecke")
    ax.set_xlabel("Wellenlaenge (nm)")
    ax.set_ylabel(r"$\sigma_\mathrm{ext}$ (nm$^2$)")
    ax.set_title(f"Au-Kugel, d = {2 * a.unit:g} nm, n = {a.nbg}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(a.plot, dpi=150)
    print("Abbildung:", a.plot)
