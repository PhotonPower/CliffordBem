"""Fluoreszenzverstaerkung eines fest, aber zufaellig orientierten Emitters vor einer Kugel: BEM (apps/dipole mit
--lambda-exc, CSV) gegen Mie (Nahfeld tools/mie_nearfield.py bei der Anregung, Raten tools/mie_dipole.py bei der Emission).
F/F0 = sum_a |E_a|^2 q_a / (|E0|^2 q0), q_a = gamma_rad,a / (gamma_tot,a + (1 - q0)/q0); die Raten haengen nicht von q0 ab,
deshalb wird F fuer mehrere q0 aus einer Rechnung gewonnen.
Aufruf: python3 tools/plot_fluorescence.py results/fluor_sphere_550_570.csv docs/fig_fluorescence_sphere.png
        [--radius 20 --lexc 550 --lem 570 --nbg 1.33 --extra scratch/fluor_d10.csv]
"""
import csv, os, sys
from collections import defaultdict
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mie_dipole as md, mie_nearfield as mn
from mie_spectrum import eps_tab

fn, out = sys.argv[1], sys.argv[2]
arg = lambda k, d: (sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d)
R, Le, Lm, nb = float(arg('--radius', 20)), float(arg('--lexc', 550)), float(arg('--lem', 570)), float(arg('--nbg', 1.33))
Q0 = [1.0, 0.1, 0.01]
rows = list(csv.DictReader(open(fn)))
if '--extra' in sys.argv:
    hdr = rows[0].keys() if rows else None
    for line in open(arg('--extra', '')):
        v = line.strip().split(',')
        if len(v) == 12: rows.append(dict(zip(list(hdr), v)))
by = defaultdict(dict)
for r in rows: by[float(r['dist_nm'])][r['orientation']] = r
mat = os.path.join(os.path.dirname(__file__), '..', 'data', 'materials', 'Au_Johnson.yml')
ee, em = eps_tab(mat, Le), eps_tab(mat, Lm); oe, om = 2 * np.pi * R / Le, 2 * np.pi * R / Lm


def F_of(exc, rates, q0):                                   # exc[a], rates[a] = (tot, rad), a = x, y, z
    return sum(exc[a] * rates[a][1] / (rates[a][0] + (1 - q0) / q0) for a in range(3)) / q0


def mie(d):
    E, _ = mn.fields(oe, [1.0], [ee], nb ** 2, np.array([1 + d / R, 0, 0]), s=0)
    rr = md.rates_mp(om, 1.0, em, nb ** 2, 1 + d / R, 'radial'); rt = md.rates_mp(om, 1.0, em, nb ** 2, 1 + d / R, 'tangential')
    return list(np.abs(E) ** 2), [rr, rt, rt]


ds = sorted(by)
bem = {d: ([float(by[d][a]['exc']) for a in 'xyz'], [(float(by[d][a]['gamma_tot']), float(by[d][a]['gamma_rad'])) for a in 'xyz']) for d in ds}
ref = {d: mie(d) for d in ds}
dg = np.geomspace(1.5, 40, 16); curve = [mie(x) for x in dg]
fig, ax = plt.subplots(figsize=(5.8, 4.0))
for i, q0 in enumerate(Q0):
    ax.semilogx(dg, [F_of(e, r, q0) for e, r in curve], '-', color=f'C{i}', lw=1, label=f'Mie, q₀ = {q0:g}')
    ax.semilogx(ds, [F_of(*bem[d], q0) for d in ds], 'o', color=f'C{i}', mfc='none')
ax.axhline(1, color='k', lw=0.5); ax.set_xlabel('Abstand zur Oberfläche (nm)'); ax.set_ylabel(r'Fluoreszenzverstärkung $F/F_0$')
ax.set_yscale('log'); ax.legend(fontsize=7); ax.grid(True, which='both', alpha=0.3)
ax.set_title(f'Goldkugel {R:.0f} nm, Wasser, Anregung {Le:.0f} nm, Emission {Lm:.0f} nm (Kreise: BEM)', fontsize=9)
fig.tight_layout(); fig.savefig(out, dpi=130)
print("d (nm) | Anregung x BEM/Mie | F/F0 fuer q0 = 1, 0,1, 0,01 (BEM / Mie)")
for d in ds:
    eb, rb = bem[d]; er, rr = ref[d]
    print(f"  {d:5.1f} | {eb[0]:7.3f} {er[0]:7.3f} | " + " | ".join(f"{F_of(eb, rb, q):8.3f} {F_of(er, rr, q):8.3f}" for q in Q0))
