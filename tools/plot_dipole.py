"""Zerfallsraten eines Emitters vor einer Kugel gegen den Abstand: BEM (apps/dipole, CSV) gegen die Reihenloesung.
Aufruf: python3 tools/plot_dipole.py results/dipole_au20_650.csv docs/fig_dipole_au20.png --radius 20 --lambda 650 --nbg 1.33
"""
import csv, os, sys
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mie_dipole as md
from mie_spectrum import eps_tab
fn, out = sys.argv[1], sys.argv[2]
arg = lambda k, d: float(sys.argv[sys.argv.index(k) + 1]) if k in sys.argv else d
R, L, nb = arg('--radius', 20.0), arg('--lambda', 650.0), arg('--nbg', 1.33)
rows = [r for r in csv.DictReader(open(fn)) if r['orientation'] == 'average']
d = np.array([float(r['dist_nm']) for r in rows]); gt = np.array([float(r['gamma_tot']) for r in rows]); gr = np.array([float(r['gamma_rad']) for r in rows])
e = eps_tab(os.path.join(os.path.dirname(__file__), '..', 'data', 'materials', 'Au_Johnson.yml'), L); om = 2 * np.pi * R / L
df = np.geomspace(1.5, 40, 14); T, Rr = [], []
for x in df:
    a = [md.rates_mp(om, 1.0, e, nb ** 2, 1 + x / R, o) for o in ['radial', 'tangential']]
    T.append((a[0][0] + 2 * a[1][0]) / 3); Rr.append((a[0][1] + 2 * a[1][1]) / 3)
T, Rr = np.array(T), np.array(Rr)
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.8))
ax[0].loglog(df, T, 'k-', lw=1, label='gesamt (Reihe)'); ax[0].loglog(df, Rr, 'k--', lw=1, label='strahlend (Reihe)')
ax[0].loglog(d, gt, 'C3o', label='gesamt (BEM)'); ax[0].loglog(d, gr, 'C0s', label='strahlend (BEM)')
ax[0].set_xlabel('Abstand zur Oberfläche (nm)'); ax[0].set_ylabel(r'$\gamma/\gamma_0$ (orientierungsgemittelt)'); ax[0].legend(fontsize=7); ax[0].grid(True, which='both', alpha=0.3)
ax[1].semilogx(df, Rr / T, 'k-', lw=1, label='Reihe'); ax[1].semilogx(d, gr / gt, 'C3o', label='BEM (n = 12, verdichtet)')
ax[1].set_xlabel('Abstand zur Oberfläche (nm)'); ax[1].set_ylabel(r'Quantenausbeute $q$ ($q_0 = 1$)'); ax[1].legend(fontsize=7); ax[1].grid(True, which='both', alpha=0.3)
fig.suptitle(f'Emitter vor einer Goldkugel ({R:.0f} nm), Wasser, {L:.0f} nm', fontsize=10)
fig.tight_layout(); fig.savefig(out, dpi=130)
for x, a, b in zip(d, gt, gr):
    a0 = [md.rates_mp(om, 1.0, e, nb ** 2, 1 + x / R, o) for o in ['radial', 'tangential']]
    t0 = (a0[0][0] + 2 * a0[1][0]) / 3; r0 = (a0[0][1] + 2 * a0[1][1]) / 3
    print(f"d = {x:4.0f} nm: gesamt {a:9.3f} (Reihe {t0:9.3f}, {a / t0 - 1:+.1%}), strahlend {b:7.3f} ({r0:7.3f}, {b / r0 - 1:+.1%}), q {b / a:.4f} ({r0 / t0:.4f})")
