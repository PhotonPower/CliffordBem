"""Nahfeldkarten aus apps/nearfield (CSV): Feldverstaerkung |E|^2/|E_0|^2 (logarithmisch) und optische Chiralitaet C/C_0.
Punkte innerhalb der Koerper bzw. zu nah am Netz sind ausgeblendet.
Aufruf: python3 tools/plot_nearfield.py datei.csv ausgabe.png [--plane xz] [--title "..."]
"""
import csv, sys
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

fn, out = sys.argv[1], sys.argv[2]
plane = sys.argv[sys.argv.index('--plane') + 1] if '--plane' in sys.argv else 'xz'
title = sys.argv[sys.argv.index('--title') + 1] if '--title' in sys.argv else ''
R = list(csv.DictReader(open(fn)))
ax1, ax2 = {'xz': ('x_nm', 'z_nm'), 'xy': ('x_nm', 'y_nm'), 'yz': ('y_nm', 'z_nm')}[plane]
a = np.array([float(r[ax1]) for r in R]); b = np.array([float(r[ax2]) for r in R])
ua, ub = np.unique(a), np.unique(b)
mask = np.array([r['inside'] == '1' or r.get('too_close', '0') == '1' for r in R])
enh = np.array([float(r['enhancement']) for r in R]); chi = np.array([float(r['chirality']) for r in R])
enh[mask] = np.nan; chi[mask] = np.nan
E = enh.reshape(len(ub), len(ua)); C = chi.reshape(len(ub), len(ua))
fig, axs = plt.subplots(1, 2, figsize=(10, 3.9))
ext = [ua[0], ua[-1], ub[0], ub[-1]]
im = axs[0].imshow(np.log10(E), origin='lower', extent=ext, cmap='inferno', aspect='equal')
fig.colorbar(im, ax=axs[0], label=r'$\log_{10}\,|E|^2/|E_0|^2$')
cl = np.nanmax(np.abs(C))
im2 = axs[1].imshow(C, origin='lower', extent=ext, cmap='RdBu_r', vmin=-cl, vmax=cl, aspect='equal')
fig.colorbar(im2, ax=axs[1], label=r'$C/C_0$')
for x in axs: x.set_xlabel(f'{ax1[0]} (nm)'); x.set_ylabel(f'{ax2[0]} (nm)')
axs[0].set_title('Feldverstärkung', fontsize=10); axs[1].set_title('optische Chiralität', fontsize=10)
if title: fig.suptitle(title, fontsize=10)
fig.tight_layout(); fig.savefig(out, dpi=130)
print(f"max |E|^2/|E0|^2 = {np.nanmax(E):.1f}, C/C0 in [{np.nanmin(C):.2f}, {np.nanmax(C):.2f}]")
