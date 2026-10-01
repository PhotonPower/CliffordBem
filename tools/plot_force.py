"""Kraftkarten aus apps/nearfield --particle (CSV): Kraftquerschnitt |sigma_F| und chirale Kraft F(A_c) - F(-A_c), dazu der
Vergleich der chiralen Kraft mit dem Gradienten der optischen Chiralitaet C (aus derselben Karte, zentrale Differenzen).
Aufruf: python3 tools/plot_force.py results/force_dimer_580.csv docs/fig_force_dimer.png [--title "..."]
"""
import csv, sys
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
fn, out = sys.argv[1], sys.argv[2]
title = sys.argv[sys.argv.index('--title') + 1] if '--title' in sys.argv else ''
R = list(csv.DictReader(open(fn)))
x = np.array([float(r['x_nm']) for r in R]); z = np.array([float(r['z_nm']) for r in R])
ux, uz = np.unique(x), np.unique(z); nx, nz = len(ux), len(uz)
G = lambda k: np.array([float(r[k]) for r in R]).reshape(nz, nx)
ok = G('force_ok') > 0.5
Fx, Fz, Cx, Cz, C = G('Fx_nm2'), G('Fz_nm2'), G('Fcx_nm2'), G('Fcz_nm2'), G('chirality')
F = np.hypot(Fx, Fz); Fc = np.hypot(Cx, Cz)
# Gradient von C auf dem Gitter (nur wo die Nachbarn gueltig sind)
dx, dz = ux[1] - ux[0], uz[1] - uz[0]
gCx = np.full_like(C, np.nan); gCz = np.full_like(C, np.nan)
gCx[:, 1:-1] = (C[:, 2:] - C[:, :-2]) / (2 * dx); gCz[1:-1, :] = (C[2:, :] - C[:-2, :]) / (2 * dz)
inner = ok.copy(); inner[:, 1:-1] &= ok[:, 2:] & ok[:, :-2]; inner[1:-1, :] &= ok[2:, :] & ok[:-2, :]
inner[:, 0] = inner[:, -1] = False; inner[0, :] = inner[-1, :] = False
a = np.stack([Cx[inner], Cz[inner]], 1); b = np.stack([gCx[inner], gCz[inner]], 1)
cosang = np.sum(a * b, 1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + 1e-300)
coef = np.sum(a * b) / np.sum(b * b); resid = np.linalg.norm(a - coef * b) / np.linalg.norm(a)
w = np.linalg.norm(a, axis=1)
print(f"{inner.sum()} Punkte mit Gradient: chirale Kraft gegen grad C -- Kosinus gewichtet {np.sum(w * cosang) / np.sum(w):+.4f}, "
      f"Median {np.median(cosang):+.4f}; Faktor {coef:+.4e} nm^3, Rest nach Anpassung {resid:.1%}")
F[~ok] = np.nan; Fc[~ok] = np.nan
fig, ax = plt.subplots(1, 2, figsize=(10, 3.9)); ext = [ux[0], ux[-1], uz[0], uz[-1]]
im = ax[0].imshow(np.log10(F), origin='lower', extent=ext, cmap='viridis', aspect='equal'); fig.colorbar(im, ax=ax[0], label=r'$\log_{10}|\sigma_F|$ (nm$^2$)')
im2 = ax[1].imshow(np.log10(Fc), origin='lower', extent=ext, cmap='magma', aspect='equal'); fig.colorbar(im2, ax=ax[1], label=r'$\log_{10}|F_{chir}|$ (nm$^2$)')
s = max(1, nx // 25)
for A_, B_, axx, col in [(Fx, Fz, ax[0], 'w'), (Cx, Cz, ax[1], 'c')]:
    U = np.where(ok, A_, np.nan); V = np.where(ok, B_, np.nan); n_ = np.hypot(U, V)
    axx.quiver(x.reshape(nz, nx)[::s, ::s], z.reshape(nz, nx)[::s, ::s], (U / n_)[::s, ::s], (V / n_)[::s, ::s], color=col, scale=30, width=0.003)
ax[0].set_title('Kraft auf das Teilchen (Richtung, Betrag)', fontsize=10); ax[1].set_title('chirale Kraft F(+) − F(−)', fontsize=10)
for axx in ax: axx.set_xlabel('x (nm)'); axx.set_ylabel('z (nm)')
if title: fig.suptitle(title, fontsize=10)
fig.tight_layout(); fig.savefig(out, dpi=130)
