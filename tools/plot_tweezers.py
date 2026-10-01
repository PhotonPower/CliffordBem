"""Kennlinien der optischen Pinzette aus apps/tweezers (CSV): axial Q_z(z) und seitlich Q_x(x) in der Gleichgewichtslage.
Aufruf: python3 tools/plot_tweezers.py axial.csv lateral_x.csv lateral_y.csv out.png
"""
import csv, sys
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
ax_f, lx_f, ly_f, out = sys.argv[1:5]
rd = lambda f: list(csv.DictReader(open(f)))
A, X, Y = rd(ax_f), rd(lx_f), rd(ly_f)
z = np.array([float(r['z_um']) for r in A]); Qz = np.array([float(r['Qz']) for r in A])
i = np.where((Qz[:-1] > 0) & (Qz[1:] <= 0))[0][0]; zeq = z[i] + (z[i + 1] - z[i]) * Qz[i] / (Qz[i] - Qz[i + 1]); kz = (Qz[i + 1] - Qz[i]) / (z[i + 1] - z[i])
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].plot(z, Qz, 'o-', ms=3); ax[0].axhline(0, color='k', lw=0.6); ax[0].axvline(zeq, color='g', ls='--', lw=0.8, label=f'Gleichgewicht z = {zeq:.2f} µm')
ax[0].set_xlabel('Lage des Teilchens hinter dem Fokus z (µm)'); ax[0].set_ylabel('$Q_z$'); ax[0].set_title('axial', fontsize=10); ax[0].legend(fontsize=8)
for D, lab in [(X, 'parallel zur Polarisation'), (Y, 'senkrecht zur Polarisation')]:
    x = np.array([float(r['x_um']) for r in D]); Q = np.array([float(r['Qx']) for r in D])
    ax[1].plot(x, Q, 'o-', ms=3, label=f'{lab} (Steifigkeit {-(Q[1] - Q[0]) / (x[1] - x[0]):.2f}/µm)')
ax[1].axhline(0, color='k', lw=0.6); ax[1].set_xlabel('seitliche Auslenkung (µm)'); ax[1].set_ylabel('$Q$ (rückstellend)'); ax[1].set_title(f'seitlich bei z = {zeq:.2f} µm', fontsize=10); ax[1].legend(fontsize=8)
fig.suptitle('Polystyrolkugel R = 0,25 µm in Wasser, λ = 1064 nm, NA 1,2', fontsize=10); fig.tight_layout(); fig.savefig(out, dpi=130)
print(f"Gleichgewicht z = {zeq:.3f} um, axiale Steifigkeit {-kz:.3f}/um")
