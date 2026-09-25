"""Abbildungen fuer docs/results_spectra.md."""
import csv, sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, 'tools'); import mie; from mie_spectrum import eps_tab
r = list(csv.DictReader(open('results/au40_water.csv')))
L = np.array([float(x['lambda_nm']) for x in r]); S = np.array([float(x['sigma_nm2']) for x in r])
Lf = np.linspace(440, 660, 200); a = 40.0
M = [mie.qext(2*np.pi*a/l, eps_tab('data/materials/Au_Johnson.yml', l), 1.33**2)*np.pi*a*a for l in Lf]
fig, ax = plt.subplots(figsize=(5.5, 3.6)); ax.plot(Lf, np.array(M)/1e3, 'k-', label='Mie'); ax.plot(L, S/1e3, 'o', label='BEM (1 280 Dreiecke)')
ax.set_xlabel('Wellenlänge (nm)'); ax.set_ylabel(r'$\sigma_{ext}$ (10$^3$ nm$^2$)'); ax.set_title('Goldkugel, Radius 40 nm, in Wasser'); ax.legend(); fig.tight_layout(); fig.savefig('docs/fig_au40_spectrum.png', dpi=130)
r = sorted(csv.DictReader(open('results/bornkuhn60_orient6.csv')), key=lambda x: float(x['lambda_nm']))
L = np.array([float(x['lambda_nm']) for x in r]); S = np.array([float(x['sigma_nm2']) for x in r]); CD = np.array([float(x['CD_nm2']) for x in r])
fig, ax = plt.subplots(2, 1, figsize=(5.5, 5), sharex=True)
ax[0].plot(L, S/1e3, 'o-'); ax[0].set_ylabel(r'$\sigma_{ext}$ (10$^3$ nm$^2$)'); ax[0].set_title('Born-Kuhn-Dimer (Au, 60°), Orientierungsmittel', fontsize=10)
ax[1].plot(L, CD/1e3, 's-', color='C3'); ax[1].axhline(0, color='k', lw=0.6); ax[1].set_ylabel(r'CD = $\sigma_+ - \sigma_-$ (10$^3$ nm$^2$)'); ax[1].set_xlabel('Wellenlänge (nm)')
fig.tight_layout(); fig.savefig('docs/fig_bornkuhn_cd.png', dpi=130)
