import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
from summarize import load
plt.rcParams.update({'font.size': 9, 'figure.figsize': (6.2, 3.1)})
E = load('eps')
fig, ax = plt.subplots()
for n, lab, st in [('maxwell','Original, Maxwell-Raum','k-'), ('dirac','Original, voller Raum','C0--'),
                   ('dirac_dual','Dual, $W=J^{-1}$','C3-'), ('lossy_dual','Dual, $W=J^{-1}$, Im$\\,\\varepsilon_1=0{,}1$','C3:'),
                   ('dirac_dual_W1','Dual, $W=1$','C2-'), ('dirac_energy','Original, Energiewahl','C1-.')]:
    x, v = E[n]; ax.semilogy(x, v, st, label=lab, lw=1.2)
for l in range(1, 6):
    ax.axvline(-(l+1)/l, color='0.6', lw=0.6, ls=':'); ax.axvline(-l/(l+1), color='0.6', lw=0.6, ls='-.')
ax.set_xlabel(r'$\varepsilon_1$ (verlustfrei, $\varepsilon_2=1$, $\omega a = 0{,}5$)'); ax.set_ylabel(r'$\sin\theta_{\min}$')
ax.legend(fontsize=7, loc='lower left', ncol=2); ax.set_ylim(5e-4, 1.5); fig.tight_layout(); fig.savefig('fig_eps.pdf')
F = load('freq')
fig, ax = plt.subplots()
for n, lab, st in [('maxwell','Original, Maxwell-Raum','k-'), ('dirac','Original, voller Raum','C0--'),
                   ('dirac_dual','Dual, $W=J^{-1}$','C3-'), ('dirac_dual_W1','Dual, $W=1$','C2-'),
                   ('dirac_bad','Original, $b=-\\sqrt{r_\\varepsilon}$','C7-')]:
    x, v = F[n]; ax.plot(x, v, st, label=lab, lw=1.2)
ax.set_xlabel(r'$\omega a$ (Glas, $\varepsilon_1=2{,}25$)'); ax.set_ylabel(r'$\sin\theta_{\min}$'); ax.set_ylim(0, 1.02)
ax.legend(fontsize=7, loc='upper right'); fig.tight_layout(); fig.savefig('fig_freq.pdf')
L = load('lowfreq')
fig, ax = plt.subplots(figsize=(6.2, 2.6))
for mat, c in [('glas','C0'),('gold','C1')]:
    for n, ls, lab in [('maxwell','-','Maxwell'),('dirac','--','voll'),('dirac_dual',':','dual')]:
        x, v = L[f'{mat}_{n}']; ax.semilogx(x, v, color=c, ls=ls, label=f'{mat.capitalize()}, {lab}')
ax.set_xlabel(r'$\omega a$'); ax.set_ylabel(r'$\sin\theta_{\min}$'); ax.set_ylim(0, 1); ax.legend(fontsize=7, ncol=2)
fig.tight_layout(); fig.savefig('fig_lowfreq.pdf')

import json
r = json.load(open('sweep_torus.json'))
fig, ax = plt.subplots(figsize=(6.4, 3.3))
for mat, c in [('glas','C0'), ('gold','C1')]:
    for n, ls, lab in [('maxwell','-','Maxwell'), ('dirac','--','voll, Standard'), ('dirac_energy','-.','voll, Energie'), ('dual_W1',':','dual, $W=1$')]:
        pts = sorted((float(k.split('|')[1]), v[0][0]) for k, v in r.items() if k.startswith(f'{mat}_{n}|'))
        ax.semilogx([p[0] for p in pts], [p[1] for p in pts], color=c, ls=ls, label=f'{mat.capitalize()}, {lab}')
ax.set_yscale('log'); ax.set_ylim(5e-3, 1.2)
ax.set_xlabel(r'$\omega R$ (Torus $R=1$, $r=0{,}4$)'); ax.set_ylabel(r'$\sin\theta_{\min}$'); ax.legend(fontsize=6.5, ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.28), frameon=False)
fig.tight_layout(); fig.savefig('fig_torus.pdf')
