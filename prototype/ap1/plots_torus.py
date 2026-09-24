import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np, json
plt.rcParams.update({'font.size': 9})
r = json.load(open('bor_eps.json'))
es = sorted(set(float(k.split('|')[1]) for k in r))
fig, ax = plt.subplots(figsize=(6.2, 3.0))
for v, lab, st in [('maxwell','Original, Maxwell-Raum','k-o'), ('dirac','Original, voller Raum','C0--x'),
                   ('dual_J','Dual, $W=J^{-1}$','C3-s'), ('dual_W1','Dual, $W=1$','C2-^')]:
    y = [max(min(r[f'{v}|{e:g}'].values()), 1e-5) for e in es]
    ax.semilogy(es, y, st, label=lab, ms=3, lw=1)
ax.set_xlabel(r'$\varepsilon_1$ (verlustfrei, $\varepsilon_2=1$, $\omega=10^{-3}$)'); ax.set_ylabel(r'$\min_m \sin\theta_{\min}(m)$')
ax.set_ylim(1e-5, 1.5); ax.legend(fontsize=7, loc='lower left'); fig.tight_layout(); fig.savefig('fig_torus_eps.pdf')
L = json.load(open('bor_lf.json'))
fig, ax = plt.subplots(figsize=(6.2, 2.7))
for mat, c in [('glas','C0'), ('gold','C1')]:
    for v, ls, lab in [('maxwell','-','Maxwell'), ('dirac','--','voll'), ('dual_W1',':','dual, $W=1$')]:
        oms = [1e-6,1e-4,1e-2,0.1,0.3,1.0]
        y = [min(L[f'{mat}_{v}|{om:g}'].values()) for om in oms]; y0 = min(L[f'{mat}_{v}|0'].values())
        ax.semilogx(oms, y, color=c, ls=ls, label=f'{mat.capitalize()}, {lab}')
        ax.plot([3e-7], [y0], marker='o', color=c, ms=3)
ax.set_xlabel(r'$\omega$ (Punkte links: $\omega=0$)'); ax.set_ylabel(r'$\min_m\sin\theta_{\min}(m)$'); ax.set_ylim(0, 1)
ax.legend(fontsize=7, ncol=2, loc='center left'); fig.tight_layout(); fig.savefig('fig_torus_lf.pdf')
