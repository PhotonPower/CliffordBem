import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np, json
from corner_weights import S
plt.rcParams.update({'font.size': 9})
fig, ax = plt.subplots(figsize=(6.2, 3.4))
G = json.load(open('corner_winding.json')); Re, Im, W = np.array(G['Re']), np.array(G['Im']), np.array(G['W'])
W = np.nan_to_num(np.round(W)); m = Im >= 0
ax.contourf(Re, Im[m], W[m], levels=[0.5, 1.5], colors=['#f4c7c3'])
t = np.linspace(-12, 12, 4001); al = np.pi/2
for c, st in [(0.5, 'C3-'), (0.35, 'C1--'), (0.2, 'C0-.'), (0.05, 'C2:')]:
    s = S(al, c + 1j*t)
    for sg in [1, -1]:
        k = (sg*s - 1)/(sg*s + 1); k = k[k.imag >= 0]
        ax.plot(k.real, k.imag, st, lw=1.2, label=f'Re$\\,a={c}$' if sg == 1 else None)
for name, k, mk in [('Gold/Vak.', -11+1.2j, 'ks'), ('Gold/Wasser', (-11+1.2j)/1.77, 'kD'), ('$-3+0{,}3i$', -3+0.3j, 'ko'), ('$-1{,}5+0{,}2i$', -1.5+0.2j, 'ko')]:
    ax.plot(k.real, k.imag, mk, ms=4); ax.annotate(name, (k.real, k.imag), textcoords='offset points', xytext=(3, 4), fontsize=7)
ax.axvspan(-3, -1/3, ymin=0, ymax=0.02, color='k')
ax.set_xlim(-12, 0.5); ax.set_ylim(0, 2.6); ax.set_xlabel(r'Re $\kappa=\varepsilon_1/\varepsilon_2$'); ax.set_ylabel(r'Im $\kappa$')
ax.legend(fontsize=7, loc='upper left'); fig.tight_layout(); fig.savefig('fig_corner.pdf')
