import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np, json
d = json.load(open('cone_curves.json')); plt.rcParams.update({'font.size': 9})
fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.0))
for ax, th, xl in [(axs[0], '20.0', (-35, 0.5)), (axs[1], '41.41', (-8, 0.5))]:
    for key, st, lab in [(f'{th}|0|0.0', 'C3-', '$m=0$, $L^2$'), (f'{th}|1|0.0', 'C0--', '$m=1$, $L^2$'), (f'{th}|0|-0.5', 'k-', '$m=0$, Energie')]:
        k = np.array([complex(a, b) for a, b in d[key]]); k = np.concatenate([k, np.conj(k)])
        k = k[k.imag >= -1e-9]; o = np.argsort(k.real)
        ax.plot(k.real, np.abs(k.imag), st, lw=1.2, label=lab) if 'Energie' not in lab else ax.plot(k.real, 0*k.real + 0.02, st, lw=3, label=lab)
    for name, kv, mk in [('Gold/Vak.', -11+1.2j, 'ks'), ('Gold/Wasser', (-11+1.2j)/1.77, 'kD'), ('Silber', -15+0.5j, 'k^')]:
        if kv.real > xl[0]: ax.plot(kv.real, kv.imag, mk, ms=4); ax.annotate(name, (kv.real, kv.imag), xytext=(3, 3), textcoords='offset points', fontsize=7)
    ax.set_xlim(*xl); ax.set_ylim(0, None); ax.set_title(f'Kegel, halber Öffnungswinkel {float(th):.1f}°', fontsize=8)
    ax.set_xlabel(r'Re $\kappa$'); ax.set_ylabel(r'Im $\kappa$')
axs[0].legend(fontsize=6.5, loc='upper left'); fig.tight_layout(); fig.savefig('fig_cone.pdf')
