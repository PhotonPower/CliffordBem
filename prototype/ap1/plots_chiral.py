import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
from corner_weights import S
plt.rcParams.update({'font.size': 9})
def crit_eps(chi, tau, u=1.0, al=np.pi/2):
    s2 = S(al, 0.5 + 1j*tau)**2; c2 = chi**2
    # A - B s2 + C s2^2 = 0 als Polynom in e (u fest)
    P = np.polynomial.Polynomial
    e = P([0, 1])
    A = ((1 + e)*(1 + u) - c2)**2; C = ((1 - e)*(1 - u) - c2)**2
    B = 2*(1 + e**2)*(1 + u**2) - 8*e*u + 4*c2*(3 - e*u) + 2*c2**2
    return (A - B*s2 + C*s2**2).roots()
fig, ax = plt.subplots(figsize=(6.2, 3.3))
taus = np.linspace(-15, 15, 6001)
for chi, st in [(0.0, 'k-'), (0.3, 'C0--'), (0.6, 'C3-.')]:
    pts = np.concatenate([crit_eps(chi, t) for t in taus])
    pts = pts[pts.imag >= 0]
    ax.plot(pts.real, pts.imag, st[:2], ls='none', marker='.', ms=0.8, label=f'$\\chi_1={chi}$')
    ax.plot([-1 + chi**2/2], [0], st[:2]+'o', ms=5)
ax.set_xlim(-8, 0.5); ax.set_ylim(0, 2.6); ax.set_xlabel(r'Re $\varepsilon_1$ ($\mu_1=1$, Außenraum Vakuum)'); ax.set_ylabel(r'Im $\varepsilon_1$')
ax.legend(fontsize=7, markerscale=10, loc='upper left'); fig.tight_layout(); fig.savefig('fig_chiral_corner.pdf')
for chi in [0.0, 0.3, 0.6]:
    r = crit_eps(chi, 0.0); print(chi, np.round(np.sort_complex(r), 4))
