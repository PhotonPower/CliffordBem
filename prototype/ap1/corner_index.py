"""Windungszahl von det sigma_Mellin(T_1)(tau; kappa) ueber tau in R (L^2-Gerade c = 1/2),
alpha = pi/2. Nichtverschwindende Windungszahl = Index != 0 => T_1 nicht invertierbar."""
import numpy as np, json
from corner import M_E, Jfull
alpha = np.pi/2
taus = np.linspace(-9, 9, 721)
Ms = np.array([M_E(alpha, t) for t in taus]); I = np.eye(16)
def winding(kappa, a=None, b=None):
    J = Jfull(alpha, kappa, 1.0, a, b)
    T = 0.5*(I + Ms) + 0.5*(I - Ms)@J
    d = np.linalg.det(T)
    ph = np.unwrap(np.angle(d))
    return (ph[-1] - ph[0])/(2*np.pi), abs(d[-1]/d[0] - 1), np.min(np.abs(d))/np.max(np.abs(d))
if __name__ == '__main__':
    pts = {'Gold/Vakuum -11+1,2i': -11+1.2j, 'Gold/Wasser (-11+1,2i)/1,77': (-11+1.2j)/1.77,
           'Silber/Vakuum -15+0,5i': -15+0.5j, 'Resonanzbereich -3+0,3i': -3+0.3j, '-2+0,3i': -2+0.3j,
           '-1,5+0,2i': -1.5+0.2j, '-0,5+0,1i': -0.5+0.1j, 'Glas 2,25': 2.25, '-2+3,5i': -2+3.5j, '-5+0,5i': -5+0.5j}
    for name, k in pts.items():
        w, cl, mn = winding(k); we = winding(k, a=np.conj(k))[0]
        print(f"{name:30s} kappa={k:.3f}: Windung {w:+.3f} (Energiewahl {we:+.3f}), Schluss {cl:.1e}, min|det|/max {mn:.1e}")
    Re = np.linspace(-8, 0.5, 69); Im = np.linspace(-3, 3, 49)
    W = np.array([[winding(x + 1j*y)[0] for x in Re] for y in Im])
    json.dump({'Re': Re.tolist(), 'Im': Im.tolist(), 'W': W.tolist()}, open('corner_winding.json', 'w'))
    print('Windungszahlen im Gitter:', np.unique(np.round(W)))
