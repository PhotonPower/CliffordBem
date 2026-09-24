"""Konvergenz einer Ausgangsgroesse: Streufeld F^s(x) = -C h^s(x) an x = (5,3), h^s = h - h_inc,
Loesung per direktem Loeser auf geometrischen Netzen."""
import numpy as np
from cl3 import MV
from square_geo import build_geo
def far_field(ng, kappa, x=np.array([5.0, 3.0]), q=0.5):
    A, X, rr = build_geo(ng, 0.0, kappa, q=q); Lg = build_geo.last['L']; d = np.repeat(np.sqrt(Lg), 8)
    Ns = len(X); h_inc = np.zeros((Ns, 8), complex); h_inc[:, 1] = 1.0
    g = np.linalg.solve(A, d*h_inc.ravel()); h = (g/d).reshape(Ns, 8)
    # Normalen der Elemente (Seite k)
    Vv = np.array([(1, -1), (1, 1), (-1, 1), (-1, -1)], float)
    M = Ns//4; n = []
    for k in range(4):
        u = (Vv[(k+1) % 4] - Vv[k])/2; n += [np.array([u[1], -u[0]])]*M
    Fs = np.zeros(8, complex)
    for j in range(Ns):
        z = np.array([x[0]-X[j][0], x[1]-X[j][1], 0.]); K = z/(2*np.pi*(z@z))
        Fs += -(-(MV.vec(K)*MV.vec([n[j][0], n[j][1], 0])*MV(h[j]-h_inc[j])).c)*Lg[j]
    # Hilfsanteil der Loesung
    aux = np.sqrt(np.sum(Lg[:, None]*np.abs(h[:, [0, 7]])**2)/np.sum(Lg[:, None]*np.abs(h)**2))
    return Fs, aux
if __name__ == '__main__':
    for name, k in [('innen -3+0,3i', -3+0.3j), ('-4', -4+0j), ('Gold', -11+1.2j), ('Glas', 2.25+0j)]:
        prev = None; line = []
        for ng in [8, 16, 24, 32]:
            Fs, aux = far_field(ng, k)
            v = Fs[1]  # e1-Komponente
            line.append(f"ng={ng}: F1={v.real:+.6f}{v.imag:+.6f}i" + (f" (Diff {abs(v-prev):.1e})" if prev is not None else "") + f" Hilfs {aux:.0e}")
            prev = v
        print(name, '|', '; '.join(line), flush=True)
