"""Auswertung der Kugel-Streurechnungen: Konvergenz gegen Mie und Extrapolation Q(n) = Q_inf - C n^-p."""
import csv, sys, math
from scipy.optimize import brentq
import mie
for fn in sys.argv[1:]:
    rows = sorted(csv.DictReader(open(fn)), key=lambda r: int(r['n']))
    om = float(rows[0]['omega']); e1 = complex(float(rows[0]['eps1_re']), float(rows[0]['eps1_im'])); M = mie.qext(om, e1)
    print(f"{fn}: omega = {om}, eps1 = {e1}, Mie Q_ext = {M:.6f}")
    for r in rows:
        Q = float(r['Qext']); print(f"   n={int(r['n']):3d} N={int(r['N']):6d}: Q = {Q:.6f}  rel. Abw. {abs(Q-M)/M:.2e}  It. {r['iterations']}  H {float(r['MB_H']):.0f} MB")
    for tri in zip(rows, rows[1:], rows[2:]):
        n = [int(r['n']) for r in tri]; q = [float(r['Qext']) for r in tri]
        rat = (q[2]-q[1])/(q[1]-q[0])
        f = lambda p: (n[1]**-p - n[2]**-p)/(n[0]**-p - n[1]**-p) - rat
        try:
            p = brentq(f, 0.2, 6); C = (q[2]-q[1])/(n[1]**-p - n[2]**-p); Qi = q[2] + C*n[2]**-p
            print(f"   Extrapolation aus n={n}: p = {p:.2f}, Q_inf = {Qi:.6f}, rel. Abw. zu Mie {abs(Qi-M)/M:.1e}")
        except ValueError:
            print(f"   Extrapolation aus n={n}: keine monotone Folge")
