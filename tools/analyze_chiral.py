"""Chirale Kugel: BEM gegen chirale Mie-Loesung (Q_+, Q_-, CD), mit Extrapolation Q(n) = Q_inf - C n^-p."""
import csv, sys
import numpy as np
from scipy.optimize import brentq
import mie_chiral as mc
def extrap(ns, qs):
    rat = (qs[2]-qs[1])/(qs[1]-qs[0])
    f = lambda p: (ns[1]**-p - ns[2]**-p)/(ns[0]**-p - ns[1]**-p) - rat
    try:
        p = brentq(f, 0.2, 6); C = (qs[2]-qs[1])/(ns[1]**-p - ns[2]**-p); return p, qs[2] + C*ns[2]**-p
    except ValueError: return float('nan'), float('nan')
for fn in sys.argv[1:]:
    rows = sorted(csv.DictReader(open(fn)), key=lambda r: int(r['n']))
    r0 = rows[0]; om = float(r0['omega']); eps = complex(float(r0['eps1_re']), float(r0['eps1_im'])); chi = complex(float(r0['chi_re']), float(r0['chi_im']))
    Mp = mc.cross_sections(om, eps, chi=chi, s=+1)[0]/np.pi; Mm = mc.cross_sections(om, eps, chi=chi, s=-1)[0]/np.pi
    print(f"{fn}: omega={om}, eps={eps}, chi={chi};  Mie: Q_+ = {Mp:.6f}, Q_- = {Mm:.6f}, CD = {Mp-Mm:.6f}")
    for r in rows:
        qp, qm = float(r['Qplus']), float(r['Qminus'])
        print(f"   N={int(r['N']):6d}: Q_+ {qp:.6f} ({abs(qp-Mp)/Mp:.1e})  Q_- {qm:.6f} ({abs(qm-Mm)/Mm:.1e})  CD {qp-qm:.6f} ({abs(qp-qm-(Mp-Mm))/abs(Mp-Mm):.1e})  It. {r['its_plus']}/{r['its_minus']}")
    ns = [int(r['n']) for r in rows[-3:]]
    for name, key, M in [('Q_+', 'Qplus', Mp), ('Q_-', 'Qminus', Mm), ('CD', 'CD', Mp - Mm)]:
        p, qi = extrap(ns, [float(r[key]) for r in rows[-3:]])
        print(f"   Extrapolation {name} aus n={ns}: p = {p:.2f}, {qi:.6f}, rel. Abw. zu Mie {abs(qi-M)/abs(M):.1e}")
