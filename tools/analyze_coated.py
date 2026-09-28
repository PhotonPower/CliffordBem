"""Auswertung von scatter_coated-CSV-Dateien (Kugel, Kernradius 1, Schichten nach aussen) gegen Aden-Kerker.

Fuer jede Netzfeinheit n und Schichtdicke: sigma_ext und arg S(0) der beschichteten, der neutralen (Schicht aus
Aussenmedium) und der unbeschichteten Rechnung sowie die Schichtwirkung als Differenz
  Delta_bare    = beschichtet - ohne Schicht,
  Delta_neutral = beschichtet - neutral           (empfohlen, siehe docs/results_coated.md),
jeweils gegen die exakte Differenz aus tools/mie_coated.py.
Aufruf: python3 tools/analyze_coated.py results/coated_thin.csv [--eps2 1.0]
"""
import csv, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import mie_coated as mc

fn = sys.argv[1]; eps2 = float(sys.argv[sys.argv.index('--eps2') + 1]) if '--eps2' in sys.argv else 1.0
rows = list(csv.DictReader(open(fn)))
data = defaultdict(dict)                       # (n, coat) -> {'coated'/'neutral': row};  (n, None) -> bare
coats = []
last = None
for r in rows:
    n = int(r['n']); c = r['coat']
    if c == 'none': data[(n, None)]['bare'] = r; continue
    if c == 'neutral': data[(n, last)]['neutral'] = r; continue
    data[(n, c)]['coated'] = r; last = c
    if c not in coats: coats.append(c)
ns = sorted({k[0] for k in data})
for c in coats:
    d, er, ei = [float(v) for v in c.split(',')]
    om = float(next(r for r in rows if r['coat'] == c)['omega']); core = complex(float(rows[0]['core_re']), float(rows[0]['core_im']))
    S_c = mc.forward_amplitude(om, [1.0, 1.0 + d], [core, complex(er, ei)], eps2); S_0 = mc.forward_amplitude(om, [1.0], [core], eps2)
    k2 = om ** 2 * eps2; sig = lambda S: 4 * np.pi / k2 * S.real
    dsig = sig(S_c) - sig(S_0); dph = np.angle(S_c) - np.angle(S_0)
    print(f"\nSchicht d = {d} (eps {er}{ei:+}i), omega = {om}: Mie sigma {sig(S_0):.5f} -> {sig(S_c):.5f} (Delta {dsig:+.5f}), "
          f"arg S {np.angle(S_0):.5f} -> {np.angle(S_c):.5f} (Delta {dph:+.5f})")
    print(f"{'n':>4} {'d/h':>6} | {'sigma':>9} {'neutral':>9} {'ohne':>9} | {'Delta_bare':>11} {'Delta_neut':>11} {'Fehler':>8} | {'dArg_bare':>10} {'dArg_neut':>10} {'Fehler':>8}")
    for n in ns:
        e = data.get((n, c), {}); b = data.get((n, None), {}).get('bare')
        if 'coated' not in e: continue
        h = 1.05 / n                                           # Kantenlaenge des Ikosaedernetzes (Radius 1)
        sc = float(e['coated']['sigma_ext']); ac = float(e['coated']['S_arg'])
        sn = float(e['neutral']['sigma_ext']) if 'neutral' in e else np.nan; an = float(e['neutral']['S_arg']) if 'neutral' in e else np.nan
        sb = float(b['sigma_ext']) if b else np.nan; ab = float(b['S_arg']) if b else np.nan
        print(f"{n:4d} {d / h:6.2f} | {sc:9.5f} {sn:9.5f} {sb:9.5f} | {sc - sb:+11.5f} {sc - sn:+11.5f} {((sc - sn) - dsig) / dsig:+8.1%} |"
              f" {ac - ab:+10.5f} {ac - an:+10.5f} {((ac - an) - dph) / dph:+8.1%}")
