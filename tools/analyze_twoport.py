"""Auswertung der Zweitor-Formulierung (scatter_coated --twoport, results/twoport.csv, results/twoport_chiral.csv).

Bewertet werden Absolutwerte gegen Aden-Kerker bzw. die chirale Schicht-Mie-Loesung, weil das Zweitor E_2 auf der
Aussenflaeche (Radius 1 + d) rechnet: Differenzen zu einer Rechnung ohne Schicht auf Radius 1 waeren durch den
unterschiedlichen Diskretisierungsfehler verfaelscht. Zum Vergleich: Fehler der Kugel ohne Schicht auf demselben Netz.
Methoden: twoport (Zweitor), thin-dirac2fit (Duennschicht 2. Ordnung), zwei Flaechen (exakt, LayeredScatteringProblem).
Aufruf: python3 tools/analyze_twoport.py [results/twoport.csv]   (v0.20: results/twoport2.csv, mit Unterteilung dicker
Schichten und Mehrfachschichten; Schichtangabe "d,re,im;d,re,im" von innen nach aussen)
"""
import csv, os, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mie_coated as mc, mie_chiral_layered as mcl

om = 0.5; core = complex(-11, 1.2)
fn = sys.argv[1] if len(sys.argv) > 1 else 'results/twoport.csv'


def mie_sigma_arg(spec):
    radii = [1.0]; eps = [core]
    for part in spec.split(';'):
        v = [float(x) for x in part.split(',')]; radii.append(radii[-1] + v[0]); eps.append(complex(v[1], v[2] if len(v) > 2 else 0.0))
    S = mc.forward_amplitude(om, radii, eps); return 16 * np.pi * S.real, np.angle(S)


if os.path.exists(fn):
    rows = list(csv.DictReader(open(fn)))
    res = defaultdict(dict); bare = {}
    for r in rows:
        n = int(r['n']); c = r['coat']
        if c == 'none': bare[n] = r; continue
        if c.startswith('twoport:'): m, spec = 'Zweitor', c.split(':', 1)[1]
        elif c.startswith('thin'): m, spec = 'Dünnschicht', c.split(':', 1)[1]
        else: m, spec = 'zwei Flächen', c
        res[(spec, m)][n] = r
    S0 = mc.forward_amplitude(om, [1.0], [core])
    print("Kugel ohne Schicht: rel. Fehler sigma / Fehler arg S je Netz")
    for n in sorted(bare):
        s = float(bare[n]['sigma_ext']); print(f"  n={n:2d}: {s / (16 * np.pi * S0.real) - 1:+7.2%}  {float(bare[n]['S_arg']) - np.angle(S0):+.4f} rad")
    print(f"\n{'d':>6} {'Methode':>13} | " + " | ".join(f"n={n:2d}: sigma / arg S / It." for n in [4, 6, 8, 12]))
    for (spec, m) in sorted(res, key=lambda k: (';' in k[0], float(k[0].split(',')[0]), k[0], k[1])):
        sx, ax = mie_sigma_arg(spec)
        cells = []
        for n in [4, 6, 8, 12]:
            r = res[(spec, m)].get(n)
            cells.append(f"{float(r['sigma_ext']) / sx - 1:+6.2%} {float(r['S_arg']) - ax:+.4f} {int(r['iterations']):4d}" if r else " " * 24)
        print(f"{spec:>30} {m:>13} | " + " | ".join(cells))

if os.path.exists('results/twoport_chiral.csv'):
    g = defaultdict(dict)
    for r in csv.DictReader(open('results/twoport_chiral.csv')):
        g[(int(r['n']), r['coat'])][int(r['pol'])] = float(r['sigma_ext'])
    gold = (core, 1.0, 0.0)
    thin = defaultdict(dict)
    if os.path.exists('results/thin_chiral.csv'):
        for r in csv.DictReader(open('results/thin_chiral.csv')):
            if float(r['core_re']) == -11.0: thin[(int(r['n']), r['coat'])][int(r['pol'])] = float(r['sigma_ext'])
    print("\nchirale Huelle (chi = 0,1): rel. Fehler des CD gegen Mie -- Zweitor / Duennschicht 2. Ordnung / zwei Flaechen")
    for d in [0.02, 0.05, 0.1, 0.2]:
        cdm = mcl.cross_sections(om, [1, 1 + d], [gold, (2.25, 1.0, 0.1)], +1)[0] - mcl.cross_sections(om, [1, 1 + d], [gold, (2.25, 1.0, 0.1)], -1)[0]
        out = []
        for n in [6, 8, 12]:
            def e(v): return f"{(v[+1] - v[-1]) / cdm - 1:+6.1%}" if v and +1 in v and -1 in v else "   –  "
            out.append(f"n={n:2d}: {e(g.get((n, f'twoport:{d},2.25,0,0.1')))} / {e(thin.get((n, f'thin-dirac2fit:{d},2.25,0,0.1')))} / {e(thin.get((n, f'{d},2.25,0,0.1')))}")
        print(f"  d={d:<5} CD Mie {cdm:+.4e} | " + " | ".join(out))
