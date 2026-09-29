"""Auswertung der Zweitor-Formulierung (scatter_coated --twoport, results/twoport.csv, results/twoport_chiral.csv).

Bewertet werden Absolutwerte gegen Aden-Kerker bzw. die chirale Schicht-Mie-Loesung, weil das Zweitor E_2 auf der
Aussenflaeche (Radius 1 + d) rechnet: Differenzen zu einer Rechnung ohne Schicht auf Radius 1 waeren durch den
unterschiedlichen Diskretisierungsfehler verfaelscht. Zum Vergleich: Fehler der Kugel ohne Schicht auf demselben Netz.
Methoden: twoport (Zweitor), thin-dirac2fit (Duennschicht 2. Ordnung), zwei Flaechen (exakt, LayeredScatteringProblem).
Aufruf: python3 tools/analyze_twoport.py
"""
import csv, os, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mie_coated as mc, mie_chiral_layered as mcl

om = 0.5; core = complex(-11, 1.2)
if os.path.exists('results/twoport.csv'):
    rows = list(csv.DictReader(open('results/twoport.csv')))
    res = defaultdict(dict); bare = {}
    for r in rows:
        n = int(r['n']); c = r['coat']
        if c == 'none': bare[n] = r; continue
        if c.startswith('twoport:'): m, spec = 'Zweitor', c.split(':', 1)[1]
        elif c.startswith('thin'): m, spec = 'Dünnschicht', c.split(':', 1)[1]
        else: m, spec = 'zwei Flächen', c
        d = float(spec.split(',')[0]); res[(d, m)][n] = r
    S0 = mc.forward_amplitude(om, [1.0], [core])
    print("Kugel ohne Schicht: rel. Fehler sigma / Fehler arg S je Netz")
    for n in sorted(bare):
        s = float(bare[n]['sigma_ext']); print(f"  n={n:2d}: {s / (16 * np.pi * S0.real) - 1:+7.2%}  {float(bare[n]['S_arg']) - np.angle(S0):+.4f} rad")
    print(f"\n{'d':>6} {'Methode':>13} | " + " | ".join(f"n={n:2d}: sigma / arg S / It." for n in [4, 6, 8, 12]))
    for (d, m) in sorted(res, key=lambda k: (k[0], k[1])):
        S = mc.forward_amplitude(om, [1.0, 1.0 + d], [core, 2.25]); sx = 16 * np.pi * S.real
        cells = []
        for n in [4, 6, 8, 12]:
            r = res[(d, m)].get(n)
            cells.append(f"{float(r['sigma_ext']) / sx - 1:+6.2%} {float(r['S_arg']) - np.angle(S):+.4f} {int(r['iterations']):4d}" if r else " " * 24)
        print(f"{d:6.3f} {m:>13} | " + " | ".join(cells))

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
