"""Auswertung chiraler Schichten (scatter_coated --cd, results/thin_chiral.csv) gegen tools/mie_chiral_layered.py.

Je Rechnung (Netz n, Schichtangabe, Kern) werden die Zeilen beider Helizitaeten (Spalte pol = +1 / -1) gepaart:
sigma_+, sigma_-, CD = sigma_+ - sigma_- und die Schichtwirkung Delta sigma_+ = sigma_+ - sigma_+(ohne Schicht, gleiches Netz).
Die Chiralitaet des Kerns steht nicht in der CSV: --core-chi re=chi ordnet sie ueber den Realteil von eps_Kern zu
(Standard fuer die Studie: 2.25=0.2). Kugel mit Kernradius 1, Schichten nach aussen, Vakuum aussen.
Aufruf: python3 tools/analyze_chiral_thin.py [results/thin_chiral.csv] [--core-chi 2.25=0.2]
"""
import csv, sys, os
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mie_chiral_layered as mcl

args = [a for a in sys.argv[1:]]
fn = 'results/thin_chiral.csv'; cmap = {2.25: 0.2}
i = 0
while i < len(args):
    if args[i] == '--core-chi': k, v = args[i + 1].split('='); cmap[float(k)] = float(v); i += 2
    else: fn = args[i]; i += 1
rows = list(csv.DictReader(open(fn)))
groups = defaultdict(dict)
for r in rows:
    key = (int(r['n']), r['coat'], float(r['core_re']), float(r['core_im']))
    groups[key][int(r['pol'])] = r
om = float(rows[0]['omega'])


def mie(core, layers):
    radii = [1.0]; media = [core]
    for d, e, chi in layers: radii.append(radii[-1] + d); media.append((e, 1.0, chi))
    sp = mcl.cross_sections(om, radii, media, +1)[0]; sm = mcl.cross_sections(om, radii, media, -1)[0]
    return sp, sm


print(f"{'n':>3} {'Kern':>14} {'Schicht':>22} {'Verf.':>6} | {'sigma+':>10} {'Mie':>10} | {'CD':>11} {'CD Mie':>11} {'Fehler':>7} | {'dSig+':>9} {'Mie':>9} {'Fehler':>7}")
for (n, coat, cre, cim), g in sorted(groups.items(), key=lambda t: (t[0][2], t[0][1], t[0][0])):
    if coat == 'none' or +1 not in g or -1 not in g: continue
    chic = cmap.get(cre, 0.0); core = (complex(cre, cim), 1.0, chic)
    spec = coat.split(':', 1)[1] if ':' in coat else coat
    method = 'Dünn' if coat.startswith('thin') else 'zwei'
    layers = []
    for part in spec.split(';'):
        v = [float(x) for x in part.split(',')]; layers.append((v[0], complex(v[1], v[2] if len(v) > 2 else 0.0), v[3] if len(v) > 3 else 0.0))
    sp, sm = float(g[+1]['sigma_ext']), float(g[-1]['sigma_ext'])
    mp, mm = mie(core, layers); m0 = mie(core, [])[0]
    b = groups.get((n, 'none', cre, cim), {})
    dsp = sp - float(b[+1]['sigma_ext']) if +1 in b else float('nan')
    cd, cdm = sp - sm, mp - mm
    print(f"{n:3d} {str(core[0]) + ' chi ' + str(chic):>14} {spec:>22} {method:>6} | {sp:10.6f} {mp:10.6f} | {cd:+11.4e} {cdm:+11.4e} {cd / cdm - 1:+7.1%} |"
          f" {dsp:+9.5f} {mp - m0:+9.5f} {dsp / (mp - m0) - 1:+7.1%}")
