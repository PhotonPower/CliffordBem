"""Inkrementelle Torus-Studien mit Modenzerlegung (BoR). Speichert nach jedem Punkt."""
import sys, json, os, time, numpy as np
from bor import torus_bor, angles
g = torus_bor(); MS = [0, 1, 2, 3, 4, 5, 6, 16]
VAR = {'maxwell': dict(kind='maxwell'), 'dirac': dict(kind='dirac'),
       'dual_W1': dict(kind='dirac', dual=True, W='one'), 'dual_J': dict(kind='dirac', dual=True, W='J')}
def tasks(job):
    if job == 'lf':
        for mat, e1 in [('glas', 2.25), ('gold', -11+1.2j)]:
            for v in ['maxwell', 'dirac', 'dual_W1']:
                for om in [0.0, 1e-6, 1e-4, 1e-2, 0.1, 0.3, 1.0]:
                    yield f'{mat}_{v}|{om:g}', (mat, e1, v, om)
    if job == 'eps':
        for e in [-30,-25,-20,-16,-13,-11,-9,-7.5,-6,-5,-4,-3.3,-2.7,-2.2,-1.8,-1.5,-1.3,-1.15,-1.05,-0.95,-0.85,-0.7,-0.5,-0.3]:
            for v in ['maxwell', 'dirac', 'dual_W1', 'dual_J']:
                yield f'{v}|{e:g}', ('x', complex(e), v, 1e-3)
job = sys.argv[1]; budget = float(sys.argv[2]); T0 = time.time()
fn = f'bor_{job}.json'; res = json.load(open(fn)) if os.path.exists(fn) else {}
for key, (mat, e1, v, om) in tasks(job):
    if key in res: continue
    if time.time() - T0 > budget: break
    c = VAR[v]; r = angles(g, MS, c['kind'], om, e1, 1, 1.0, 1, dual=c.get('dual', False), W=c.get('W', 'J'))
    res[key] = {str(m): r[m] for m in MS}; json.dump(res, open(fn, 'w'))
print(job, len(res), '/', sum(1 for _ in tasks(job)))
