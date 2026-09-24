"""Inkrementelle Parameterstudien; jeder Aufruf rechnet bis zum Zeitbudget weiter und
speichert nach jedem Punkt (sweep_<job>.json)."""
import sys, json, os, time, numpy as np
from mfs import *
s = Sphere(16, 32, 50)
if sys.argv[1].startswith('torus'):
    from torus import Torus
    s = Torus(nu=48, nv=24, su=18, sv=8)
T0 = time.time(); budget = float(sys.argv[2]) if len(sys.argv) > 2 else 250
def tasks(job):
    if job == 'freq':
        for name, c in [('maxwell', dict(kind='maxwell')), ('dirac', dict(kind='dirac')),
                        ('dirac_dual', dict(kind='dirac', dual=True)), ('dirac_b1', dict(kind='dirac', b=1.0)),
                        ('dirac_bad', dict(kind='dirac', b=-1.5)), ('dirac_dual_W1', dict(kind='dirac', dual=True, W='one')), ('dirac_energy', dict(kind='dirac', a='energy'))]:
            for om in np.linspace(0.05, 4.5, 60):
                yield (name, float(om)), (lambda om=om, c=c: min_angle(s, c['kind'], om, 2.25, 1, 1.0, 1, a=c.get('a'), b=c.get('b'), dual=c.get('dual', False), W=c.get('W','J')))
    if job == 'eps':
        for name, kind, dual, W in [('maxwell','maxwell',False,'J'), ('dirac','dirac',False,'J'), ('dirac_dual','dirac',True,'J'),
                                    ('dirac_dual_W1','dirac',True,'one'), ('maxwell_dual','maxwell',True,'J')]:
            for e in np.linspace(-6, -0.2, 59):
                yield (name, float(e)), (lambda e=e, kind=kind, dual=dual, W=W: min_angle(s, kind, 0.5, complex(e), 1, 1.0, 1, dual=dual, W=W))
        for e in np.linspace(-6, -0.2, 59):
            yield ('dirac_energy', float(e)), (lambda e=e: min_angle(s, 'dirac', 0.5, complex(e), 1, 1.0, 1, a='energy'))
        for name, W in [('lossy_dual', 'J'), ('lossy_dual_W1', 'one')]:
            for e in np.linspace(-6, -0.2, 59):
                yield (name, float(e)), (lambda e=e, W=W: min_angle(s, 'dirac', 0.5, complex(e, 0.1), 1, 1.0, 1, dual=True, W=W))
    if job == 'lowfreq':
        for mat, e1 in [('glas', 2.25), ('gold', -11+1.2j)]:
            for name, kind, dual in [('dirac','dirac',False), ('dirac_dual','dirac',True), ('maxwell','maxwell',False)]:
                for om in np.logspace(-3, 0, 13):
                    yield (mat+'_'+name, float(om)), (lambda om=om, e1=e1, kind=kind, dual=dual: min_angle(s, kind, om, e1, 1, 1.0, 1, dual=dual))
    if job == 'torus_lf':
        for name, kind, dual, W in [('maxwell','maxwell',False,'J'), ('dirac','dirac',False,'J'), ('dirac_dual_W1','dirac',True,'one')]:
            for om in [1e-3, 1e-2, 1e-1, 0.5, 1.0]:
                yield ('glas_'+name, om), (lambda om=om, kind=kind, dual=dual, W=W: min_angle(s, kind, om, 2.25, 1, 1.0, 1, dual=dual, W=W))
    if job == 'torus_eps':
        for name, kind in [('maxwell','maxwell'), ('dirac','dirac')]:
            for e in [-2.0, -4.0, -6.0, -8.0, -11.0, -14.0, -20.0, -30.0]:
                yield (name, e), (lambda e=e, kind=kind: min_angle(s, kind, 1e-3, complex(e, 1.2), 1, 1.0, 1))
job = sys.argv[1]; fn = f'sweep_{job}.json'
res = json.load(open(fn)) if os.path.exists(fn) else {}
n_new = 0
for (name, p), f in tasks(job):
    key = f'{name}|{p:.6g}'
    if key in res: continue
    if time.time() - T0 > budget: break
    res[key] = f(); n_new += 1
    json.dump(res, open(fn, 'w'))
total = sum(1 for _ in tasks(job))
print(f'{job}: {len(res)}/{total} fertig (+{n_new})')
