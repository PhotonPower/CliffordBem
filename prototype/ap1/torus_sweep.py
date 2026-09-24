import sys, json, os, time, numpy as np
from torus_sym import TorusSym, angles_sym
g = TorusSym(nt=48, L=3, n=40, mt=16)
T0 = time.time(); budget = float(sys.argv[1]) if len(sys.argv) > 1 else 250
fn = 'sweep_torus.json'; res = json.load(open(fn)) if os.path.exists(fn) else {}
cfgs = [('dirac', dict(kind='dirac')), ('dirac_energy', dict(kind='dirac', a='energy')),
        ('maxwell', dict(kind='maxwell')), ('dual_W1', dict(kind='dirac', dual=True, W='one'))]
tasks = [(mat, e1, n, c, om) for mat, e1 in [('glas', 2.25), ('gold', -11+1.2j)] for n, c in cfgs for om in np.logspace(-3, 0, 10)]
for mat, e1, n, c, om in tasks:
    key = f'{mat}_{n}|{om:.6g}'
    if key in res: continue
    if time.time() - T0 > budget: break
    r = angles_sym(g, c['kind'], om, e1, 1, 1.0, 1, a=c.get('a'), dual=c.get('dual', False), W=c.get('W', 'J'), nang=3)
    res[key] = [[float(v), int(m)] for v, m in r]
    json.dump(res, open(fn, 'w'))
print(f'torus: {len(res)}/{len(tasks)}')
