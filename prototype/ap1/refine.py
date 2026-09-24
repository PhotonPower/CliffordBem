import numpy as np, json, sys
from scipy.optimize import minimize_scalar
from scipy.special import spherical_jn
from mfs import *
res = int(sys.argv[1]) if len(sys.argv)>1 else 0
s = [Sphere(16,32,50), Sphere(20,40,70)][res]
f = lambda om, kind='dirac', dual=False: min_angle(s, kind, om, 2.25, 1, 1.0, 1, dual=dual)
out = {}
for c in [1.94, 2.69, 3.52, 4.27]:
    r = minimize_scalar(f, bounds=(c-0.08, c+0.08), method='bounded', options={'xatol':1e-4})
    out[c] = (r.x, r.fun); print(f"Dirac lok. Min nahe {c}: om*={r.x:.4f}, sin(theta)={r.fun:.4f}", flush=True)
json.dump(out, open(f'refine_{res}.json','w'))
