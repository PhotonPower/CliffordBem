import numpy as np
from scipy.special import spherical_jn
from scipy.optimize import brentq
from mfs import *
import sys
res = int(sys.argv[1]); s = [Sphere(16,32,50), Sphere(20,40,70), Sphere(24,48,100)][res]
# Dirichlet- und Neumann-Eigenwerte der Kugel (k R), l = 0,1,2
cand = []
for l in range(3):
    g = lambda z: spherical_jn(l, z); gd = lambda z: spherical_jn(l, z, derivative=True)
    zs = np.linspace(0.5, 7, 400)
    for fn, nm in [(g,'D'),(gd,'N')]:
        v = fn(zs)
        for i in np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]:
            cand.append((brentq(fn, zs[i], zs[i+1]), f'{nm},l={l}'))
out = []
for z, nm in sorted(cand):
    for which, om in [('k1', z/1.5), ('k2', z)]:
        if om > 4.5: continue
        out.append((om, nm, which, min_angle(s, 'dirac', om, 2.25, 1, 1.0, 1), min_angle(s, 'dirac', om, 2.25, 1, 1.0, 1, dual=True)))
for o in sorted(out): print(f"om={o[0]:.4f} ({o[1]}, {o[2]} R): sin theta = {o[3]:.4f}, dual {o[4]:.4f}")
