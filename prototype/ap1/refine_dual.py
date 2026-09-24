from scipy.optimize import minimize_scalar
from mfs import *
for r in [(16,32,50),(20,40,70)]:
    s = Sphere(*r)
    f = lambda e: min_angle(s, 'dirac', 0.5, complex(e), 1, 1.0, 1, dual=True)
    o = minimize_scalar(f, bounds=(-0.78, -0.62), method='bounded', options={'xatol': 1e-5})
    print(r, f"eps*={o.x:.5f}, sin theta={o.fun:.2e}", flush=True)
