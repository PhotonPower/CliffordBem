import numpy as np, sys
from torus_sym import TorusSym, angles_sym
for p in [dict(nt=48, L=3, n=40, mt=12), dict(nt=48, L=3, n=40, mt=16), dict(nt=48, L=3, n=40, mt=20),
          dict(nt=40, L=2, n=40, mt=16), dict(nt=64, L=4, n=40, mt=16)]:
    g = TorusSym(**p)
    r = angles_sym(g, 'maxwell', 0.5, 2.25, 1, 1.0, 1, nang=200)
    per = {}
    for v, m in r: per.setdefault(m, v)
    print(p, 'min je Mode m=0..3:', [round(per.get(m, np.nan), 4) for m in range(4)], flush=True)
