import time, numpy as np
from torus_sym import TorusSym, angles_sym
for p in [dict(nt=40, L=3, n=32, mt=12), dict(nt=48, L=3, n=40, mt=16), dict(nt=64, L=4, n=48, mt=20)]:
    g = TorusSym(**p); t = time.time()
    out = []
    for kind in ['dirac', 'maxwell']:
        for e1 in [2.25, -11+1.2j]:
            v, m = angles_sym(g, kind, 0.5, e1, 1, 1.0, 1, nang=1)[0]; out.append(f"{kind[:1]}/{'G' if e1==2.25 else 'Au'}: {v:.4f}(m={m})")
    v, m = angles_sym(g, 'dirac', 0.5, 2.25, 1, 1.0, 1, dual=True, W='one', nang=1)[0]; out.append(f"dualW1/G: {v:.4f}(m={m})")
    print(p, ' | '.join(out), f"{time.time()-t:.0f}s", flush=True)
