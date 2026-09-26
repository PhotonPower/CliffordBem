import numpy as np, sys
import pml2d_galerkin as pg
from absorber2d import panel_r, profile, solve
kappa = complex(sys.argv[1]); q = float(sys.argv[2])
cases = [tuple(float(v) for v in c.split(':')) for c in sys.argv[4].split(',')]   # r_abs:dekaden:dmax
for depth in [float(v) for v in sys.argv[3].split(',')]:
    ng = int(np.ceil(np.log(depth/0.5)/np.log(q)))
    parts = pg.build(ng, kappa, 0.0, q=q, parts_only=True); r = panel_r(parts[2])
    out = [f"ohne {solve(parts, [kappa]*len(r)):.6f}"]
    for r_abs, dec, dmax in cases:
        out.append(f"({r_abs:.0e},{dec:g},{dmax:g}) {solve(parts, profile(r, kappa, r_abs, dec*np.log(10), dmax)):.6f}")
    print(f"Tiefe {depth:.0e}: " + '; '.join(out), flush=True)
