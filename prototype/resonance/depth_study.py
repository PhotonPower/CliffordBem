import numpy as np, sys
from pml2d_galerkin import far_field
k = complex(sys.argv[1]); q = float(sys.argv[2]); ths = [float(v) for v in sys.argv[3].split(',')]; depths = [float(v) for v in sys.argv[4].split(',')]
for th in ths:
    line = []
    for dep in depths:
        ng = int(np.ceil(np.log(dep/0.5)/np.log(q))); F = far_field(ng, k, th, q=q)
        line.append(f"Tiefe {dep:.0e} (ng={ng}): {F[1].real:+.7f}{F[1].imag:+.7f}i")
    print(f"Galerkin kappa={k}, q={q}, theta={th:+.1f}: " + '; '.join(line), flush=True)
