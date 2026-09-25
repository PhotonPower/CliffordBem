import numpy as np, sys
from pml2d import far_field
k = complex(sys.argv[1]); ths = [float(x) for x in sys.argv[2].split(',')]
for th in ths:
    line = []
    for q in [float(x) for x in sys.argv[3].split(",")]:
        ng = int(np.ceil(np.log(float(sys.argv[4]))/np.log(q)))
        F, smin, smax = far_field(ng, k, th, q=q)
        line.append(f"q={q} (ng={ng}): {F[1].real:+.7f}{F[1].imag:+.7f}i")
    print(f"kappa={k}, theta={th:+.1f}: " + '; '.join(line), flush=True)
