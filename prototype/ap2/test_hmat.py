import numpy as np, sys
from mesh import icosphere, Mesh
from assemble import assemble_E
from hmat import HMatrix
n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
m = Mesh(*icosphere(n)); N = len(m.T); k = 1.5
E = assemble_E(m, k); rng = np.random.default_rng(0); x = rng.normal(size=8*N) + 1j*rng.normal(size=8*N); y0 = E@x
for mode in ['comp', 'joint']:
    for eps in [1e-3, 1e-5]:
        H = HMatrix(m, k, eps=eps, mode=mode)
        y = H.matvec(x); print(f"     rel. Fehler Matrix-Vektor-Produkt: {np.linalg.norm(y - y0)/np.linalg.norm(y0):.1e}", flush=True)
