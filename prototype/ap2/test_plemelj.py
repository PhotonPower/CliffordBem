import numpy as np
from mesh import icosphere, Mesh
from assemble import assemble_E, project
from fields import dirac_kernel, gp
rng = np.random.default_rng(1); c = rng.normal(size=8) + 1j*rng.normal(size=8)
for k in [0.0, 1.5, 1.2 + 0.4j]:
    for n in [2, 3, 4]:
        m = Mesh(*icosphere(n)); E = assemble_E(m, k, verbose=False)
        fi = lambda x: gp(dirac_kernel(x, np.array([2.3, -0.8, 0.9]), k), c)      # innere Loesung
        fe = lambda x: gp(dirac_kernel(x, np.array([0.1, 0.25, -0.2]), k), c)     # aeussere Loesung
        gi, ge = project(m, fi), project(m, fe)
        ei = np.linalg.norm(E@gi - gi)/np.linalg.norm(gi); ee = np.linalg.norm(E@ge + ge)/np.linalg.norm(ge)
        s = np.linalg.svd(E@E - np.eye(len(E)), compute_uv=False)
        print(f"k={k}, N={len(m.T)}: |E h_i - h_i|/|h_i| = {ei:.2e}, |E h_e + h_e|/|h_e| = {ee:.2e}", flush=True)
