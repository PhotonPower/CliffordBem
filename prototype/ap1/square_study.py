import sys, json, os, time, numpy as np
from square_bem import square_mesh, E_matrix, T1_matrix, weighted_svals
Ns = [int(x) for x in sys.argv[1].split(',')]; fn = 'square_study.json'
res = json.load(open(fn)) if os.path.exists(fn) else {}
cases = {'innen -3+0,3i': -3+0.3j, 'aussen Gold -11+1,2i': -11+1.2j, 'Glas 2,25': 2.25+0j, 'verlustfrei -4': -4+0j}
for N in Ns:
    P0, P1, Nr, V = square_mesh(N); E, X, L = E_matrix(P0, P1, Nr)
    for name, k in cases.items():
        T = T1_matrix(E, Nr, k)
        for beta in [0.0, 0.35]:
            key = f'{name}|{beta}|{N}'
            if key in res: continue
            t = time.time(); s, smax = weighted_svals(T, X, L, V, beta)
            res[key] = [float(x) for x in s] + [float(smax)]; json.dump(res, open(fn, 'w'))
            print(f"N={N:3d} {name:22s} beta={beta}: kleinste SW {np.round(s[:6],4)}  max {smax:.2f}  ({time.time()-t:.0f}s)", flush=True)
