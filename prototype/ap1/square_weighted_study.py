import sys, json, os, time, numpy as np
from square_weighted import build
fn = 'square_weighted.json'; res = json.load(open(fn)) if os.path.exists(fn) else {}
Ns = [int(x) for x in sys.argv[1].split(',')]; betas = [float(x) for x in sys.argv[2].split(',')]
cases = {'innen -3+0,3i': -3+0.3j, 'verlustfrei -4': -4+0j, 'Gold -11+1,2i': -11+1.2j, 'verlustfrei -2': -2+0j}
for N in Ns:
    for beta in betas:
        for name, k in cases.items():
            key = f'{name}|{beta}|{N}'
            if key in res: continue
            t = time.time(); A = build(N, beta, k)
            s = np.linalg.svd(A, compute_uv=False); ev = np.sort(np.abs(np.linalg.eigvals(A)))
            res[key] = dict(s=s[-6:][::-1].tolist(), ev=ev[:6].tolist()); json.dump(res, open(fn, 'w'))
            print(f"N={N:3d} beta={beta:.2f} {name:15s} SW {np.round(s[-6:][::-1][:5],4)}  |EW| {np.round(ev[:5],4)} ({time.time()-t:.0f}s)", flush=True)
