import numpy as np, json, sys, time
from pml2d import far_field
kappas = {'Gold': -11+1.2j, '-2,8+0,3i': -2.8+0.3j, '-2,8+0,05i': -2.8+0.05j, '-2 (verlustfrei)': -2.0+0j}
res = {}
for name, k in kappas.items():
    for th in [0.0, -0.5, -1.0]:
        vals = []
        for ng in [8, 16, 24, 32]:
            F, smin, smax = far_field(ng, k, th); vals.append(F[1])
        d = [abs(vals[i+1]-vals[i]) for i in range(3)]
        print(f"{name:16s} theta={th:+.1f}: " + ' '.join(f"{v.real:+.6f}{v.imag:+.6f}i" for v in vals) + "  | Diff " + ' '.join(f"{x:.1e}" for x in d), flush=True)
        res[f"{name}|{th}"] = [[v.real, v.imag] for v in vals]
json.dump(res, open('study1.json', 'w'))
