"""Absorbierender Eckabschluss: kuenstlicher Verlust in den innersten L_abs Schichten,
kappa_loc = kappa + i*gamma*(l/L_abs)^2, l = Schichten unterhalb r_abs = h_min*2^L_abs."""
import numpy as np, json, sys
from square_enrich import far_field
def make_fn(ng, L_abs, gamma, q=0.5):
    hmin = 0.5*q**ng; rabs = hmin*(1/q)**L_abs
    def fn(r, kappa):
        if r >= rabs: return kappa
        l = np.log(rabs/r)/np.log(1/q)
        return kappa + 1j*gamma*min(l/L_abs, 1.0)**2
    return fn
if __name__ == '__main__':
    k = complex(sys.argv[1]) if len(sys.argv) > 1 else -2.8+0.3j
    out = {}
    for L_abs, gamma in [(0, 0.0), (6, 2.0), (6, 5.0), (10, 5.0)]:
        vals = []
        for ng in range(12, 33, 4):
            fn = None if L_abs == 0 else make_fn(ng, L_abs, gamma)
            v = far_field(ng, k, kappa_fn=fn)[1]; vals.append(v)
        d = [abs(vals[i+1]-vals[i]) for i in range(len(vals)-1)]
        out[f'{L_abs}|{gamma}'] = [[v.real, v.imag] for v in vals]
        print(f"L_abs={L_abs:2d} gamma={gamma}: F1(ng=12..32) " + ', '.join(f'{v.real:.4f}{v.imag:+.4f}i' for v in vals) + "  Diffs " + ' '.join(f'{x:.1e}' for x in d), flush=True)
    json.dump(out, open(f'absorber_{k}.json', 'w'))
