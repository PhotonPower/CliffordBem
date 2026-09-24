import json, numpy as np, sys
def load(job):
    r = json.load(open(f'sweep_{job}.json')); d = {}
    for k, v in r.items():
        n, p = k.split('|'); d.setdefault(n, []).append((float(p), v))
    return {n: (np.array(sorted(l))[:,0], np.array(sorted(l))[:,1]) for n, l in d.items()}
if __name__ == '__main__':
    for job in sys.argv[1:]:
        print('==', job)
        for n, (x, v) in load(job).items():
            loc = [(round(x[j],3), round(v[j],4)) for j in range(1,len(v)-1) if v[j]<v[j-1] and v[j]<v[j+1]]
            print(f"{n:16s} min {v.min():.4f} @ {x[np.argmin(v)]:.3g}  max {v.max():.3f}  lokMin {loc[:6]}")
