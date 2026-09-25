"""Mie-Spektrum einer Kugel mit tabellierten optischen Konstanten (gleiche lineare Interpolation in n, k wie
der Kern) und Vergleich mit einer spectrum-CSV.  python3 mie_spectrum.py ../results/au40_water.csv Au"""
import csv, sys, numpy as np
import mie
def load_nk(path):
    l, n, k = [], [], []; data = False
    for line in open(path):
        if 'data: |' in line: data = True; continue
        if data:
            p = line.split()
            if len(p) == 3: l.append(1000*float(p[0])); n.append(float(p[1])); k.append(float(p[2]))
            elif line.strip(): data = False
    return np.array(l), np.array(n), np.array(k)
def eps_tab(path, lam):
    l, n, k = load_nk(path); return (np.interp(lam, l, n) + 1j*np.interp(lam, l, k))**2
if __name__ == '__main__':
    fn, mat = sys.argv[1], sys.argv[2]
    rows = list(csv.DictReader(open(fn)))
    for r in rows:
        L = float(r['lambda_nm']); a = float(r['unit_nm']); nb = float(r['nbg'])
        e = eps_tab(f'../data/materials/{mat}_Johnson.yml', L); om = 2*np.pi*a/L
        sM = mie.qext(om, e, nb**2)*np.pi*a*a; sB = float(r['sigma_nm2'])
        print(f"  lambda {L:6.1f} nm (N={r['N']}): BEM {sB:10.1f} nm^2, Mie {sM:10.1f} nm^2, rel. Abw. {abs(sB-sM)/sM:.2e}")
