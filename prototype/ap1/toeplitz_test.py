"""Pruefung der Endlicher-Abschnitt-Hypothese: kleinster |Eigenwert| der Systemmatrix und des
isolierten Eckblocks ueber n_g, verglichen mit den Ausschlaegen von F_1(n_g)."""
import numpy as np, json
from square_geo import build_geo
from square_precond import corner_blocks
k = -2.8+0.3j; F = json.load(open('square_ng_scan.json')); out = {}
for ng in range(6, 31, 2):
    A, X, rr = build_geo(ng, 0.0, k)
    ev = np.linalg.eigvals(A); lam = ev[np.argmin(np.abs(ev))]
    b = corner_blocks(X, rr, 0.5)[0]; evb = np.linalg.eigvals(A[np.ix_(b, b)]); lb = evb[np.argmin(np.abs(evb))]
    f = complex(*F[str(ng)])
    out[ng] = dict(lam=[lam.real, lam.imag], lam_block=[lb.real, lb.imag], F=[f.real, f.imag])
    print(f"n_g={ng:2d}: min|EW| gesamt {abs(lam):.4f} ({lam:.3f}), Eckblock {abs(lb):.4f}   |F1|={abs(f):.4f}", flush=True)
json.dump(out, open('toeplitz_test.json', 'w'))
