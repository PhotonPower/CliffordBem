"""Auswertung von curved_galerkin.py: Extrapolation m -> oo je (Geometrie, Dichte, grobes Netz) und Tabelle der Fehler gegen
Mie. Modell Q(m) = Q_oo + c/m^2 (Restfehler der feinen Polyeder bzw. der feinen Darstellung linearer Dichten); bei mindestens
drei Punkten zusaetzlich Q_oo + c/m^2 + d/m^4, die Differenz beider Schaetzungen ist die angegebene Unsicherheit.

    python3 prototype/curved/analyze.py results/curved_sphere_glass.csv results/curved_sphere_gold.csv
(auch die JSON-Ausgabe von curved_galerkin.py)
"""
import csv
import json
import sys

import numpy as np

MIE = {"glass": 0.2150978, "gold": 0.5900182}


def extrapolate(ms, qs, m_min):
    ms, qs = np.asarray(ms, float), np.asarray(qs, float)
    sel = ms >= m_min
    ms, qs = ms[sel], qs[sel]
    if len(ms) < 2:
        return None, None
    A2 = np.column_stack([np.ones_like(ms), ms ** -2.0])
    q2 = np.linalg.lstsq(A2[-2:], qs[-2:], rcond=None)[0][0]           # aus den zwei feinsten Punkten
    if len(ms) >= 3:
        A4 = np.column_stack([np.ones_like(ms), ms ** -2.0, ms ** -4.0])
        q4 = np.linalg.lstsq(A4[-3:], qs[-3:], rcond=None)[0][0]
        return q4, abs(q4 - q2)
    return q2, None


def main(files):
    for f in files:
        if f.endswith(".csv"):
            rows = [dict(geometry=r["geometry"], space=r["space"], n=int(r["n"]), m=int(r["m"]), Q=float(r["Q_ext"]),
                         err=float(r["rel_error_vs_mie"])) for r in csv.DictReader(open(f))]
        else:
            rows = json.load(open(f))
        mat = "gold" if "gold" in f else "glass"
        q0 = MIE[mat]
        print(f"\n{mat}: Mie Q_ext = {q0}")
        groups = {}
        for r in rows:
            groups.setdefault((r.get("geometry", "curved"), r.get("space", "pc"), r["n"]), []).append(r)
        ns = sorted({k[2] for k in groups})
        cols = [("flat", "pc"), ("curved", "pc"), ("flat", "p1"), ("curved", "p1")]
        heads = ["eben/konst (heute)", "gekruemmt/konst", "eben/linear", "gekruemmt/linear"]
        if any(k[0] == "quadratic" for k in groups):                      # Stufe 1b
            cols.append(("quadratic", "p1")); heads.append("quadratisch/linear")
        print(f"{'N':>6} | " + " | ".join(f"{h:>22}" for h in heads))
        for n in ns:
            cells = []
            for g, s in cols:
                rs = sorted(groups.get((g, s, n), []), key=lambda r: r["m"])
                if g == "flat" and s == "pc":
                    rs = sorted(groups.get(("curved", "pc", n), []), key=lambda r: r["m"])[:1]   # m = 1: heutige Rechnung
                    cells.append(f"{rs[0]['err']:+.4%}" if rs else "-")
                    continue
                if not rs:
                    cells.append("-")
                    continue
                q, du = extrapolate([r["m"] for r in rs], [r["Q"] for r in rs], 2 if s == "p1" else 1)
                if q is None:
                    cells.append("-")
                    continue
                txt = f"{q / q0 - 1:+.4%}"
                txt += f" +- {du / q0:.1e}" if du is not None else " (2 Pkt.)"
                cells.append(txt + f" [m<={rs[-1]['m']}]")
            N = 20 * n * n
            print(f"{N:6d} | " + " | ".join(f"{c:>22}" for c in cells))
        print("  Rohwerte (Fehler gegen Mie je m):")
        for (g, s, n), rs in sorted(groups.items()):
            rs = sorted(rs, key=lambda r: r["m"])
            print(f"    {g:6s}/{s}  N = {20 * n * n:5d}: " + "  ".join(f"m={r['m']}: {r['err']:+.4%}" for r in rs))


if __name__ == "__main__":
    main(sys.argv[1:])
